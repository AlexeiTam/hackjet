from pathlib import Path
import re

import numpy as np

# reading 4-momenta and constituents from EP0,EP1,EP2,EP3 and CIX columns in output .cut files of AEACuS SIFT implementation
_MOMENTUM_LABEL = re.compile(r"^EP([0-3])_(\d+)$")
_CONSTITUENT_LABEL = re.compile(r"^CIX_(\d+)$")


def read_cut_four_momenta(
    cut_filename,
    event_index=0,
    drop_undefined_jets=True,
    read_constituents=True,
    zero_based=True
):
    """
    Read jet four-momenta from an AEACuS .cut file.

    The component mapping is

        EP0_NNN -> E
        EP1_NNN -> px
        EP2_NNN -> py
        EP3_NNN -> pz

    Parameters
    ----------
    cut_filename : str or pathlib.Path
        Path to the AEACuS .cut file.

    event_index : int, optional
        Which event-data row to read, starting from zero.
        The default is the first event.

    drop_undefined_jets : bool, optional
        If True, omit jet slots for which all four components are
        UNDEF. If False, preserve those slots as np.nan.

    read_constituents: bool
        If True, will also read constituents
        defaults to True, but includes option to be False for back-compatibility with previous version that did not include it

    zero_based: bool
        If True, adjusts read constituent indices to be compatible with assignments by naive.JetReader from naive.py,
            which begins by assigning particles starting at index 0, compared to .lhco which starts at 1 (as first particle, with index 0, is nonphysical).

    Returns
    -------
    E_array : numpy.ndarray
        Jet energies.

    px_array : numpy.ndarray
        Jet x-momenta.

    py_array : numpy.ndarray
        Jet y-momenta.

    pz_array : numpy.ndarray
        Jet z-momenta.
    """

    cut_path = Path(cut_filename)
    lines = cut_path.read_text(encoding="utf-8").splitlines()

    # Locate the output-label row containing EP0_001, EP1_001, etc.
    header_index = None
    labels = None

    for line_index, line in enumerate(lines):
        candidate_labels = line.split()

        if any(
            _MOMENTUM_LABEL.fullmatch(label)
            for label in candidate_labels
        ):
            header_index = line_index
            labels = candidate_labels
            break

    if header_index is None or labels is None:
        raise ValueError(
            f"No EP0_NNN–EP3_NNN label row was found in "
            f"'{cut_path}'."
        )

    # Locate all data rows corresponding to this label row.
    data_rows = []

    for line in lines[header_index + 1:]:
        stripped_line = line.strip()

        if not stripped_line or stripped_line.startswith("#"):
            continue

        values = stripped_line.split()

        # A corresponding data row must have exactly one value
        # for every label.
        if len(values) == len(labels):
            data_rows.append(values)

    if not data_rows:
        raise ValueError(
            f"No data row corresponding to the output labels was "
            f"found in '{cut_path}'."
        )

    if not -len(data_rows) <= event_index < len(data_rows):
        raise IndexError(
            f"event_index={event_index} is invalid. The file "
            f"contains {len(data_rows)} event-data row(s)."
        )

    selected_values = data_rows[event_index]
    label_to_value = dict(zip(labels, selected_values))

    # Find every numbered jet slot appearing in the label row.
    jet_numbers = sorted(
        {
            int(match.group(2))
            for label in labels
            if (match := _MOMENTUM_LABEL.fullmatch(label))
        }
    )
    
    #CIX_jjj columns: {idx} of each jet's constituents
    #jet N stored in JET_001 has constituents in CIX(001+N)
    #1 in SIFT-IDX.dat card
    constituents = []

    if read_constituents:
        cix_numbers = sorted(
            int(match.group(1))
            for label in labels
            if (match := _CONSTITUENT_LABEL.fullmatch(label))
        )

        if not cix_numbers:
            raise ValueError(
                f"No CIX_NNN label found in {cut_path}"
            )

        shelf_number = cix_numbers[0] - 1
        index_offset = 1 if zero_based else 0
        

    E_values = []
    px_values = []
    py_values = []
    pz_values = []

    for jet_number in jet_numbers:
        suffix = f"{jet_number:03d}"

        component_labels = [
            f"EP0_{suffix}",
            f"EP1_{suffix}",
            f"EP2_{suffix}",
            f"EP3_{suffix}",
        ]

        missing_labels = [
            label
            for label in component_labels
            if label not in label_to_value
        ]

        if missing_labels:
            raise ValueError(
                f"Jet {suffix} does not have all four-momentum labels. Missing: {missing_labels}"
            )

        raw_components = [
            label_to_value[label]
            for label in component_labels
        ]

        #an empty jet slot normally gets just four UNDEFs (from having to pre-set shelf size without knowing ahead of time how many final-state jets there would be
        if (
            drop_undefined_jets
            and all(value.upper() == "UNDEF"
                    for value in raw_components)
        ):
            continue

        components = [
            np.nan if value.upper() == "UNDEF" else float(value)
            for value in raw_components
        ]

        E, px, py, pz = components

        E_values.append(E)
        px_values.append(px)
        py_values.append(py)
        pz_values.append(pz)

    
        if read_constituents:
            slot = f"{shelf_number + jet_number:03d}"    #SIFT output jet index 
            raw_ids = label_to_value.get(f"CIX_{slot}")

            if raw_ids is None:
                raise ValueError(
                    f"Jet {suffix} has no CIX_{slot} column; widen the "
                    f"card's shelf range to match the OUT ranges."
                )

            ids = (
                []
                if raw_ids.upper() == "UNDEF"
                else [int(value) - index_offset
                      for value in raw_ids.split(",")]
            )

            #cross check with AEACuS native constituent count
            raw_count = label_to_value.get(f"JET_{slot}")

            if raw_count is not None and int(raw_count) != len(ids):
                raise ValueError(
                    f"Jet {suffix}: CIX_{slot} lists {len(ids)} constituents but JET_{slot} = {raw_count}."
                )

            constituents.append(ids)

    four_momenta = (
        np.asarray(E_values, dtype=np.float64),
        np.asarray(px_values, dtype=np.float64),
        np.asarray(py_values, dtype=np.float64),
        np.asarray(pz_values, dtype=np.float64),
    )

    if read_constituents:
        # Cross-check the jet count against the shelf's own count column.
        raw_njets = label_to_value.get(f"JET_{shelf_number:03d}")

        if (
            drop_undefined_jets
            and raw_njets is not None
            and int(raw_njets) != len(E_values)
        ):
            raise ValueError(
                f"JET_{shelf_number:03d} = {raw_njets} jets, but {len(E_values)} were read; widen the card's OUT ranges."
            )

        return (*four_momenta, constituents)

    return four_momenta



E_array, px_array, py_array, pz_array, constituents = (
    read_cut_four_momenta("testing-output_006.cut", read_constituents=True)
)    #collecting read 4-momenta components and constituent lists for writing to .root

four_momenta = np.column_stack(
    [E_array, px_array, py_array, pz_array]
)

print(four_momenta)

print(f"Total # of jets:{len(four_momenta)}")

import naive

outjets = [naive.MJet(E_array[i], px_array[i], py_array[i], pz_array[i], constituents[i]) for i in range(len(E_array))]    
#collection of read final-state jets

naive.JetWriter("aeacus_SIFT_withconstituents.root", outjets, "finalJets")




