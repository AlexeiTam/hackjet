from pathlib import Path
import re

import numpy as np


_MOMENTUM_LABEL = re.compile(r"^EP([0-3])_(\d+)$")


def read_cut_four_momenta(
    cut_filename,
    event_index=0,
    drop_undefined_jets=True,
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
                f"Jet {suffix} does not have all four component "
                f"labels. Missing: {missing_labels}"
            )

        raw_components = [
            label_to_value[label]
            for label in component_labels
        ]

        # A nonexistent jet slot is normally represented by four
        # UNDEF values.
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

    return (
        np.asarray(E_values, dtype=np.float64),
        np.asarray(px_values, dtype=np.float64),
        np.asarray(py_values, dtype=np.float64),
        np.asarray(pz_values, dtype=np.float64),
    )



E_array, px_array, py_array, pz_array = (
    read_cut_four_momenta("testing-output_003.cut")
)

four_momenta = np.column_stack(
    [E_array, px_array, py_array, pz_array]
)

print(four_momenta)

print(f"Total # of jets:{len(four_momenta)}")

import naive

outjets = [naive.MJet(E_array[i], px_array[i], py_array[i], pz_array[i], i) for i in range(len(E_array))]

naive.JetWriter("aeacus_SIFT.root", outjets, "finalJets")




