"""
A python mini-module to take .root files with one TTree (named "pdEventTree" by default"),
reads the four-momenta, in TBranches "E", "px", "py", "pz" for each input particle in the event, 
then uses naive module to construct MJets out of four-momenta, 
and extracts eta, phi, pt, and jet mass, then produces .lhco file with this data,
for input to SIFT clustering in aeacus.pl.

Output .lhco file includes "#", "TYP", "ETA", "PHI", "PT", "JMAS", "NTRK", "BTAG", "HAD/EM", and "DUM1" and "DUM2"


"#" = particle index for tracking input particles throughout clustering
"TYP" = still figuring it out; maybe it's PDG IDs?
"ETA", "PHI", "PT", "JMAS" = aforementioned extracted values from 4-momenta in .root file
"NTRK", "BTAG", "HAD/EM" = detector information; from expected Delphes event record
"DUM1" and "DUM2" = dummy columns 

"""

import naive
import uproot

input_filename = "PythiaEventsBatchTest.root"

output_filename = "input_event.lhco"

tree_name = "pdEventTree"

"""
with uproot.open(input_filename) as root_file:
  tree = root_file[tree_name]

  arrays = tree.arrays(
    ["E", "px", "py", "pz"],
    library = "np",
  )

arraystotal = (arrays["E"])*(arrays["px"])
print(len(arrays["E"]))
print(len(arraystotal))

print(arrays["E"])
print(arrays["px"])
print(arraystotal)

# Past attempt to use arrays(), useful but in this case another easier method is available
"""

def root_to_lhco_(input_root, output_lhco, treename):
  """
  documentation coming soon...
  """
  inset = naive.JetReader(input_root, treename)

  eta_array = [jet.eta for jet in inset]
  phi_array = [jet.phi for jet in inset]
  pt_array = [jet.pT for jet in inset]
  jmas_array = [jet.m for jet in inset]

  with open(output_lhco, "w") as output:
    #column-name header
    output.write(
        "#  TYP    ETA    PHI       PT     JMAS"
        "  NTRK  BTAG  HAD/EM  DUM1  DUM2\n\n"
    )

    # This file contains one event, numbered 1.
    output.write("0 1 0\n")

    numparticles = len(inset)
  
    for i in range(numparticles):
      typ = 4
      eta = eta_array[i]
      phi = phi_array[i]
      pt = pt_array[i]
      jmas = jmas_array[i]

      output.write(
            f"{i + 1:d} "
            f"{typ:d} "
            f"{eta:+.3f} "
            f"{phi:+.3f} "
            f"{pt:.3f} "
            f"{jmas:.3f} "
            f"0.0 "       # NTRK
            f"0.0 "       # BTAG
            f"0.000 "     # HAD/EM
            f"0.0 "       # DUM1
            f"0.0\n"      # DUM2
        )

    output.write(
            f"{len(inset) + 1:d} "
            f"6 "
            f"0.0 "
            f"0.0 "
            f"0.0 "
            f"0.0 "
            f"0.0 "       # NTRK
            f"0.0 "       # BTAG
            f"0.000 "     # HAD/EM
            f"0.0 "       # DUM1
            f"0.0\n"      # DUM2
        )
  print(f"Wrote {len(inset)} particles to lhco.")
  



  
#doing conversion

root_to_lhco_("PythiaEventsBatchTest.root", "testing-output.lhco", "pdEventTree")





