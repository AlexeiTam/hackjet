"""
Similar to FourMomentumCompare, but uses a much more simple comparison method: directly comparing constituent indices
to test if same particles ended up in same jets.

"""

import uproot
import numpy as np
import matplotlib.pyplot as plt
import naive

#collecting final-state jets from naive, aeacus outputs

#assuming this script is being run right in the hackjet repository:
naive_jets = naive.JetReader("naive/naive_SIFT.root", "finalJets", auto_index=False)
aeacus_jets = naive.JetReader("SIFT-AEACuS/Cuts/aeacus_SIFT_withconstituents.root", "finalJets", auto_index=False)

#First key test: do we even end up with same number of jets
print(f"# of naive jets: {len(naive_jets)}")
print(f"# of aeacus jets: {len(aeacus_jets)}")


naive_idx = [jet.id for jet in naive_jets]
aeacus_idx = [jet.id for jet in aeacus_jets]


# print("naive:")
# print(naive_idx)

# print("aeacus:")
# print(aeacus_idx)

def match(idx1,idx2):
  if set(idx1) == set(idx2):
    return True
  else:
    return False


matched_idx = []

breakcount = 0
while len(aeacus_idx) > 0:
  oldlen = len(aeacus_idx)

  for id_n in naive_idx:
    for id_a in aeacus_idx:

      if match(id_n, id_a):
        n_can = id_n
        a_can = id_a

  matched_idx.append((n_can, a_can))
  naive_idx.remove(n_can)
  aeacus_idx.remove(a_can)

  if oldlen == len(aeacus_idx):
    print(f"while loop continuing: oldlen = {oldlen}; jets remaining = {len(aeacus_idx)}")
    breakcount += 1

  if breakcount > 10:
    print("no matches being made; breaking out")
    break


print("matching pairs:")
for i,matchingpair in enumerate(matched_idx):
  print(i, "...", matchingpair)





  

