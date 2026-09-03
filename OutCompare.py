"""
OutCompare: to handle comparisons of .root outputs between naive and SIFT (from aeacus.pl);
each containing the 4-momenta of final-state clustered jets

current mechanism:

--first test: same number of final state jets?

--second test:
  **note to self: treat aeacus outputs as TRUTH**
  *separate into two 'active' pools: active_aeacus and active_naive
  *for s_i jets in active_aeacus, compare to all naive jets n_j in active_naive;
  *determine a "gap" g_ij quantity between s_i and n_j, equal to dE**2 + dpx**2 + dpy**2 + dpz**2
  *for min g_ij, 'pair' the two,
    place tuple pair (s_i, n_j) into a 'matched' list: matched_jets, remove s_i and n_j from 'active' pools
  *keep going until 'active' pools are empty

--third test:
  *for each pair in matched_jets, determine percent difference between each component of each 4-momenta,
    [((E1-E2)/(E1))*100, ((px1-px2)/(px1))*100, ((py1-py2)/(py1))*100, ((pz1-pz2)/(pz1))*100],
    with (E1, px1, py1, pz1) = (s_i.E(), s_i.px(), s_i.py(), s_i.pz()) and (E2, px2, py2, pz2) = (n_j.E(), n_j.px(), n_j.py(), n_j.pz())

--final test: histogram of percent differences for each component
"""



import uproot
import numpy as np
import matplotlib.pyplot as plt
import naive

#collecting final-state jets from naive, aeacus outputs

#assuming this script is being run right in the hackjet repository:
naive_jets = naive.JetReader("naive/naive_SIFT.root", "finalJets")
aeacus_jets = naive.JetReader("SIFT-AEACuS/Cuts/aeacus_SIFT.root", "finalJets")

print(f"# of naive jets: {len(naive_jets)}")
print(f"# of aeacus jets: {len(aeacus_jets)}")
def gap(j1, j2):
  """
  returns the <gap> mentioned in the top of the script, for the purposes of trying to find matches for each jet in the <active> pools
  """

  dE = j1.E - j2.E
  dpx = j1.px - j2.px
  dpy = j1.py - j2.py
  dpz = j1.pz - j2.pz

  output_gap = dE**2 + dpx**2 + dpy**2 + dpz**2

  return output_gap


#iterate through s_i, n_j
matched_jets = []
while len(aeacus_jets) > 0:
  min_gap = float(np.inf)
  for s_i in naive_jets:
    for n_j in aeacus_jets:
      gap_found = gap(s_i, n_j)
      if gap_found < min_gap:
        s_can = s_i
        n_can = n_j
  matched_jets.append((s_can, n_can))
  naive_jets.remove(s_can)
  aeacus_jets.remove(n_can)
  
  #remove the matches:
#  print(f"Match found: {str(s_i)} & {str(n_j)}")
#  matched_jets.append((s_can, n_can))
    
#  naive_jets.remove(s_can)
#  aeacus_jets.remove(n_can)


def percent_gap(j_1, j_2):

  dE = j_1.E - j_2.E
  dpx = j_1.px - j_2.px
  dpy = j_1.py - j_2.py
  dpz = j_1.pz - j_2.pz

  dE_percent = (dE/(j_1.E))*100
  dpx_percent = (dpx/(j_1.px))*100
  dpy_percent = (dpy/(j_1.py))*100
  dpz_percent = (dpz/(j_1.pz))*100

  return [float(dE_percent), float(dpx_percent), float(dpy_percent), float(dpz_percent)]


print("investigating matches...")

percent_gap_matrix = []
for j1, j2 in matched_jets:
#  print(f"Jets: {str(j1)} & {str(j2)}")
  print(percent_gap(j1, j2))
  percent_gap_matrix.append(percent_gap(j1, j2))

E_percentgap = [p[0] for p in percent_gap_matrix]
px_percentgap = [p[1] for p in percent_gap_matrix]
py_percentgap = [p[2] for p in percent_gap_matrix]
pz_percentgap = [p[3] for p in percent_gap_matrix]

plt.hist(E_percentgap, bins=10)
plt.title("E Percent Gap")
plt.xlabel("Percent Difference")
plt.ylabel("Frequency")

plt.tight_layout()
plt.savefig("E_percent_gap.png", format='png')
plt.close()

plt.hist(px_percentgap, bins=10)
plt.title("$P_x$ Percent Gap")
plt.xlabel("Percent Difference")
plt.ylabel("Frequency")

plt.tight_layout()
plt.savefig("px_percent_gap.png", format='png')
plt.close()

plt.hist(py_percentgap, bins=10)
plt.title("$P_y$ Percent Gap")
plt.xlabel("Percent Difference")
plt.ylabel("Frequency")

plt.tight_layout()
plt.savefig("py_percent_gap.png", format='png')
plt.close()

plt.hist(pz_percentgap, bins=10)
plt.title("$P_z$ Percent Gap")
plt.xlabel("Percent Difference")
plt.ylabel("Frequency")

plt.tight_layout()
plt.savefig("pz_percent_gap.png", format='png')
plt.close()


print("basic comparison:")


for j1, j2 in matched_jets:

  print(f"E: {j1.E:.2f} , {j2.E:.2f} -- px: {j1.px:.2f} , {j2.px:.2f} -- py: {j1.py:.2f} , {j2.py:.2f} -- pz: {j1.pz:.2f} , {j2.pz:.2f}")











