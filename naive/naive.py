#Dependencies

import numpy as np

import uproot

import pandas as pd


#important phi-related functions,
# to wrap phi in the [0,2pi) range and take differences of phi accounting for 2*pi-periodicity

pi = np.pi
two_pi = 2*pi

def wrap_phi(phi):
  """
  Wrap phi azimuthal angle to the range [0, 2*pi).

  Parameters:
      phi (float): The input angle in radians.

  Returns:
      float: The wrapped angle in radians.
  """
  return phi % (2 * pi)

def delta_phi(phi1,phi2):
  """
  Calculate the difference between two azimuthal angles (phi1, phi2) in radians,
  accounting for the 2*pi periodicity.

  Parameters:
      phi1 (float): The first angle in radians.
      phi2 (float): The second angle in radians.

  Returns:
      float: The difference between the two angles in radians.
  """
  phi1 = wrap_phi(phi1)
  phi2 = wrap_phi(phi2)
  dphi = phi1 - phi2
  conditionlist = [dphi < -1.0*np.pi, dphi > np.pi]
  choicelist = [np.abs(dphi + 2*np.pi), np.abs(dphi - 2*np.pi)]
  dphi = np.select(conditionlist, choicelist, default=np.abs(dphi))
  return dphi


#custom class for jets:
class MJet:
  """
  class for jets, containing 4-momenta and other useful information for jet clustering;
  argument is simply 4-momenta (E,px,py,pz)
  makes use of NumPy for vectorized operations, and other math operations

  Attributes:
      E (float): Energy of the jet.
      px (float): x-component of the jet's momentum.
      py (float): y-component of the jet's momentum.
      pz (float): z-component of the jet's momentum.
      id (list(int)): Jet ID number. meant to track input particles and clustering history.

      pT (float): Transverse momentum of the jet.
      p (float): Momentum of the jet.
      phi (float): Azimuthal angle of the jet.
      eta (float): Pseudorapidity of the jet.
      m (float): Mass of the jet.

      ghost_list (list(float)): List of ghost jets contained in this MJet, in form of y, phi coordinates where the original ghost jet was formed.
  """

  def __init__(self, E, px, py, pz, id_num,ghost_list=[]):
    self.E = E
    self.px = px
    self.py = py
    self.pz = pz
    self.id = id_num



    self.pT = float(np.sqrt(px**2 + py**2))
    self.p = float(np.sqrt(px**2 + py**2 + pz**2))

    self.phi = wrap_phi(float(np.arctan2(self.py, self.px)))

    # Calculate mass, add tolerance (1e-6)for floating point issues
    m_squared = self.E**2 - self.p**2
    if m_squared < 0 and np.isclose(m_squared, 0, atol=1e-6): # Treat small negative as zero
        self.m = 0.0
    else:
        self.m = float(np.sqrt(m_squared))


    # Calculate rapidity y, pseudorapidity eta, handling edge cases for rapidity
    # maybe also use np.isclose()

    #FastJet makes use of MaxRap = 1e5 to stand in for infinity in edge cases
      #may use MaxRap += pz, to adjust slightly and lift degeneracy in extra-low pz cases

    #pseudorapidity, eta
    if self.pT == 0:
      self.eta = 1e5
    elif self.pz == 0:
      self.eta = 0
    else:
      self.theta = float(np.arctan(self.pT / self.pz))
      if self.theta <= 0:
        self.theta += float(np.pi)
      self.eta = float(-np.log(np.tan(self.theta/2)))

    #rapidity, y
    if self.E + self.pz == 0:
        self.y = -1e5
    elif self.E - self.pz == 0:
        self.y = 1e5
    else:
        self.y = (1/2)*float( np.log( (self.E + self.pz) / (self.E - self.pz) ) )


    #ghost_list: a set of y, phi coordinates for ghost jets contained in this MJet
    self.ghost_list = ghost_list



  def __str__(self):

    output = {
        "id": self.id,
        "E": self.E,
        "px": self.px,
        "py": self.py,
        "pz": self.pz,
        "pT": self.pT,
        "p": self.p,
        "ghost_list": self.ghost_list,
        "phi": self.phi,
        "eta": self.eta,
        "m": self.m
    }
    return str(output)


  def __add__(self, other):
    return MJet(self.E + other.E, self.px + other.px, self.py + other.py, self.pz + other.pz, self.id + other.id, self.ghost_list + other.ghost_list)

  def __eq__(self, other):
    return (self.E == other.E) and (self.px == other.px) and (self.py == other.py) and (self.pz == other.pz) and set(self.id) == set(other.id)


#clustering distances
# in general,

# pairwise:

#   d_ij = min(pT_i^(2*m), pT_j^(2*m));

# beam:

#   d_iB = min(pT_i^(2*m))

#need to set a mode m:
# m = 1: kT-algorithm
# m = 0: Cambridge-Aachen
# m = -1: anti-kT algorithm

def d_ij(mode, j1, j2, R_0, epsilon=1e-6):
  """
  pairwise distance between two jets, using generalized kT algorithm:
  min(pT_i^(2*m), pT_j^(2*m));

  Parameters:
      mode (int): The mode of the distance calculation:
          - 1: kT-algorithm
          - 0: Cambridge-Aachen
          - (-1): anti-kT algorithm
      j1 (MJet): The first jet.
      j2 (MJet): The second jet.
      R_0 (float): The jet radius parameter.
      epsilon (float): A small positive value to avoid division by zero from floating point inaccuracies

  Returns:
      float: The distance between the jets.
  """

  pT1_sq = (j1.pT)**(2*mode)
  pT2_sq = (j2.pT)**(2*mode)

  min_pT_sq = min( pT1_sq , pT2_sq )

  delta_R_sq = (j1.y - j2.y)**2 + (delta_phi(j1.phi,j2.phi))**2

  output = min_pT_sq * (delta_R_sq / (R_0**2))

  return output

def d_iB(mode, j, epsilon=1e-6):
  """
  beam distance. Just pT^(2*mode)

  Parameters:
      j (MJet): The jet.

  Returns:
      float: The beam distance.
  """
  # The beam distance for kT is kT^2
  if mode == -1:
    return (j.pT + epsilon)**(2*mode)
  else:
    return j.pT**(2*mode)

def gen_clustering(mode, j_initial, R_0, thresholding=False, pT_threshold=0.0, eta_threshold=0.0):
  """
  Function for clustering an input set of jets, using anti-kT clustering algorithm,
  to return a set of final-state jets.

  Parameters:
    mode: The clustering mode to use.
      1: kT-algorithm
      0: Cambridge-Aachen
      (-1): anti-kT algorithm

    j_initial: List of initial jets,
    R_0: the jet radius parameter,

    thresholding: Boolean, configurable by user in case they wish to apply jet pT and eta thresholds; off by default
    pT_threshold: minimum jet pT,
    eta_threshold: minimum jet eta,

  Returns: clustering_history, final_jets

    clustering_history: list of tuples--('final'/'merge', ids of jets being promoted)
    final_jets: List of final-state jets

  """

  active = list(j_initial)  #list of jets being clustered, make a copy to avoid modifying original
  final_jets = []  #to store final state jets

  clustering_history = []

  # Main loop
  while active:

    d_min = float('inf')
    indices_min = None  #indices of jets to either be clustered or promoted to final state
    type_min = None #determine whether d_iB, d_ij is minimum, and thus whether to cluster/promote jet

    #Calculate d_iB, and see if anything can be promoted to final jet

    for i, j_i in enumerate(active):

      diB = d_iB(mode, j_i)

      if diB < d_min:
        d_min = diB
        indices_min = (i,)
        type_min = 'beam'

    #Calculate d_ij and see if anything can be clustered

    for i in range(len(active)):

      for j in range(i+1, len(active)):

        dij = d_ij(mode, active[i], active[j],R_0)

        if dij < d_min:
          d_min = dij
          indices_min = (i, j)
          type_min = 'pairwise'

    #Now, either cluster/promote

    if type_min == 'beam':
      j_candidate = active.pop(indices_min[0])
      # final_jets.append(j_candidate)
      # clustering_history.append(('final',j_candidate.id))   #OLD VERSION: reverses the ordering of the fastjet implementation
      final_jets = [j_candidate] + final_jets
      clustering_history = [('final',j_candidate.id)] + clustering_history

    elif type_min == 'pairwise':
      idx1, idx2 = sorted(indices_min, reverse=True)
      j_i = active.pop(idx1)
      j_j = active.pop(idx2)

      active.append(j_i + j_j)
      clustering_history.append(('merge', (j_i.id, j_j.id)))

    else:
      raise ValueError("Invalid type_min value: No valid clustering or beam distance found.")

  # Apply thresholding after clustering is complete
  if thresholding:
    filtered_final_jets = []
    for jet in final_jets:
      # Assuming eta_threshold is a maximum absolute eta value
      if jet.pT >= pT_threshold and abs(jet.eta) <= eta_threshold:
        filtered_final_jets.append(jet)
    filtered_final_jets = sorted(filtered_final_jets, key=lambda x: x.pT, reverse=True) #final sorting in descending pt
    return filtered_final_jets
  else:
    final_jets = sorted(final_jets, key=lambda x: x.pT, reverse=True) #final sorting in descending pt
    return clustering_history, final_jets

############################################################
#################### Ghost Jets and Jet Area ####################

#still preliminary: let's use two versions, swappable by version parameter
def GhostJetPopulator(jet_input, y_min, y_max, ny, phi_min, phi_max, nphi, epsil=1e-6, version=1):
  """
  Populates a set of input jets with ghosts jets. In particular, takes jet_input set,
  and adds a series of ny*nphi many ghost jets in a uniform (y,phi) grid.

  Parameters:
      jet_input (list): List of input jets.
      y_min, y_max (float): rapidity (y) bounds of the uniform (y,phi) grid where ghost jets are populated.
      phi_min, phi_max (float): phi-bounds of the uniform (y,phi) grid where ghost jets are populated.
      ny, nphi: number of ghost jets populated along y,phi axis respectively, resulting in N = ny*nphi many ghost jets.

      epsil (float): Energy of ghost jet, by default a soft, massless jet, with E = |p|= epsil.

      version (-1 or 1): alternating between versions of ghost jet indexing.

  Returns:
      list: List of input jets with added ghost jets.

  """

  y_axis = np.linspace(y_min, y_max, ny)
  phi_axis = np.linspace(phi_min, phi_max, nphi)


  y_grid, phi_grid = np.meshgrid(y_axis, phi_axis)

  y_phi_coordinates = np.column_stack((y_grid.ravel(), phi_grid.ravel()))

  # creating soft jets, using E = |p|= epsil,
  # and obtaining p-components through y, phi with
  # pz = E*tanh(y); pT = E**2 - pz**2,
  # px = pT*cos(phi); py = pT*sin(phi)

  pz_values = epsil * np.tanh(y_grid.ravel())
  pT_values = np.sqrt(epsil**2 - pz_values**2)
  px_values = pT_values * np.cos(phi_grid.ravel())
  py_values = pT_values * np.sin(phi_grid.ravel())

  #assigns indices to ghost jets sequentially, as len(jet_input), len(jet_input)+1, len(jet_input)+2, and so on.
  #input particles in jet_input assumed to be indexed as integers, from 0, 1, 2, ..., len(jet_input)-1.

  if version == 1:
    ghost_jets = [MJet(epsil, px, py, pz, [i+len(jet_input)]) for i, (px, py, pz) in enumerate(zip(px_values.ravel(), py_values.ravel(), pz_values.ravel()))]
  elif version == -1:
    ghost_jets = [MJet(epsil, px, py, pz, [-(i+1)]) for i, (px, py, pz) in enumerate(zip(px_values.ravel(), py_values.ravel(), pz_values.ravel()))]
  else:
    raise ValueError("Invalid version value. Must be -1 or 1.")
  #alternate possibility: index as negative numbers [i+len(jet_input)] --> [-(i+1)] (as i starts at 0), so one can more easily distinguish ghost jets by index


  jet_input += ghost_jets

  return jet_input





############################################################
#################### SIFT ####################
# makes use of a very different distanc metric

def d_ij_SIFT(j1, j2):
  """
  SIFT (Scale-Invariant Filtered Tree) distance metric between two jets,
  (2* p^{\mu}_i * p_{\mu}_j) / ((E_T^i)^2 + (E_T^j)^2)
  with E_T^2 = E^2 - p_z^2.

  Parameters:
      j1 (MJet): The first jet.
      j2 (MJet): The second jet.

  Returns:
      float: The distance between the jets.
  """

  dm = ( (j1.E)*(j2.E) ) - ( (j1.px)*(j2.px) ) - ( (j1.py)*(j2.py) ) - ( (j1.pz)*(j2.pz))

  E_T_i_squared = j1.E**2 - j1.pz**2
  E_T_j_squared = j2.E**2 - j2.pz**2

  denom = E_T_i_squared + E_T_j_squared

  output = dm/denom


  return output


#defining SIFT clustering function
#based heavily on general clustering function above

# future: incorporating clustering, drop, isolate conditions
# otherwise, we begin with just a version that brings everything to one big jet

def SIFT(j_initial,thresholding=False, pT_threshold=0.0, eta_threshold=0.0):

  active = list(j_initial)
  final_jets = []

  clustering_history = []

  #main loop
  while active:

    #if only one jet remaining:

    if len(active) == 1:
      final_jets.append(active.pop(0))
      break

    #initializing minimum distance variables
    d_min = float('inf')
    indices_min = None

    for i, j_i in enumerate(active):

      for j in range(i+1, len(active)):

        dij = d_ij_SIFT(j_i, active[j])

        if dij < d_min:
          d_min = dij
          indices_min = (i, j)

    if indices_min is not None:
      idx1, idx2 = sorted(indices_min, reverse=True)
      j_i = active.pop(idx1)
      j_j = active.pop(idx2)

      active.append(j_i + j_j)
      clustering_history.append(('merge', (j_i.id, j_j.id)))
    else:
      raise ValueError("No valid clustering found.")

  if thresholding:
    filtered_final_jets = []
    for jet in final_jets:
      if jet.pT >= pT_threshold and abs(jet.eta) <= eta_threshold:
        filtered_final_jets.append(jet)

    filtered_final_jets = sorted(filtered_final_jets, key=lambda x: x.pT, reverse=True)
    return clustering_history, filtered_final_jets
  else:
    final_jets = sorted(final_jets, key=lambda x: x.pT, reverse=True)
    return clustering_history, final_jets


############################################################
#################### I/O ####################

#Jet reader: use uproot to make output files and read them
# functions for reading and writing: ROOT files <-> naive.MJets or other objects made from 4-momenta


def JetReader(f1, treename):
  """
  takes a .root file input of PYTHIA particles,
  returns a list of naive.MJet objects,
  meant to then be used in clustering.

  inputs:
  f1 = name of .root file;
    must contain E, px, py, pz columns for MJet inputs.
  treename = name of tree in .root file.

  returns a list of naive.MJet objects.
  """

  infile = uproot.open(f1)
  tree = infile[treename]

  E_array = tree["E"].array().to_numpy()
  px_array = tree["px"].array().to_numpy()
  py_array = tree["py"].array().to_numpy()
  pz_array = tree["pz"].array().to_numpy()


  return [MJet(E_array[i], px_array[i], py_array[i], pz_array[i], [i]) for i in range(len(E_array))]


def JetWriter(fname, outjets, outtree):
  """
  Takes an input list of MJet objects, or objects with 4-momentum attributes
  (e.g., for jets j, has j.E, j.px, j.py, j.pz)
  and writes them to a .root file.

  inputs:
  fname - name for output file
  outjets - list of clustered jets
  outtree - name for tree in output file
  """

  outfile = uproot.recreate(fname)
  output_dict = {
    "E": [jet.E for jet in outjets],
    "px": [jet.px for jet in outjets],
    "py": [jet.py for jet in outjets],
    "pz": [jet.pz for jet in outjets],
    "id": [jet.id for jet in outjets],
    "pT": [jet.pT for jet in outjets],
    "eta": [jet.eta for jet in outjets],
    "phi": [jet.phi for jet in outjets]
    }

  output_df = pd.DataFrame(output_dict)

  outfile.mktree(outtree, output_df)

  outfile[outtree]