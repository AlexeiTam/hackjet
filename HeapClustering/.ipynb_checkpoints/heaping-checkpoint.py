"""
heaping.py

a custom module for implementing heap into jet clustering,
  particularly handles storing pairwise distances between clusters and comparing distances

new class: DistHeap
  -- takes matrix of distances D[i][j], built from pairs j_i and j_j from <active_jets> list
  -- 1-indexed; runs from H[s], H[s+1], ..., H[e]
  -- needs to satsify H[i] <= H[2i]; H[i] <= H[2i+1] for i = s, ..., e//2

  Pair: 
        Pair[dh_i] = tuple of indices a_i and a_j, 
                      indicating indices of two jets in <active_jets>, to which 
                      the distance in DistHeap slot <dh_i> corresponds

  Where:
        Where[a_i, a_j] - list of indices b, such that d_ij(j_i,j_j) = H[b]



shift_down: x,p,i

  wants to insert x distance corresponding to pair p, into heap H slot i (at the root H[1] by default)
  need to H heap above slot s to be valid!
  keep evaluating prospective "children" values (H[2*i], H[2*i+1])and compare to x: 
      if heap conditions satisfied ==> can just break and place x at slot i
      if not, move to next generation (i,j = j, 2*j)
        --am I convinced that this is the most efficient way to get to a good heap slot for x...?

shift_up: x,p.i

  same as shift_down, but in reverse;
  shifting "up" the heap instead of "down"; 
  need to compare x against prospective "parent"; 
    if heap conditions satisfied ==> can just break and place x at slot i
    if not, move up a generation (i,j = j, j//2)



"""

#! important dependency: needs naive unless we plan on just importing the fuctions here...

import naive
import time 

def jet_key(a,b):
  """
  Will use sets as identifiers for jets and pairs of jets.  Using frozensets as keys for <Where> dict;
  immutability of frozenset allows for use as key in dicts
  """
  return frozenset([a.id + b.id])

  
class DistHeap:

  def __init__(self, clustering_objects, pairwise_distance):
    """
    A heap-based alternative to the matrix-based method of earlier iterations in <naive.py>,
      based on Kurita paper, proposing an efficient, agglomerative clustering algorithm with a heap
      <https://www.sciencedirect.com/science/article/pii/003132039190062A>
    
    Accepts:

    clustering_objects (list/array) = list or array of objects that are to be clustered;

    pairwise_distance (function) = a function f(j_a, j_b) where j_a, j_b are in <clusering_objects>, 
                                    and returns the "pairwise distance" between clusters j_a, j_b
    
    DistHeap holds two key arrays:

      H:  -- the "heap" itself
          -- 1-indexed; from H[1], ..., H[e] (going to stick a [None] at the beginning, that shouldn't be considered)
          -- i.e, H[i] = pairwise_distance(j_a, j_b) = d_ij (j_a,j_b) where j_a, j_b are jets (naive.MJet)
          -- must satisfy:
              * H[i] <= H[2i]; H[i] <= H[2i+1] for i = 1,...,e//2

      Pairs:  -- the "Pairs" list;
              -- Pairs[i] identifies the pair of jets whose pairwise distance is in H slot [i]
              -- i.e., Pairs[i] = [j_a, j_b] (same j_a, j_b of H[i] above)

    In addition:
      
      Where: 
              -- "Pairs" and "Where" are to serve as reverse maps of each other;
              -- "Where" will return the heap slot index i, to which pair/jet p corresponds
                  -- if trying to locate jet, needs to check if set(jet.id) is a subset of corresponding Where element,
                  as opposed to if the two sets match
              -- i.e., if Pairs[i] = p, then Where(p) returns i.
    """

    distlists = [] 
    pairlists = []

    for i in range(len(clustering_objects)):
      for j in range(i+1, len(clustering_objects)):
        pairlists.append([clustering_objects[i], clustering_objects[j]])
        distlists.append(pairwise_distance(clustering_objects[i], clustering_objects[j]))

    self.e = len(distlists)
    self.H = [None] + distlists
    self.Pairs = [None] + pairlists
    self.Where = {jet_key(*self.Pairs[i]: i for i in range(1, len(distlists) + 1))}
    self.active = {frozenset(j.id): j for j in clustering_objects}#a full list of all of the 'active' jets still being clustered; equivalent to similarly named list from naive.SIFT

    s = (len(distlists) ) // 2 + 1
    while (s > 1):
      s -= 1
      self.shiftdown(s,len(distlists)) #start evaluating from the lowest end of the tree back upwards
    print(f"heap constructed, with {len(clustering_objects)} inputs, {len(self.H)} many heap entries, {len(self.Pairs)} many pairs")

  def where_locator(self, p):
    """
    serves role of "Where" from Kurita's heap paper;

    input: p -- a list of two jets/clusters j_a, j_b = p[0], p[1]

    must take in <p>, determine the heap slot <i> its pairwise distance is associated
    i.e., if Pair[i] = p; then Where[p] = i

    will identify where p = [j_a, j_b] are by matching naive.Mjet.ids (*assuming j_a, j_b are naive.MJet objects)
    """
    j_a = p[0]
    j_b = p[1]
    pair_identifier = set(j_a.id + j_b.id)

    pairs_id_listing = [set(j_1.id+j_2.id) for [j_1,j_2] in self.Pairs[1:]]

    output = pairs_id_listing.index(pair_identifier)
    return output

  def jet_locator(self, jet, setting):
    """
    the same as where_locator, but instead returns a list of indices where <jet> is located in heap
    i.e., if jet = j_a, then returns indices = [i where H[i] = d_ij(j_a, ...)]

    will return different things based on <setting> input:

    setting "all":
      designed for "remove" action; will just return <heap_slot>,
      a list of indices from self.H heap where self.Pairs and self.H elements need to be removed;
      can be done by just iterating backwards from list.

    setting "single":
      designed for "update" action; will just return <>
    """

    if setting == "single":
      for i in range(1, len(self.Pairs)):
        a, b = self.Pairs[i]
        if a is jet: return i, b
        if b is jet: return i, a
      return None, None

    jet_ID = set(jet.id)
    pairs_id_listing = [set(j_1.id + j_2.id) for [j_1,j_2] in self.Pairs[1:]]

    # output = [i for i,x in enumerate(pairs_id_listing) if jet_ID.issubset(x)]
    # return output

    if setting == "all":
      #print("searching for all indices of given jet...")
      return [i for i,x in enumerate(pairs_id_listing, start=1) if jet_ID.issubset(x)]

    #old, defunct method for "single" setting
    # if setting == "single":
    #   #print("will return heap_slot, other_jet: the slot where the first instance of this jet occurs in heap, and the other jet in the pairing where the input jet first occurs")
    #   #is the jet even still there?
    #     #if not, return None

    #   if len([i for i,x in enumerate(pairs_id_listing, start=1) if jet_ID.issubset(x)]) == 0:
    #     return None, None

      
    #   #jet is still there, let's keep going
    #   for i,x in enumerate(pairs_id_listing, start=1):

    #     if jet_ID.issubset(x):
    #       heap_slot = i #the heap slot where jet's first associated distance, first pairing in self.H and self.Pairs

    #       for jet_cand in self.Pairs[heap_slot]:
    #         if jet_cand.id != jet.id:
    #           return heap_slot, jet_cand

          
        

    
  
  def place_in_i(self, i, x, p):
    """
    inputs: i (int), x (float), p ( list( list, list ) )
    places item <x> from pair <p> into heap slot <i>
    """
    self.H[i] = x
    self.Pairs[i] = p
    self.Where[jet_key(*p)] = i

  def shiftdown(self, s, e):
    """
    works from root H[1] down the generations;
    evaluates whether <x> = H[s] belongs at slot s, 
      by comparing to its "children" H[2s] or H[2s+1] (whichever is lower)
      
      if heap conditions violated (by H[s] being bigger than its children), 
        need to go further down a generation
      once heap conditions met, <x> allowed to stay in slot s


    inputs:
      s - starting point of heap H, within which the starting heap H is valid up to e
      e - ending point of heap H, within which the starting heap H is valid up to s

    action: <self.H> becomes a fully valid heap
    """

    H = self.H
    Pairs = self.Pairs
    
    x = H[s]
    p = Pairs[s]

    i = s
    j = 2*s

    while (j <= e):
      #will keep evaluating validity of heap slot <i> for x, from s -- down the generations -- all the way to end slot <e>
      if (j < e) and H[j] > H[j+1]:
        j += 1  # will compare x to H[j]; need to compare against smallest child to ensure both heap conditions satisfied
        #j = e: we're at the end; can't make swap between possible children to compare to <x>
      
      #condition met: can just break out
      if x <= H[j]:
        break

      
      #not break -- x needs to move down, so child takes its place as it is smaller
      #child moves up to take up hole left once swap is made
      self.place_in_i(i, H[j], Pairs[j])

      #move down a generation:
      i = j
      j = 2*j

    self.place_in_i(i, x, p) # break: can place x in slot i while satisfying heap conditions



  def shiftup(self, s):
    """
    same as shiftdown, but in reverse;
    try and place <x> = H[s] into heap slot <s>;
    instead of looking "down" and comparing to children, we look "up" and compare to parents,
      see if <x> needs to shift up the heap.
    """

    H = self.H
    Pairs = self.Pairs

    x = H[s]
    p = Pairs[s]

    i = s
    j = s // 2

    while (j >= 1):
      # keep looking up all the way to the top of the heap, H[1]
      if H[j] <= x:
        break

      #not breaking: heap condition still not met; parent needs to swap places as x "shifts up"
      self.place_in_i(i, H[j], Pairs[j])

      #shifting up generations:
      i = j
      j = j // 2

    self.place_in_i(i, x, p) #break: can place x in slot i while satisfying heap conditions


  def updating_heap(self, x, s, pair):
    """
    Want to update the heap entry H[s] corresponding to <pair>, to value <x>.
    Depending on whether <x> is bigger or smaller than value it is replacing in H[s],
      will need to shiftdown or shiftup.
    """

    old_entry, self.H[s] = self.H[s], x

    # H = self.H
    # Pairs = self.Pairs

    # x = H[s]
    # p = Pairs[s]
    
    self.H[s] = x
    del self.Where[jet_key(*self.Pairs[s])]  #remove old pair
    self.Pairs[s] = pair #plugging in new pair into appropriate slot
    self.Where[jet_key(*pair)] = s #setting 'Where' of new pair

    if (x > old_entry):
      self.shiftdown(s, self.e)
    elif (x < old_entry):
      self.shiftup(s)



  def remove_from_heap(self,s):
    """
    want to remove an element H[s] from heap, by exchanging it with last element of heap and
    shifting up or down as necessary
    """

    leavingpairkey = jet_key(*self.Pairs[s]) #ID of pair that is being removed
    e = self.e
    x = self.H[s]

    self.place_in_i(s, self.H[e], self.Pairs[e]) # will move last element to take the place of removed element

    self.e = e - 1

    #remove the elements to be removed
    self.H.pop()
    self.Pairs.pop()
    del self.Where[leavingpairkey]

    #shift original end-element around to fix the heap
    if s <= self.e:

      if x < self.H[s]:
        self.shiftdown(s, self.e)
      if x > self.H[s]:
        self.shiftup(s)

    

  # new simpler method to remove jet using remove_from_heap
  def remove_jet(self, jet):
    """
    Removes all instances of <jet> from heap;
    from the pairwise distances it is involved in, in the H;
    the pairs it is a part of in Pairs,
    and from the Where and active dicts

    """
    # while True:
    #   i, blank = self.jet_locator(jet, "single")
    #   if i is None: 
    #     return
    #   self.remove_from_heap(i)  #past iteration with inefficient jet_locator method
    
    for k in list(self.active.values()):
      if k is not jet:
        self.remove_from_heap(self.Where[jet_key(jet, k)])
    del self.active(frozenset(jet.id))
    

  def top(self):
    """
    easily access minimal distance and its corresponding pair, H[1] and Pairs[1]
    """

    return self.H[1], self.Pairs[1]



def hac(event_objects, pair_dist_func, timing=True):

  """
  Heap-based Agglomerative Clustering: again based on Kurita outline for agglomerative clustering algorithm

  inputs:

    event_objects (list) - list of input particles to be clustered; expected to be naive.MJet objects

    pair_dist_func (function) - function of two jets j_a, j_b from event_objects, that determines pairwise distance btwn them

    timing (bool) - a new setting, to ask function to also measure time it takes to complete clustering process;

  returns:

    final_jets (list) - list of final-state jets once all clustering is done

    dt (float) - time in seconds it takes to complete clustering process

    
  
  """


  #start the clock 
  t_start = time.perf_counter()
  
  #construct heap with given list of objects for clustering, and a function to determine distance
  distHeap = DistHeap(event_objects, pair_dist_func)

  

  #final state jets
  final_jets = []
  while (distHeap.e >= 1):
    #need to find "top" of the heap, pluck out its distance and the pair within it
    MinDist, MinPairing = distHeap.top()

    #status message print
    print(f"next iteration; {distHeap.e} heap slots remaining...")
  
    #first, determine course of action using naive.epsil_param(j1,j2)
    # 0 -- Merge; 1 -- Drop (the softer jet); 2 -- Isolate (both jets to final state)
  
    cdi_condition = naive.epsil_param(MinPairing[0], MinPairing[1])
  
    if cdi_condition == 2:
      #isolate: need to append both jets into the final_state, remove all instances of the jet
      print("isolating two jets to final state...")
      #remove first jet:

      #stale, involved jet removal method
      # j1_locations = distHeap.jet_locator(MinPairing[0], "all")
      # j1_locations.reverse()

      # for j1_idx in j1_locations:
      #   distHeap.remove_from_heap(j1_idx)

      distHeap.remove_jet(MinPairing[0])

      #remove second jet

      #stale removal method
      # j2_locations = distHeap.jet_locator(MinPairing[1], "all")
      # j2_locations.reverse()

      # for j2_idx in j2_locations:
      #   distHeap.remove_from_heap(j2_idx)

      distHeap.remove_jet(MinPairing[1])

      #append the promoted jets to final state
      final_jets.append(MinPairing[0])
      final_jets.append(MinPairing[1])
      
    elif cdi_condition == 1:
      #drop: need to drop softest jet (in terms of E_T)
      print("dropping one jet...")
      ET1 = (MinPairing[0].E)**2 - (MinPairing[0].pz)**2 #really E_T^1 squared;
      ET2 = (MinPairing[1].E)**2 - (MinPairing[1].pz)**2

      if ET1 < ET2:
        dropjet = MinPairing[0]
      else:
        dropjet = MinPairing[1]

      #stale removal method
      # dropjet_locations = distHeap.jet_locator(dropjet, "all")
      # dropjet_locations.reverse()

      # for dj_idx in dropjet_locations:
      #   distHeap.remove_from_heap(dj_idx)
      distHeap.remove_jet(dropjet)

      
    elif cdi_condition == 0:
      #cluster: merge the two jets, remove one of them, then update elems of another
      print("clustering two jets...")
      new_merged = MinPairing[0] + MinPairing[1]

      #will remove j2:
      # stale removal method
      # j2_locations = distHeap.jet_locator(MinPairing[1], "all")
      # j2_locations.reverse()

      # for j2_idx in j2_locations:
      #   distHeap.remove_from_heap(j2_idx)
      distHeap.remove_jet(MinPairing[1])

      #will update j1 to new_merged, update distance elements
      #another deprecated method;
      # j1_locations, pairingjet = distHeap.jet_locator(MinPairing[0], "single")

      # while j1_locations is not None:
          
      #     #now, find first instance in heap where j1 occurs (in <j1_locations>)

      #     #go to heap, update:

      #   distHeap.updating_heap(naive.d_ij_SIFT(new_merged, pairingjet), j1_locations, [new_merged, pairingjet])
          
      #     #at end, do j1_locations recalc.
      #   j1_locations, pairingjet = distHeap.jet_locator(MinPairing[0], "single")

      #will update j1 to new_merged, update distance elements
      for k in list(self.active.values()):
        if k is MinPairing[0]:
          continue
        s = distHeap.Where[jet_key(MinPairing[0], k)] #found k that isn't j1; need to find where this unique j1,k pairing is
        distHeap.updating_heap(naive.d_ij_SIFT(new_merged, k), s, [new_merged, k]) #replacing distance d_ij(j1, k) --> d_ij

      #replace j1 with new_merged in active pool
      del distHeap.active[frozenset(MinPairing[0].id)]
      distHeap.active[new_merged.id] = new_merged


  #clustering complete: can safely return the final state jet list

  #end the clock
  t_end = time.perf_counter()
  dt = t_end - t_start

  if timing:
    return final_jets, dt
  else:
    return final_jets
        
      

        
      
    
  

  

  
    

    

      
        
        


    
    


      



















