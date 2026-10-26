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
        Where[a_i, a_j] - 


############!!!! Decide what to do for Pair#################


Pair: do we still want to do indices?
      do we want Pair elements to be tuples of the jets themselves?
      do we want to treat this as a second index? We end up limiting communication with <naive>...
        --for now, let's do this but index them the same way we index the input particles;
        --so what, it's a tuple of jet index lists?
        --there's a certain strength to this; not only can you identify each jet uniquely by the end,
          but at every step of the process

Where:
  --can work with preliminary method; 
  --note to self: can use list.index(value), even if value is a list like naive.MJet.id
############################################################

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

 
import naive
  
class DistHeap:

  def __init__(self, dists, ids):
    """
    A heap-based alternative to the matrix-based method of earlier iterations in <naive.py>,
      based on Kurita paper, proposing an efficient, agglomerative clustering algorithm with a heap
      <https://www.sciencedirect.com/science/article/pii/003132039190062A>
    
    Accepts:

    dists (list) -- list of pairwise distances, between pairs of "active" clusters
                  i.e., dist[i] = d_ij(j_a, j_b)
    ids (list) -- list of naive.MJet.id of each pair, with ids[i] = [j_a.id, j_b.id]; same j_a & j_b of dist[i]
    
    DistHeap holds three key arrays:

      H:  -- the "heap" itself
          -- 1-indexed; from H[1], ..., H[e] (going to stick a [None] at the beginning, that shouldn't be considered)
          -- must satisfy:
              * H[i] <= H[2i]; H[i] <= H[2i+1] for i = 1,...,e//2

      Pairs:  -- the "Pairs" list;
              -- Pairs[i] identifies the pair of jets whose pairwise distance is in H slot [i]
              -- does identifying by list of naive.MJet.ids
      
      Where:  -- the "Where" list:
              -- Where([idx1, idx2]) identifies the slot i in H, whose pairwise distance comes from d_ij(j_a,j_b)
              where j_a.id = idx1, j_b.id = idx2
              -- need to somehow make it order independent, both within idx's, and within the input pair [idx1,idx2]
    """
    
    self.e = len(dists) #length of the heap
    self.H = [None] + list(dists)
    self.Pairs = [None] + ids
    self.Where = [None] + list(range(1, len(dists) + 1))

    s = (len(dists) ) // 2 + 1
    while (s > 1):
      s -= 1
      self.shiftdown(s,len(dists)) #start evaluating from the lowest end of the tree back upwards

  def place_in_i(self, i, x, p):
    """
    inputs: i (int), x (float), p ( list( list, list ) )
    places item <x> from pair <p> into heap slot <i>
    """
    self.H[i] = x
    self.Pairs[i] = p
    self.Where[i] = i

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
      self.place_in_i(i, H[j], Pair[j])

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

    

      
        
        


    
    


      



















