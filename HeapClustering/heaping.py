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























