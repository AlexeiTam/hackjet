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

"""
