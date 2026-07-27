print('we work on heaps here')

class PlacePoint():

    """
    a class of 3D points we want to cluster, with 3D coordinates (x,y,z)
    """

    def __init__(self, x_coord, y_coord, z_coord):

        self.x = x_coord
        self.y = y_coord
        self.z = z_coord


def PointDistance(p1,p2):
    """
    metric in this case is just the 3D euclidian distance between two points
    """
    return (p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2

def DistanceMatrix(ClusterSet, symmetric=True):
    """
    returns a matrix of the distances btwn all of the clusters in input set.

    Inputs:
    
    ClusterSet (list): input set of clusters of length N, whose distances must be evaluated.

    Returns:

    D (list) : an NxN matrix, such that D[i][j] = PointDistance(ClusterSet[i], ClusterSet[j])

    """

    N = len(ClusterSet)

    for i in range(N):

        for j in range(i+1,N):

            D[i][j] = PointDistance( ClusterSet[i], ClusterSet[j] )

    if not symmetric:
        for i in range(1,N):
            for j in range(i):
                D[i][j] = PointDistance( ClusterSet[i], ClusterSet[j] )
    else:
        for i in range(1,N):
            for j in range(i):
                D[j][i] = D[i][j]

    return D

def D_Flattener(D):
    """
    Flattens D into a starting proto-heap D_flat, that has all of the pairwise distances but doesn't yet
    satisfy the heap conditions

    ! Assumes D is symmetric. 
    """

    D_flat = []

    for i in range(len(N)):

        D_flat += D[i][i:] #start from diagonal element, all the way to end

    return D_flat



def shift_down(HeapInput, start, end, Where, PointPair):

    """
    Takes a HeapInput list, and arranges the elements into a heap, which starts with index start, ends with index end.

    Inputs:

    HeapInput (list) = proto-heap, with all of the elements to be organized in a heap.

    start, end = starting and ending indices of output heap.

    Where = array of indices of where heap elements are stored in heap.

    PointPair = array of indices for the pair of clusters associated with the given values in heap.

    """

    #to store starting point, and index where we will consider inserting element into heap

    i = start

    child = 2*i

    #temporarily store





def heapifier():