import numpy as np

def H_matrix(e):

    n = len(e)  # number of rows
    m = len(e[0]) + 1 # number of columns

    H = np.ones((n,m))

    for i in range(0,n):
        for j in range(0,m-1):
            H[i,j] = -e[i][j]
    
    return H

def H_terms(H):
    """
    returns (H^T H)^-1 H^T
    """

    H_T = np.transpose(H)
    HTH = np.matmul(H_T, H)
    HTH_inv = np.linalg.inv(HTH)
    HTH_invHT = np.matmul(HTH_inv, H_T)

    return HTH_invHT

## Tests
if __name__ == "__main__":
    print(H_matrix([[-1,1],[0,1]]))
    