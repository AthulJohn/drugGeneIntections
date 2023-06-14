import numpy as np

def kronrls(k1, k2, y, lmbda=1):
    l_a, q_a = np.linalg.eig(k1)
    l_b, q_b = np.linalg.eig(k2)
    l = np.kron(np.asmatrix(l_a).T,l_b)
    inverse = l / (l + lmbda)
    m1 = q_a.T @ y @ q_b
    m2 = m1 * np.asarray(inverse)
    A = q_a @ m2 @ q_b.T
    return A
         
