from Bio import SeqIO
import numpy as np

def spectrum_kernel(seq1, seq2, k):
    """Calculate the spectrum kernel of two protein sequences"""
    kmers1 = set([seq1[i:i+k] for i in range(len(seq1)-k+1)])
    kmers2 = set([seq2[i:i+k] for i in range(len(seq2)-k+1)])
    return len(kmers1 & kmers2)/min(len(kmers1), len(kmers2))

def compute_spectrum_matrix(seqs, k):
    """Compute the spectrum kernel matrix for a list of protein sequences"""
    n = len(seqs)
    kernel_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            kernel_matrix[i][j] = spectrum_kernel(seqs[i], seqs[j], k)
            print(i,j,kernel_matrix[i][j])
            kernel_matrix[j][i] = kernel_matrix[i][j]
    return kernel_matrix


