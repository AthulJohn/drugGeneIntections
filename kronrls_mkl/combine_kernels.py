import numpy as np

def combine_kernels(weights, kernels):
    # length of weights should be equal to length of matrices
    n = len(weights)
    result = np.zeros(kernels[0,:,:].shape)
    
    for i in range(n):
        result = result + weights[i] * kernels[i,:,:]
    
    return result