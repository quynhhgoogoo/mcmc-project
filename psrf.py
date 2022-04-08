import numpy as np

def compute_psrf(samples):
    
    """
    param samples: 3 dimensional numpy array of shape (M, N, P)
    returns: R-hat values for each model parameter
    """
    b = 2.0
    M, N, P = samples.shape
    
    within_chain_means = samples.mean(axis=1)
    within_chain_vars = W = samples.var(axis=1)
    between_chain_means = samples.mean(axis=(0, 1))
    
    B = N * ((within_chain_means - between_chain_means) ** 2).sum(axis=(0)) / (M - 1)
    
    W = within_chain_vars.mean(axis=0)
    
    V = (N - 1) * W / N + (M + 1) * B / (M * N)
    
    psrf = np.sqrt(V / W)
    
    return psrf