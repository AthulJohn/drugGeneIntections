import numpy as np
from scipy.io import savemat, loadmat
from sklearn.model_selection import KFold
from kronrls_mkl.evaluate_performance import aupr_score
from kronrls_mkl.kronrls import kronrls

def optimize_lambda(K1, K2, y):
    exponents = np.arange(-15, 31, 5)
    bestaupr = 0
    best_lambda = 1
    nfolds = 5
    
    crossval_idx = KFold(n_splits=nfolds, shuffle=True).split(y)
    
    for e in exponents:
        predictions = np.zeros_like(y)
        for train_idx, test_idx in crossval_idx:
            y_train = np.copy(y)
            y_train[test_idx] = 0
            
            lambda_ = (2.0)**e
            
            y2 = kronrls(K1, K2, y_train, lambda_)
            predictions[test_idx] = y2[test_idx]
        
        yy = np.copy(y)
        yy[yy==0] = -1
        
        aupr = aupr_score(yy.flatten(), predictions.flatten())
        
        
        if aupr > bestaupr:
            best_lambda = lambda_
            bestaupr = aupr
            
    return best_lambda, bestaupr
