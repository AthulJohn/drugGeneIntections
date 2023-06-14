import numpy as np
from sklearn.model_selection import KFold
from kronrls_mkl.kronrls_mkl_fun import kronrls_mkl

from kronrls_mkl.evaluate_performance import evaluate_performance

def optimize_regcoef(K1, K2, y):
    regcoef_values = np.arange(0, 1.25, 0.25)

    best_regcoef = 0
    best_aupr = 0

    lambda_val = 1
    nfolds = 5
    cv = KFold(n_splits=nfolds, shuffle=True)

    for r in regcoef_values:
        predictions = np.zeros_like(y)

        for train_idx, test_idx in cv.split(y):
            y_train = np.copy(y)
            y_train[test_idx] = 0
            
            regcoef = r
            isinner = 1

            P, _, _ = kronrls_mkl(K1, K2, y_train, lambda_val, regcoef, '', isinner)
            predictions[test_idx] = P[test_idx]
            
        yy = np.copy(y)
        yy[yy == 0] = -1
        stats = evaluate_performance(predictions.flatten(), yy.flatten(), 'classification')

        print('(kronrls): INNER-CV, regcoef={:.2f} - AUPR={:.4f}'.format(r, stats['aupr']))

        if stats['aupr'] > best_aupr:
            best_regcoef = r
            best_aupr = stats['aupr']

    return best_regcoef
