import numpy as np
from scipy.optimize import minimize_scalar,minimize
from kronrls_mkl.combine_kernels import combine_kernels

from kronrls_mkl.kronrls import kronrls
from kronrls_mkl.optimize_lambda import optimize_lambda

from sklearn.model_selection import KFold

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

            P = kronrls_mkl(K1, K2, y_train, lambda_val, regcoef, 20, isinner)
            predictions[test_idx] = P[test_idx]
            
        yy = np.copy(y)
        yy[yy == 0] = -1
        stats = evaluate_performance(predictions.flatten(), yy.flatten(), 'classification')
        if stats['aupr'] > best_aupr:
            best_regcoef = r
            best_aupr = stats['aupr']
    return best_regcoef


def optimize_weights(x0, fun):
    """
    Optimizes the weights using fmincon.
    
    Args:
    x0: Initial weights.
    fun: Function to optimize.
    
    Returns:
    Optimized weights and the function value.
    """
    n = len(x0)
    Aineq   = None
    bineq   = None
    Aeq     = np.ones((1,n))
    beq     = 1
    LB      = np.zeros(n)
    UB      = np.ones(n)

    options = {'disp': True}
    res = minimize(fun, x0, method='SLSQP', bounds=list(zip(LB, UB)), constraints=[{'type': 'eq', 'fun': lambda x: np.dot(Aeq, x) - beq}], options=options)
    return res.x, res.fun

def obj_function(w, u, Ma, lambda_val, regcoef):
    V = combine_kernels(w, Ma)
    J = np.linalg.norm(u - V.T)/(2*lambda_val*np.size(V)) + regcoef*(np.linalg.norm(w,2))**2
    return J

def kronrls_mkl(K1, K2, y, lambda_, regcoef=0.2, maxiter=20, isinner=False):
    assert K1.ndim == 3 and K2.ndim == 3, "K1 and K2 must be order 3 tensors"
    # print(K1.shape,K2.shape)
    
    nA =len(K1)
    nB = len(K2)
    
    # initialization of kernel weights (uniform)
    alpha = np.full(nA, 1 / nA)
    beta = np.full(nB, 1 / nB)
  
    iter_ = 0
    incr = 1000
    limit_incr = 0.0100
    combFval = np.zeros(maxiter)
    
    if not isinner:
        regcoef = optimize_regcoef(K1, K2, y)
    
    # iterative steps: optimize weights
    while iter_ < maxiter and abs(incr) > limit_incr:
        iter_ += 1
        
        # Step 1
        K1_comb = combine_kernels(alpha, K1)
        K2_comb = combine_kernels(beta, K2)
        A = kronrls(K1_comb, K2_comb, y, lambda_)
        
        u = (y - lambda_ * A / 2)
        
        # Step 2
        Ma = np.zeros((len(K1),len(y[0]),len(y)))
        for i in range(nA):
            Ma[i] =K2_comb.T@A.T@K1[i, :, :]
        falpha = lambda alpha: obj_function(alpha, u, Ma, lambda_, regcoef)
        
        # Optimal Alpha
        res_alpha = minimize(falpha, x0=alpha, method='L-BFGS-B')
        x_alpha, fval_alpha = res_alpha.x, res_alpha.fun
        
        # Step 3
        Mb = np.zeros((len(K2),len(y[0]),len(y)))
        for i in range(nB):
            Mb[i] = K2[i, :, :].T@A.T@ K1_comb
        
        fbeta = lambda beta: obj_function(beta, u, Mb, lambda_, regcoef)
        
        # Optimal Beta
        res_beta =  minimize(fbeta, x0=beta, method='L-BFGS-B')
        x_beta, fval_beta = res_beta.x, res_beta.fun
        
        # update values for next iteration
        alpha = x_alpha 
        beta = x_beta 
        
        combFval[iter_ - 1] = fval_alpha + fval_beta
        if iter_ > 1:
            incr = combFval[iter_ - 2] - combFval[iter_ - 1]
        
        # if not isinner:
            # print(f"ITER={iter_} \tlambda={lambda_:.2f} \tregcoef={regcoef:.1f} - (fval_alpha={fval_alpha:.4f}, fval_beta={fval_beta:.4f})")
    
    # Perform predictions
    K1_comb = combine_kernels(alpha, K1)
    K2_comb = combine_kernels(beta, K2)
    
    # build final model with best weights/lambda
    if not isinner:
        lambda_ , best_aupr = optimize_lambda(K1_comb, K2_comb, y)
    
    A = kronrls(K1_comb, K2_comb, y,lambda_)
    
    return A
