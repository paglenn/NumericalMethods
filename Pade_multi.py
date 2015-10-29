import numpy as np
import scipy.optimize

def cdjac(f,x) :
    # Jacobian from centered differencing
    # for a scalar valued function (gradient)

    eps = 1e-10
    Df = np.zeros_like(x)
    for i in range(x.size) :
        dX = np.zeros_like(x)
        dX[i]= eps
        Df[i] = (f(x+dX)  - f(x-dX)) / (2* eps )

    return Df

def cdjac(F, X):
    # Jacobian from centered differencing
    # for a scalar valued function
    eps = 1e-10
    DF = np.zeros((F.size, X.size))
    for i in range(F.size):
        for j in range(x.size) :
            dX = np.zeros_like(X)
            dX[j] = eps
            DF[i,j] = (F[i](X+dX) - F[i](X-dX)) / (2* eps)

    return DF





def Pade(A,M,N):
    # A : input array, [x0 x1 ... xN y({x}) ]
    # (this is the "function" to be fit)
    # M : order of the numerator
    # N : order of the denominator

    # separate variable array from fn
    ncols = A.shape[1]
    nrows = A.shape[0]
    X = A[:,:-1]
    Y = A[:,-1]

    # pade function in variables
    a = np.ones((ncols -1 ) * M + 1)
    b = np.ones((ncols - 1) * N )
    ab = np.hstack(( a, b) )
    def fn(coeffs , x , M, N) :
        # Pade function will have
        # x: 1-by-n array of variables
        num = coeffs[0]
        den = 1.0
        itr = 0

        for m in range(ncols -1 ) :

            itr = 1 + m * M
            for i in range(M) :
                num += coeffs[itr] * x[m] ** (i+1)
                itr +=1

            itr = (ncols -1 ) * M + 1 + m * N
            # ens: constant term of denominator is 1
            for j in range(N):
                den += coeffs[itr] *x[m] ** (j+1)
                itr += 1

        if den == 0.0 :
            den += 1e-15
        return num / den

    def PHI(coeffs, X,Y , M, N) :
        # coeffs: Pade coefficients
        # X : matrix of all variables
        # Y : column vector of y-values
        resid2 = 0.0 # residual
        for i in range(nrows) :

            x = X[i,:] # one row of X = one set of indep. vars
            y =  Y[i]
            resid2 += (fn(coeffs, x, M, N ) - y) ** 2.

        return resid2 / nrows

    def F(x) :
        return PHI(x, X, Y, M, N)
    #F = lambda coeffs: PHI(coeffs,X,Y,M, N)

    def cdjac(x) :

        eps = 1e-10
        df = np.zeros_like(x)
        for i in range(x.size) :
            dX = np.zeros_like(x)
            dX[i]= eps
            df[i] = (F(x+dX)  - F(x-dX)) / (2* eps )

        #Df = np.outer(df,df)
        return df

    # now find optimal values in ab through CG minimization
    ab_min = scipy.optimize.minimize(PHI, ab,args=((X,Y, M , N)),method =
            'Nelder-Mead', options={'maxfev':3e3,'disp':True,'xtol':1e-6,'ftol':1e-6/nrows})
    #ab_min = scipy.optimize.minimize(F, ab ,method =
    #        'Nelder-Mead', options={'maxfev':3e3,'disp':True,'xtol':1e-6,'ftol':1e-6/nrows})
    #ab_min = scipy.optimize.minimize(F,ab, method='Newton-CG',jac=cdjac )
    # CG is 'extremely' slow for this - an option would be to modify PHI to
    # a vector function and see if that makes it better. (even though
    # that means huge Jacobian? ) since it may then be more precise.
    #if anyone wants to test this let me know the results!

    return ab_min['x']

