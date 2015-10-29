# Author: Paul Glenn
# given input data, returns Pade approximant
# p(x1,x2, ..., xN ) / q(x1, x2, ..., xN), where
# p is of degree M and q is of degree N (polynomials) .
#see Pade(...) for input data format. (column vectors) .

import numpy as np
import scipy.optimize


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

    # multivariate pade function of
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

    # residual fn.
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

    # gradient computation
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
    #ab_min = scipy.optimize.minimize(F,ab, method='Newton-CG',jac=cdjac )


    return ab_min['x']

# Remark: I use Nelder-Mead here because
# CG is 'extremely' slow for this - an option would be to modify PHI to
# a vector function and see if that makes it better. (even though
# that means huge Jacobian? ) since it may then be more precise.
#if anyone wants to test this let me know the results!
# other option, of course, is C++...

