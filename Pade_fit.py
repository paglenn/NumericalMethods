# Author: Paul Glenn
# given input data, returns Pade approximant
# p(x1,x2, ..., xN ) / q(x1, x2, ..., xN), where
# p is of degree M and q is of degree N (polynomials) .
#see Pade(...) for input data format. (column vectors) .
import numpy as np
import scipy.optimize
from scipy.misc import comb
from itertools import permutations


def accelAsc(n):
    # generate integer partitions of n
    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2*x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]

def genPowers(nvars, nmax) :
    # generate list of exponents
    partitions = accelAsc(nmax)
    plist = list(partitions)
    klist = []
    npart = len(plist)

    plist_init = list(plist)
    #print plist_init
    nterms = comb(nmax + nvars -1 , nmax )

    defects = [ nvars - len(l) for l in plist]
    for i in range(npart) :
        if defects[i] >= 0 :
            plist[i] = plist[i] + defects[i] * [0]
            perms = list( permutations(plist[i]) )
            unique_perms = list( set(perms) )
            klist += unique_perms


    #rint klist
    return klist

def product(myList) :
    p = 1.0
    for item in myList :
        p *= item
    return p

def Pade_fit(A,M,N):
    # A : input array, [x0 x1 ... xN y({x}) ]
    # (this is the "function" to be fit)
    # M : desired order of the numerator
    # N : desired order of the denominator

    # separate variable array from fn
    ncols = A.shape[1]
    nrows = A.shape[0]
    X = A[:,:-1]
    lenX = ncols - 1
    xitr = range(lenX)
    Y = A[:,-1]

    # multivariate pade function of
    aExps = [tuple([0]*lenX) ] ; bExps = []
    for n in range(1,1+ max(M,N)) :
        if n <= M :
            aExps.extend( genPowers(lenX, n) )
        if n <= N :
            bExps.extend(  genPowers(lenX, n) )

    lena = len(aExps)
    lenb = len(bExps)
    a = np.zeros(lena)
    b = np.zeros(lenb)

    c = np.hstack(( a, b) )

    def fn(coeffs , x, y) :
        # Pade function will have
        # x: 1-by-n array of variables
        xitr = range(len(x))
        num = 0.0
        den = 1.0

        for n in range(max(lena,lenb)):
            if n < lena :
                exps = aExps[n]
                mon = [x[i] ** exps[i] for i in xitr ]
                num += coeffs[n] * product(mon)

            if n < lenb:
                exps = bExps[n]
                mon = [x[i] ** exps[i] for i in xitr ]
                den += coeffs[lena + n] * product(mon)

        return num - den * y

    # residual fn.
    def F(coeffs) :
        # coeffs: Pade coefficients
        # X : matrix of all variables
        # Y : column vector of y-values
        resid2 = 0.0 # residual
        for i in range(nrows) :

            x = X[i,:] # one row of X = one set of indep. vars
            y =  Y[i]
            resid2 +=  fn(coeffs, x,y) ** 2.

        #print resid2
        return resid2 / nrows

    # gradient computation via centered differencing
    def cdjac(x) :

        eps = 1e-10
        df = np.zeros_like(x)
        for i in range(x.size) :
            dX = np.zeros_like(x)
            dX[i]= eps
            df[i] = (F(x+dX)  - F(x-dX)) / (2* eps )

        #Df = np.outer(df,df)
        return df


    # now find optimal values in ab through minimization
    ab_min = scipy.optimize.minimize(F, c, method = 'Nelder-Mead',
			options={'disp':True,'ftol':1e-5})
    #ab_min = scipy.optimize.minimize(F,c, method='Newton-CG',jac=cdjac )


    return ab_min['x'], aExps, bExps

# Remark: I use Nelder-Mead here because
# CG is 'extremely' slow for this - an option would be to modify PHI to
# a vector function and see if that makes it better. (even though
# that means huge Jacobian? ) since it may then be more precise.
#if anyone wants to test this let me know the results!
# other option, of course, is C++...

def Pade_eval(x, coeffs, aExps, bExps) :
    # Pade function will have
    # x: 1-by-n array of variables
    xitr = range(len(x))
    lena = len(aExps)
    lenb = len(bExps)
    num = 0.0
    den = 1.0

    for n in range(max(lena,lenb)):
        if n < lena :
            exps = aExps[n]
            mon = [x[i] ** exps[i] for i in xitr ]
            num += coeffs[n] * product(mon)

        if n < lenb:
            exps = bExps[n]
            mon = [x[i] ** exps[i] for i in xitr ]
            den += coeffs[lena + n] * product(mon)

    return num/den


