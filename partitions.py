import numpy as np
from scipy.misc import comb
from itertools import permutations

def accelAsc(n):
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

partitions = accelAsc(3)

def genPowers(partitions, nvars) :
    plist = list(partitions)
    klist = []
    nmax = max( len(l) for l in plist )
    npart = len(plist)

    plist_init = list(plist)
    nterms = comb(nmax + nvars -1 , nmax )

    defects = [ nvars - len(l) for l in plist]
    for i in range(npart) :
        plist[i] = plist[i] + defects[i] * [0]
        perms = list( permutations(plist[i]) )
        unique_perms = list( set(perms) )
        klist.append( unique_perms )


    return klist



genPowers(partitions,5)
