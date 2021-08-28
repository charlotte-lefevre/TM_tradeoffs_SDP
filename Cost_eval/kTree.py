## @package kTree
# Implementation of cost function of ktree and smoothed ktree algorithm

import math
import numpy as np
from utils import skkl, real_skkl


## computes costs of the kTree using Stirling approximation
# For each l value, the number of floors is optimized to have the lowest possible memory and granularty 
# @param npoints int the number of points used to compute on the graph (in other words the precision of the graph)
# @param Remin float min relative size of sub-syndrome
# @param Remax float ditto with max. If equal to 0, it will take the maximal possible value
# @return (M,T,Reused, aused) tuple of arrays, with M containing memories, T times , Reused sub-syndrome relative sizes , and  aused the number of levels 
# If n is the code size, \f$\forall\f$ i, one point corresponds to the kTree instanciated with sub-syndrome size equal to \f$n \times Reused[i]\f$ with aused[i] levels, memory cost \f$M = 2^{nM[i]}\f$, time cost \f$T= 2^{nT[i]}\f$. 
def kTree(R, W,npoints=200,Remax=0, Remin=0.01):
    
    if Remax == 0:
        Remax = W-R
    Res = np.linspace(Remin,Remax, npoints)

    complexities = [];  mems = [] ; Reused = []; aused = []
    
    for Re in Res:
        a = 1
        if (2**a/a <= math.log(2,3)*(R+Re)/Re ):
            while 2**a/a <= math.log(2,3)*(R+Re)/Re:
                a = a +1
            a = a - 1
            complexities.append( Re/a*math.log(3,2)+max(0,skkl(Re,R,W)*math.log(3,2)-Re/a*math.log(3,2)) )
            mems.append(Re/a*math.log(3,2))
            Reused.append(Re)
            aused.append(a)
            
    return (mems,complexities, Reused,aused)

## Computes cost of smoothed kTree using Stirling approximation
# For each l value, the number of floors is optimized to have lowest possible memory and granularity 
# @param npoints int the number of points to compute on the graph (in other words the precision of the graph)
# @param Remin  float minimum sub syndrome relative size
# @param Remax float ditto with max. If equal to 0, it will take the maximal possible value
# @return (M,T,Reused, aused) tuple of arrays, with M containing memories, T times ,  Reused sub-syndrome relative sizes ,  and aused the number of levels 
# If n is the code size, \f$\forall\f$ i, one point corresponds to the smoothed kTree instanciated sub-syndrome size equal to \f$n \times Reused[i]\f$ with aused[i] levels, memory cost \f$M = 2^{nM[i]}\f$, time cost \f$T= 2^{nT[i]}\f$. 
def smoothed_kTree(R, W,npoints=200,Remax=0, Remin=0.01):
    
    if Remax == 0:
        Remax = W-R
    Res = np.linspace(Remin,Remax, npoints)

    complexities = [] ; mems = [] ; Reused = [] ; aused = []
    for Re in Res:
        a = 3
        if (2**(a-1)/(a-1) <= math.log(2,3)*(R+Re)/Re ):
            while (2**(a-1)/(a-1) <= math.log(2,3)*(R+Re)/Re ):
                a = a +1
            a = a - 1
            Lambda = max( Re*math.log(3,2)/(a-2) - (R+Re)/((a-2)*2**(a-1)), (R+Re)/2**(a)) 
            
            complexities.append( Lambda +max(0,skkl(Re,R,W)*math.log(3,2)-Lambda) )
            mems.append(max( (R+Re)/2**(a), Lambda)  )
            Reused.append(Re)
            aused.append(a)
            
    return (mems,complexities, Reused,aused)




## computes real cost of kTree.
# For each l value, the number of floors is optimized to have lowest possible memory and granularity 
# @param npoints int the number of points to compute on the graph (in other words the precision of the graph)
# @param lmin int min size of syndrome subproblem
# @param lmax int ditto with max. If equal to 0, it will take the maximal possible value
# @return (M,T,lused, aused) tuple of arrays, with m containing memories, t times , l sub-syndrome sizes.
# \f$\forall\f$ i, one point corresponds to the kTree instanciated with sub-syndrome size equal to lused[i] with  memory cost \f$M = 2^{m}\f$, time cost \f$T= 2^{T[i]}\f$. 
def real_kTree(n,k, w,npoints=200,lmin = 2, lmax =0):
    if lmax == 0:
        lmax = w-k
    l = lmin
    step = int( (lmax-l)/npoints )

    complexities = [] ; mems = [] ; lused = []
    
    while l < lmax:
        a = 1
        if (2**a/a <= math.log(2,3)*(k+l)/l ):
            while 2**a/a <= math.log(2,3)*(k+l)/l:
                a = a +1
            a = a - 1
            coeff_mem = math.log(2**a*2*(l+(k+l)/2**a),2) # Number of lists times F3 coeff times number of coordinates stored
            coeff_times = math.log((2**a-1)*(2*l + 2*2*l/a) + (n-k-l)*(k+l)*2,2) #  number of considered mergings times (addition cost + merging cost + on-the-fly check)
            complexities.append(coeff_times+ l/a*math.log(3,2)+max(0,real_skkl(n,w,k,l)*n*math.log(3,2)-l/a*math.log(3,2)) )
            mems.append(l/a*math.log(3,2)+coeff_mem)
            lused.append(l)
        l = l + step
            
    return (mems,complexities, lused)


## computes real cost of smoothed kTree.
# For each l value, the number of floors is optimized to have lowest possible memory and granularity 
# @param npoints int the number of points to compute on the graph (in other words the precision of the graph)
# @param lmin int min size of syndrome subproblem
# @param lmax int ditto with max. If equal to 0, it will take the maximal possible value
# @return (M,T,lused) tuple of arrays, with m containing memories, t times , lused sub-syndrome sizes.
# \f$\forall\f$ i, one point corresponds to the smoothed kTree instanciated with sub-syndrome size equal to lused[i] with memory cost \f$M = 2^{M[i]}\f$, time cost \f$T= 2^{T [i]}\f$. 
def real_smoothed_kTree(n,k,w,npoints=200,lmin = 2,lmax =0):
    if lmax == 0:
        lmax = w-k
    l = lmin
    step = int( (lmax-l)/npoints )

    complexities = [] ; mems = [] ; lused = []
    while l < lmax:
        a = 3
        if (2**(a-1)/(a-1) <= math.log(2,3)*(k+l)/l ):
            while (2**(a-1)/(a-1) <= math.log(2,3)*(k+l)/l ):
                a = a +1
            a = a - 1
            S = real_skkl(n,w,k,l)
            Lambda = max( l*math.log(3,2)/(a-2) - (k+l)/((a-2)*2**(a-1)), (k+l)/2**(a)) 
            coeff_mem = math.log(2**a*2*(l+(k+l)/2**a),2) # Number of lists times F3 coeff times number of coordinates stored
            coeff_times = math.log((2**a-1)*(2*l + 2*2*l/a) + (n-k-l)*(k+l)*2,2) #  number of considered mergings times (addition cost + merging cost + on-the-fly check)
            complexities.append( Lambda +max(0,S*n*math.log(3,2)-Lambda) +coeff_times)
            mems.append(max( (k+l)/2**(a), Lambda) +coeff_mem )
            lused.append(l)

        l = l + step
            
    return (mems,complexities, lused)    
    