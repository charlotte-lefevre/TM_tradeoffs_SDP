## @package kTree
# Implementation of cost function of ktree and smoothed ktree algorithm

import math
import numpy as np
from utils import skkl, real_skkl

## computes cost of kTree using Stirling approximation
# For each l value, the number of floors is optimized to have lowest possible memory and granulairty 
# @param R float code rate
# @param W float relative weight
# @param npoints int the number of points to compute on the graph
# @param Remin float min coefficient of syndrome subproblem
# @param Remax float ditto with max. If equal to 0, it will take the maximal possible value
# @return (m,t,Reused, aused) tuple of arrays, with m containing memories, t times , Reused syndrome sub-problem coefficients ,  aused number of levels 
# If n is the code size, \f$\forall\f$ i, one point corresponds to kTree instanciated with target sub syndrome coefficient equal to Reused[i] with aused[i] levels memory cost \f$M = 2^{nm}\f$, time cost \f$MT= 2^{nt}\f$. 
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
# For each l value, the number of floors is optimized to have lowest possible memory and granulairty 
# @param R float code rate
# @param W float relative weight
# @param npoints int the number of points to compute on the graph
# @param Remin float min coefficient of syndrome subproblem
# @param Remax float ditto with max. If equal to 0, it will take the maximal possible value
# @return (m,t,Reused, aused) tuple of arrays, with m containing memories, t times , Reused syndrome sub-problem coefficients ,  aused number of levels 
# If n is the code size, \f$\forall\f$ i, one point corresponds to kTree instanciated with target sub syndrome coefficient equal to Reused[i] with aused[i] levels memory cost \f$M = 2^{nm}\f$, time cost \f$MT= 2^{nt}\f$. 
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
# For each l value, the number of floors is optimized to have lowest possible memory and granulairty 
# @param n int length of the code
# @param w int target weight
# @param k int code dimension
# @param l int sub-target size
# @param npoints int the number of points to compute on the graph
# @param lmin int min size of syndrome subproblem
# @param lmax int ditto with max. If equal to 0, it will take the maximal possible value
# @return (m,t,Reused, aused) tuple of arrays, with m containing memories, t times , Reused syndrome sub-problem coefficients
# \f$\forall\f$ i, one point corresponds to kTree instanciated with target sub syndrome coefficient equal to Reused[i] with  memory cost \f$M = 2^{m}\f$, time cost \f$MT= 2^{t}\f$. 
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
            coeff_mem = math.log(2**a*2*n,2) # Number of lists times vector size times F3 coeff
            coeff_times = math.log((2**a-1)*2*n,2) #  number of considered mergings times memory coeff
            complexities.append(coeff_times+ l/a*math.log(3,2)+max(0,real_skkl(n,w,k,l)*n*math.log(3,2)-l/a*math.log(3,2)) )
            mems.append(l/a*math.log(3,2)+coeff_mem)
            lused.append(l)
        l = l + step
            
    return (mems,complexities, lused)


## computes real cost of smoothed kTree.
# For each l value, the number of floors is optimized to have lowest possible memory and granulairty 
# @param n int length of the code
# @param w int target weight
# @param k int code dimension
# @param l int sub-target size
# @param npoints int the number of points to compute on the graph
# @param lmin int min size of syndrome subproblem
# @param lmax int ditto with max. If equal to 0, it will take the maximal possible value
# @return (m,t,Reused) tuple of arrays, with m containing memories, t times , Reused syndrome sub-problem coefficients 
# \f$\forall\f$ i, one point corresponds to kTree instanciated with target sub syndrome coefficient equal to Reused[i] with memory cost \f$M = 2^{m}\f$, time cost \f$MT= 2^{t}\f$. 
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
            coeff_mem = math.log(2**a*2*n,2) # Number of lists times vector size times F3 coeff
            coeff_times = math.log((2**a-1)*2*n,2) #  number of considered mergings times memory coeff
            complexities.append( Lambda +max(0,S*n*math.log(3,2)-Lambda) +coeff_times)
            mems.append(max( (k+l)/2**(a), Lambda) +coeff_mem )
            lused.append(l)

        l = l + step
            
    return (mems,complexities, lused)    
    