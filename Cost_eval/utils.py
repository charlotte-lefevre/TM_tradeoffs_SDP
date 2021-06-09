## @package utils
# Useful functions 

import math
import numpy as np
import matplotlib.pyplot as plt
from sage.functions.log import *
from sage.arith.all import *

## Useful when we have plenty of memories and times and there are many points which are strictly better than others \n
# (m,t) is strictly better than (m',t') if m<m' and t<t'
# @param Mx the array containing all memories 
# @param Ty the array containing all times so that M_x[i] is associated to T_y[i] for all i
# @param col an optional array which is strongly associated to memory and times (ie \f$\forall i\f$ col[i] is associated to M_x[i] and T_y[i])
# @returns epurated arrays (M_x, T_y) and if an array col is provided it also returns cols epurated
def clean_array(M_x,T_y, col = [], nocol = False):
    l = len(M_x)
    i = 0
    collen = len(col)
    while i < l:
        j = i+1; 
        while(j < l):
            if (M_x[i] < M_x[j] and T_y[i] < T_y[j]):
                M_x = np.append(M_x[0:j], M_x[j+1:l])
                T_y = np.append(T_y[0:j], T_y[j+1:l])
                if collen != 0:
                    col = np.append(col[0:j], col[j+1:l])
                l = l-1
            elif (M_x[i] > M_x[j] and T_y[i] > T_y[j]):
                M_x = np.append(M_x[0:i], M_x[i+1:l])
                T_y = np.append(T_y[0:i], T_y[i+1:l])
                if collen != 0:
                    col = np.append(col[0:i], col[i+1:l])
                j = i + 1
                l = l-1
            else:
                j = j + 1
        i = i + 1
    if collen == 0:
        return (M_x, T_y)
    return (M_x, T_y, col)




##  Ternary entropy
# Useful to compute approximations
def hq(x):
    if x > 1:
        raise Exception ("Ternary entropy : out ouf bound input : x = "  + str(x))
    elif x < 0:
        raise Exception("Ternary entropy : out ouf bound input : x = "  + str(x))
    if x == 1:
        return -x*math.log(x/2,3)
    if x == 0:
        return -(1-x)*math.log(1-x,3)
    return -(1-x)*math.log(1-x,3)-x*math.log(x/2,3)

## Crutial function to evaluate asymptotic time cost of the subset sum. \n 
# The computation is based on Stirling approximation
# @param Re float sub problem rate : Re is such that l = Re * n
# @param R float code rate : R is such that k = R * n
# @param W float relative rate : W such that w = W*n
# @returns float approximated log3 of needed solutions for probabilistic step divided by n
def skkl(Re,R, W):
    C = (W-R-Re)/(1-R-Re)
    top = (1-R-Re)*hq(C) + Re
    bottom = min(1-R, hq(W))
    return(-(top-bottom))


## Exact version of skkl function
# It computes the real average expected number of soluions for the probabilistic step
# @param n int length of the code
# @param w int target weight
# @param k int code dimension
# @param l int sub-target size
# @returns float log3 of needed solutions for probabilistic step divided by n
def real_skkl(n,w,k,l):
    val = min(3**(n-k-l), binomial(n,w)*2**(w)*3**(-l)) /(binomial(n-k-l,w-k-l)*2**(w-k-l)) 
    return log(val,3)/n


## Magic sequence \n 
# Important quantity for the dissection \n 
# Here, the m parameter for the dissection is assumed to be equal to 1 
# @param int x index of desired sequence (or in other words the desired gain of associated dissection)
# @returns int magic[x]
def magic(x):
    return ( 0.5*(1+x)*(2+x) + 1)

## Magic sequence with any m parameter \n 
# Every i parameter provides a sub sequence, where the gain is increased by m at each step.
# @param  m integer the memory unit. If base list cardinal is \f$3^n\f$ then maximal allowed memory is \f$3^{nm}\f$ 
# @param i integer is such that i < m and the base dissection is a i 2m+2i meet-in-the-middle (gain i)
# @param x integer index in the sub sequence with fixed i. \n 
# @returns magic[m,i,x] integer
def magic_m(m,i,x):
    return magic(x) * m + (x+2)*i

## Generalized magic sequence \n 
# @param r integer : gives the dissection degree
# @param m the memory unit
# @returns (r_magic, g) such that \n
# * r_magic is the biggest number in m magic sequence such that r_magic > r
# * g is the gain associated to r_magic. Therefore it is the biggest possible expectable gain associated to r-disssection
def find_best_magic(m,r):
    if r < 2*m:
        return (0,0)
    l = []
    g = []
    for i in range(0,m):
        #Look at the coefficients of degree 2 equation
        a = m/2
        b = i + 3*m/2
        c = 2*m+2*i-r
        delta = b*b - 4*a*c
        if delta >= 0 :
            sol = math.floor((-b + math.sqrt(delta)) /(2*a))
            mag = magic_m(m,i,sol)
            l = np.append(l, i + sol * m)
            g = np.append(g, mag)
    return (g[np.argmax(l)], np.max(l))


## Given a r-dissection and m memory parameter, returns the associated gain
def gain(r,m):
    return find_best_magic(m,r)[1]


## computes asymptotic cost of dissection
# @param r int dissection order
# @returns c float such that r-dissection with memory \f$3^n\f$ costs \f$3^{cn}\f$
def expo_complexity_dissection(r):
    _,gain  = find_best_magic(1,r) #rmagic has complexity rmagic-gain-1
    return r-int(gain)-1

    