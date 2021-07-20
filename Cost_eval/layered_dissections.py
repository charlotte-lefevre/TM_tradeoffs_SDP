## @package layered_dissections
# Computes the costs of layered dissections
# Depending on the input, the applicable algorithm changes.

import numpy as np
import math


from utils import *
from dissection import find_best_magic


""" 
Determines if we can apply the algorithms
"""

## Determines whether base layered ddissection can be applied 
# @return true iff no granularity, nor list cardinal problems
def is_classic_OK(R,Re,r,c,h, Skkl):
    m = 1/(r + (h-1)*(r-1) - c)
    return (c*m*Re  < Skkl and (m * r**h < (Re+R)/Re * math.log(2,3)) )


##  instanciated version of is_classic_OK
# @return true iff no granularity, nor list cardinal problems
def real_is_classic_OK(n,k,l,r,c,h, Skkl):
    m = 1/(r + (h-1)*(r-1) - c)
    return (c*m*l  < Skkl*n and (m * r**h < (l+k)/l * math.log(2,3)) )

## Determines whether smoothed tree can be applied
# @return true iff Smoothing can be applied and no granularity problems
def is_smoothed_OK(R,Re,r,c,h, Skkl): 
    denominateur = r + (h-2)*(r-1) - c -1
    Lambda = math.log(3,2)*Re/denominateur - (R+Re)/(denominateur*r**(h-1))
    mem = math.log(2,3)*Lambda
    m = 1/(r + (h-1)*(r-1) - c)
    return mem * c < Skkl and (m * r**h > (Re+R)/Re * math.log(2,3)) and math.log(3,2)*Re <= (Re+R)*(r-c+(h-2)*(r-1) )/r**(h-1) and Re >= (R+Re)*math.log(2,3)/r**(h-1)

##  instanciated version of is_smoothed_OK
# @return true iff Smoothing can be applied and no granularity problems
def real_is_smoothed_OK(n,k,l,r,c,h, Skkl): 
    denominateur = r + (h-2)*(r-1) - c -1
    Lambda = math.log(3,2)*l/denominateur - (k+l)/(denominateur*r**(h-1))
    mem = math.log(2,3)*Lambda
    m = 1/(r + (h-1)*(r-1) - c)
    return mem * c < Skkl*n and (m *l * r**h > (l+k) * math.log(2,3)) and math.log(3,2)*l <= (l+k)*(r-c+(h-2)*(r-1) )/r**(h-1) and l >= (k+l)*math.log(2,3)/r**(h-1)

## Determines whether calibrated tree can be applied
# @return true iff no list cardinal problems and no problem with dissection granularity
def is_calibrated_OK(R,Re,r,c,h, Skkl):
    _,gainX  = find_best_magic(1,r)
    _,gainOfgain = find_best_magic(1,gainX+1)
    MaxAlpha = r - 2*gainX - 1 - gainOfgain
    alpha = (Re*c+Skkl*(c-r-(h-1)*(r-1) ) )/(Re - (h-1)*Skkl )
    if alpha >=0 and alpha <= MaxAlpha: 
        M = 1 /(r + (h-1)*(r-1) - c - alpha * (h-1) )
        return (M * r**h < (Re+R)/Re * math.log(2,3))
    else:
        return False


## Instanciated version of is_calibrated_OK
# @return true iff no list cardinal problems and no problem with dissection granularity
def real_is_calibrated_OK(n,k,l,r,c,h, Skkl):
    _,gainX  = find_best_magic(1,r)
    _,gainOfgain = find_best_magic(1,gainX+1)
    MaxAlpha = r - 2*gainX - 1 + gainOfgain
    alpha = (l*c+Skkl*n*(c-r-(h-1)*(r-1) ) )/(l - (h-1)*Skkl*n )
    if alpha >=0 and alpha <= MaxAlpha: 
        M = 1 /(r + (h-1)*(r-1) - c - alpha * (h-1) )
        return (M *l  < (l+k)/r**h * math.log(2,3))
    else:
        return False




""" 
Applies the algorithms and returns the time and the memory 
"""

## Applies layered dissection algorithm  
# @return (m,t) log2 of memory / complexity
# @attention no check is done about the validity of the algorithm application : you mus call prior to this algorith is_classic_OK
def apply_classic(R,Re,r,c,h, Skkl):
    m = Re/(r + (h-1)*(r-1) - c)
    return (m*math.log(3,2), math.log(3,2)*max(Skkl, c*m))

## instanciated version of apply_classic
# @return (m,t) log2 of memory / complexity
# @attention no check is done about the validity of the algorithm application : you must call prior to this algorith real_is_classic_OK
def real_apply_classic(n,k,l,r,c,h, Skkl):
    m = l/(r + (h-1)*(r-1) - c)
    g = gain(r,1)
    coeff_mem = math.log(r**h*2*n,2) # Number of lists times vector size times F3 coeff
    coeff_times = math.log((g+1)*(r**h-1)/(r-1)*2*n,2) # Number of main merges done in the dissection times number of considered dissections times memory coeff
    return (m*math.log(3,2)+coeff_mem, math.log(3,2)*max(n*Skkl,c*m)+coeff_times)



## Applies smoothed algorithm  
# @return (m,t) log2 of memory / complexity
# @attention no check is done about the validity of the algorithm application : you mus call prior to this algorith is_smoothed_OK
def apply_smoothed(R,Re,r,c,h, Skkl):
    denominator = r + (h-2)*(r-1) - c -1
    Lambda = math.log(3,2)*Re/denominator - (R+Re)/(denominator*r**(h-1))
    return (Lambda, max(Skkl*math.log(3,2), c*Lambda))

## instanciated version of apply_smoothed
# @return (m,t) log2 of memory / complexity
# @attention no check is done about the validity of the algorithm application : you must call prior to this algorith real_is_smoothed_OK
def real_apply_smoothed(n,k,l,r,c,h, Skkl):
    denominator = r + (h-2)*(r-1) - c -1
    Lambda = math.log(3,2)*l/denominator - (k+l)/(denominator*r**(h-1))
    g = gain(r,1)
    coeff_mem = math.log(r**h*2*n,2) # Number of lists times vector size times F3 coeff
    coeff_times = math.log((g+1)*(r**h-1)/(r-1)*2*n,2) # Number of main merges done in the dissection times number of considered dissections times memory coeff
    return (Lambda+coeff_mem, max(n*Skkl*math.log(3,2),c*Lambda)+coeff_times)
   
## Applies calibrated algorithm 
# @param S desired number of solutions. Here it is not necessarily equal to Skkl, as this can be used by memBlowup 
# @param t the target coefficient such that target size is equal to Re*t (so t = 1 is classical calibration)
# @return (m,t) log2 of memory / complexity
# @attention no check is done about the validity of the algorithm application : you must call prior to this algorith is_calibrated_OK
def apply_calibration(R,Re,r,c,h,S,t):
    _,gainX  = find_best_magic(1,r)
    _,gainOfgain = find_best_magic(1,gainX+1)
    MaxAlpha = r - 2*gainX - 1 + gainOfgain
    Left = t*Re - S*(h-1)
    Right = -S *(r-c) - S*(h-1)*(r-1)+t*Re*c
    if Left >=0:
        alpha = max(Right/Left,0)
        if alpha > MaxAlpha:
            return -1,-1,-1
    elif Right >=0:
        return -1,-1,-1
    else:
        alpha = 0
    M = (t*Re) /(r + (h-1)*(r-1-alpha) - c )
    return (M*math.log(3,2), S*math.log(3,2),alpha)


## Instanciated version of apply_calibration
# @param S desired number of solutions. Here it is not necessarily equal to Skkl, as this can be used by memBlowup 
# @param t the target coefficient such that target size is equal to l*t (so t = 1 is classical calibration)
# @param Asymptotic if False computes memory / time cost in the binary metrics 
# (asymptotic = True is useful when this algorithm is called by memBlowup - as the binary metrics conversion is done at the very end)
# @return (m,t) log2 of memory / complexity
# @attention no check is done about the validity of the algorithm application : you must call prior to this algorith real_is_calibrated_OK
def real_apply_calibration(n,k,l,r,c,h,S,t, asymptotic=True):
    _,gainX  = find_best_magic(1,r)
    _,gainOfgain = find_best_magic(1,gainX+1)
    MaxAlpha = r - 2*gainX - 1 + gainOfgain
    Left = t*l - S*n*(h-1)
    Right = -S *n*(r-c) - S*n*(h-1)*(r-1)+t*l*c
    if Left >=0:
        alpha = max(Right/Left,0)
        if alpha > MaxAlpha:
            return -1,-1,-1
    elif Right >=0:
        return -1,-1,-1
    else:
        alpha = 0
    M = (t*l) /(r + (h-1)*(r-1-alpha) - c )
    if asymptotic:
        return (M*math.log(3,2), S*n*math.log(3,2),alpha)
    else:
        g = gain(r,1)
        coeff_mem = math.log(r**h*2*n,2) # Number of lists times vector size times F3 coeff
        coeff_times = math.log((g+1)*(r**h-1)/(r-1)*2*n,2) # Number of main merges done in the dissection times number of considered dissections times memory coeff
        return (M*math.log(3,2)+coeff_mem, S*n*math.log(3,2)+coeff_times,alpha)

    
## Applies the combined algorithm, with a mix of smoothing and solution calibration.
# @param nbtargets int the bigger, the better precision but the longer time it takes
# @return (m,t) log2 of memory / complexity or -1,-1 if algorithm is not applicable
def apply_mem_blowup(R,Re,r,c,h, Skkl, nbtargets=80):
        TrackMem = -1
        # MaxAlphabeta computation
        _,gainX  = find_best_magic(1,r)
        _,gainOfgain = find_best_magic(1,gainX+1)
        MaxAlphaBeta = r - 2*gainX - 1 + gainOfgain
        cardLeaves = math.log(2,3)*(R+Re)/r**h
        # t is the target of the first floor
        # its maximum value  corresponds to the base algorithm target 
        # We start with biggest possible target because it is the one that gives us the best results
        for t in np.linspace(0, (r-1)/(r + (h-1)*(r-1) - c) , nbtargets)[::-1]:
            cardWanted, _,alpha = apply_calibration(R,Re,r,c,h-1,Skkl,(1-t))
            if cardWanted <0: #UNSAT
                continue 
                
            cardWanted = cardWanted/(math.log(3,2))         

            # Compute bounds for beta : need beta_lower < beta < beta_upper Explanations : 
            # # MaxAlphaBeta : max dissection granularity
            # # (cardWanted+t*Re)/cardLeaves : want that card obtained from first dissection 
            # greater than cardWanted
            # #  c - cardWanted/cardLeaves * (c-alpha) : time complexity of upper part dominates 
            # Also, the bigger is beta, the better
            beta_upper = min(MaxAlphaBeta, r - (cardWanted+t*Re)/cardLeaves)
            beta_lower = c - cardWanted/cardLeaves * (c-alpha)
            
            if beta_upper < beta_lower or beta_upper <0: #UNSAT
                continue
                
            # in such case we have  beta = beta_upper
            TrackMem = cardWanted*math.log(3,2)
            
            # we exit the loop because this is the maximum t such that we have solved the problem
            break 
            

        # No point has been found for this problem
        if TrackMem == -1:
            return -1,-1
        # Every computed point is in O(1) for this Re coeff, so we just need to take the smallest memory 
        return(TrackMem, Skkl*math.log(3,2))




## Applies the combined algorithm, with a mix of smoothing and solution calibration.
# @param nbtargets int the bigger, the better precision but the longer time it takes
# @return (m,t) log2 of memory / complexity or -1,-1 if algorithm is not applicable
def real_apply_mem_blowup(n,k,l,r,c,h, Skkl, nbtargets):
        TrackMem = -1

        # MaxAlphabeta computation
        _,gainX  = find_best_magic(1,r)
        _,gainOfgain = find_best_magic(1,gainX+1)
        MaxAlphaBeta = r - 2*gainX - 1 + gainOfgain
        
        # t is the target of the first floor
        # its maximum value  corresponds to the base algorithm target 
        # We start with biggest possible target because it is the one that gives us the best results
        for t in np.linspace(0, (r-1)/(r + (h-1)*(r-1) - c) , nbtargets)[::-1]:
            cardLeaves = math.log(2,3)*(k+l)/r**h
        
            cardWanted, _,alpha = real_apply_calibration(n,k,l,r,c,h-1,Skkl,1-t)
            if cardWanted <0: #UNSAT
                continue 
             
            cardWanted = cardWanted/(math.log(3,2))         

            # Compute bounds for beta : need beta_lower < beta < beta_upper Explanations : 
            # # MaxAlphaBeta : max dissection granularity
            # # (cardWanted+t*Re)/cardLeaves : want that card obtained from first dissection 
            # greater than cardWanted
            # #  c - cardWanted/cardLeaves * (c-alpha) : time complexity of upper part dominates 
            # Also, the bigger is beta, the better
            beta_upper = min(MaxAlphaBeta, r - (cardWanted+t*l)/cardLeaves)
            beta_lower = c - cardWanted/cardLeaves * (c-alpha)
            
            if beta_upper < beta_lower or beta_upper <0: #UNSAT
                continue
                
            # in such case we have  beta = beta_upper
            
            TrackMem = cardWanted*math.log(3,2)
            
            # we exit the loop because this is the maximum t such that we have solved the problem
            break 
            

        # No point has been found for this problem
        if TrackMem == -1:
            return -1,-1

        g = gain(r,1)
        coeff_mem = math.log(r**h*2*n,2) # Number of lists times vector size times F3 coeff
        coeff_times = math.log((g+1)*(r**h-1)/(r-1)*2*n,2) # Number of main merges done in the dissection times number of considered dissections times memory coeff
        return (TrackMem+coeff_mem, Skkl*n*math.log(3,2)+coeff_times)
        # Every computed point is in O(1) for this Re coeff, so we just need to take the smallest memory 
        return(TrackMem, Skkl*math.log(3,2))



## The main function to compute the times / memories.
# @param Remin float minimum sub syndrome relative size 
# @param Remax float ditto with max
# @param nbRes int number of computed points (linked to the precision)
# @param nbtarget int number of computed targets for memBlowUp
# @return (M,T, col, res) tuple of arrayx with M[i] (resp. T[i]) log2 of memory (resp. time)) cost relative to code length
# col keeps track of the used algorithm and res the target coefficient used.
# In other words M[i] T[i] is a time/memory cost with algorithm col[i] and l coefficient res[i]
def Apply(r, h, R, W, Remin, Remax, nbRes = 100, nbtarget =80):
    c = expo_complexity_dissection(r)
    Res = np.linspace(Remin, Remax, nbRes)
    TrackMem  = [] ; TrackTimes = [] ; TrackCols =[] ; TrackRes = []
    for Re in Res:
        Skkl = skkl(Re, R,W)
        if is_classic_OK(R,Re,r,c,h, Skkl): # check if classic algo applies
            (mem,time) = apply_classic(R,Re,r,c,h, Skkl)
            TrackMem.append(mem)
            TrackTimes.append(time)
            TrackCols.append("red"); TrackRes.append(Re)
        elif is_calibrated_OK(R,Re,r,c,h, Skkl): # check if no granularity problem
            (mem,time,_) = apply_calibration(R,Re,r,c,h,Skkl,1)
            TrackMem.append(mem)
            TrackTimes.append(time)
            TrackCols.append("blue") ;TrackRes.append(Re)
        elif is_smoothed_OK(R,Re,r,c,h, Skkl): # check if smoothing applies
            (mem,time) = apply_smoothed(R,Re,r,c,h, Skkl)
            TrackMem.append(mem)
            TrackTimes.append(time)
            TrackCols.append("pink") ;TrackRes.append(Re)
        else: # try memblowup
            (mem,time) = apply_mem_blowup(R,Re,r,c,h, Skkl,nbtarget)
            if mem != -1:
                TrackMem.append(mem)
                TrackTimes.append(time)
                TrackCols.append("green") ;TrackRes.append(Re)
    return TrackMem,TrackTimes,TrackCols,TrackRes

## instanciated version of Apply
# @param lmin int minimum sub syndrome size 
# @param Remax int ditto with max
# @param nbls int number of computed points
# @param nbtarget int number of computed targets for memBlowUp
# @return (M,T, col, ls) with M[i] (resp. T[i]) log2 of memory (resp. time) cost 
# col keeps track of the used algorithm and ls the sub syndrome size used 
# In other words M[i] (resp. T[i]) is a memory (resp. time)  cost with algorithm col[i] and l coefficient res[i]
def real_Apply(n, k, w, r, h, lmin, lmax, nbls = 100, nbtarget =80):
    c = expo_complexity_dissection(r)
    step = int((lmax-lmin)/nbls)
    TrackMem = [] ; TrackTimes = [] ; TrackCols =[]; Trackls = []
    
    l = int(lmin)

    while l < lmax:
        Skkl = real_skkl(n,w,k,l)

        if real_is_classic_OK(n,k,l,r,c,h, Skkl):
            (mem,time) = real_apply_classic(n,k,l,r,c,h, Skkl)
            TrackMem.append(mem)
            TrackTimes.append(time)
            TrackCols.append("red"); Trackls.append(l)
        elif real_is_calibrated_OK(n,k,l,r,c,h, Skkl): 
            (mem,time,_) = real_apply_calibration(n,k,l,r,c,h, Skkl,1, asymptotic=False)
            TrackMem.append(mem)
            TrackTimes.append(time)
            TrackCols.append("blue") ;Trackls.append(l)
        elif real_is_smoothed_OK(n,k,l,r,c,h, Skkl):
            (mem,time) = real_apply_smoothed(n,k,l,r,c,h, Skkl)
            TrackMem.append(mem)
            TrackTimes.append(time)
            TrackCols.append("pink") ;Trackls.append(l)
        else:
            (mem,time) = real_apply_mem_blowup(n,k,l,r,c,h, Skkl,nbtarget)
            if mem != -1:
                TrackMem.append(mem)
                TrackTimes.append(time)
                TrackCols.append("green") ;Trackls.append(l)
        l = l + step
        
    return TrackMem,TrackTimes,TrackCols,Trackls
            

    