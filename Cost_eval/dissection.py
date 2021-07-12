
## @package dissection
# Implementation of cost computation of time-memory trade offs available with dissection


import math
import numpy as np

from utils import *

## Useful to plot the different regimes for the dissection. see the dissection definition for more informations
cdark = ["navy", "forestgreen", "purple", "darkcyan", "chocolate", "firebrick", "darkgreen", "darkgoldenrod", "indigo", "darkcyan"]
clight = ["blue", "lightgreen", "violet", "cyan", "darkorange", "coral", "lightgreen", "gold", "blueviolet", "cyan"]


##  @brief computes cost of dissection for the SubSetSum algorithm using stirling approximation
# @param nbMem int default to 100. Fixed the number of potentially computed points for one fixed m.
# @param optimize boolean. default is set to false. If you want to speedup the time it takes and do not mind to loose slightly in precision, set it to true. \n
# for more informations, go to section colors
# @return (m,t,colors) tuple of arrays, with m containing memories, t times colors colors, 
# \f$\forall\f$ i, one point corresponds to a dissection instanciation with memory cost \f$M = 2^{nm}\f$, time cost \f$MT= 2^{nt}\f$ with n the size of the code
# and color[i] is associated to t[i] and m[i]. Informations about the colors of the points  :
# * red : we are in a normal regime
# * black : the granularity of the dissection is not fine enough for us to have solutions in O(1). 
# To stay in this regime, it was necessary to balance the number of required solutions in the upper dissection and lower dissection
# * colored points : doing the previous balancing does not work. Here, we desesperately try to reduce the size of the meet-in the middle.
# @attention This function has been desgined such that maximum allowed  memory is the one of Wave state of the art algorithm. 
# If you are not working with Wave parameters, you should change max_m value
def dissection(r, R, Re,W, nbMem=100): 
    Skkl = skkl(Re,R, W)
    M = [];T = [];color=[]

    # BaseMem is such that the leaf cardinals is equal to 3^(BaseMem x ncode)
    max_BaseMem = math.log(2,3)*(R+Re)/r
    min_BaseMem = (Skkl+Re)/r    

    for BaseMem in np.linspace(min_BaseMem, max_BaseMem, nbMem):
        # ensures we do not have Bigger Memory than State of the art algorithm (because else it will be useless)
        max_m = int(math.log(2,3)*0.0176/BaseMem) #state of the art memory
        for m in range(1,max_m+1): # memory  parameter
            g = gain(r,m)
            if BaseMem<=Re/(g+m): # to avoid 'useless' points, but can be removed 
                a = 1/BaseMem * Skkl +max(0,- g + Re/BaseMem - m) #represents desired number of solutions on the recursive dissection (need 3^{BaseMem * a}) solutions 
                gain_gplusm_dissection= gain(g+m,m)
                gain_g_dissection= gain(g,m)
                # if complexity of g+m dissection is dominating, we enter the dissection to reduce the amortized time
                if (a < max(m,(g-gain_gplusm_dissection) )):                  
                    if ( a >= max(g - 2*m - 2*gain_g_dissection+gain_gplusm_dissection, m)): # can adjust the number of solutions returned by small dissection
                        if Skkl<=BaseMem*(a+g-gain_gplusm_dissection)/2:
                            M = np.append(M,m*BaseMem*math.log(3,2))
                            T = np.append(T,max(Skkl, BaseMem*(a+g-gain_gplusm_dissection)/2)*math.log(3,2))
                            color = np.append(color, "black")

                    elif (a >m): # try to increase the sub dissection size 
                        u = g
                        condition = True
                        while(u >=2 and condition):
                            u = u - 1
                            a = a + 1
                            if (a >= u - gain(u+m,m)): # no granularity problem
                                M = np.append(M,m*BaseMem*math.log(3,2))
                                T = np.append(T, max(Skkl, BaseMem*a)*math.log(3,2))
                                color = np.append(color, cdark[(m-1) %10])
                                condition = False
                            elif(a >= max(u - m - gain(u,m), m)): # can adjust the number of solutions returned by small dissection
                                M = np.append(M,m*BaseMem*math.log(3,2))
                                T = np.append(T,max(Skkl, BaseMem*(a+u-gain(u+m,m))/2)*math.log(3,2))
                                color = np.append(color, clight[(m-1)  %10])
                                condition = False
                else:
                    M = np.append(M,m*BaseMem*math.log(3,2))
                    T = np.append(T,max(a*BaseMem, Skkl)*math.log(3,2))
                    color = np.append(color, "red")
                        
    return (M, T, color)
