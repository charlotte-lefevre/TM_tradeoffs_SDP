## @package dissection_optimize
# Useful to compare in the case of Wave the dissection technique with one other method: it optimizes and select the best points

import numpy as np
from dissection import dissection
from utils import clean_array
from Wave_param import * # R and W

## @brief computes empirically the best time-memory trade-offs reachable by the dissection with Wave parameters and given number of lists
# @param final_treatement boolean default true. True iff we just want to keep the best points
# @return all used Re's coefficients (relative Subset Sum target size) \n
# Note that this function does not include a precomputation phase, which consits in "looking roughly at the best Re's candidates". 
# The code here is just a precise way to get the best possible points
def compute_best_complexity_dissection(r=120, final_treatement=True): 
    M_tot = [];T_tot = [];cols_tot = []

    #Res = np.append(np.linspace(0.001, 0.025, 25), np.linspace(0.025, 0.034, 10))
    Res = np.linspace(0.001,0.034, 45)


    for Re in Res:
        (M,T,_) = dissection(r, R, Re,W, nbMem=100)
        #plt.scatter(M,T, label="l" + r'$ \approx$' + str(round(Re,4)) + "n", s=20)
        #plt.scatter(M,T, label="l = " + str(round(Re,4)) + "n", s=20)
        #plt.scatter(M,T, s=20)
        M_tot = np.append(M_tot, M)
        T_tot = np.append(T_tot, T)
    if final_treatement:
        M_tot, T_tot = clean_array(M_tot,T_tot, [])

    return M_tot,T_tot

