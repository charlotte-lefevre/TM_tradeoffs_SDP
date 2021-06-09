## @package dissection_optimize
# Useful to compare in the case of Wave the dissection technique with one other method: it optimizes and select the best points

import numpy as np
from dissection import dissection
from utils import clean_array
from Wave_param import * # R and W

## @brief computes empirically the best time-memory plot for a dissection with Wave parameters and given number of lists
# @param r int the dissection order
# @param final_treatement boolean default true. True iff we just want to keep the best points
# @return all used Re's coefficients (relative subsetSum target size) \n
# Note that there was already a precomputation phase, which was rather "looking roughly at the best Re's candidates". 
# The code here is just a precise way to get the best possible points
def compute_best_complexity_dissection(r=120, final_treatement=True): 
    M_tot = [];T_tot = [];final_Res = []
    Res = np.linspace(0.025, 0.037, 30)
    for Re in Res:
        (M,T,col) = dissection(r, R, Re,W)
        for i in range(0, len(M)):
            if col[i] != "red":
                M_tot = np.append(M_tot, M[i])
                T_tot = np.append(T_tot, T[i])
                final_Res = np.append(final_Res, Re)
                
        if (final_treatement):
            M_tot, T_tot, final_Res = clean_array(M_tot,T_tot, final_Res)

    return(M_tot, T_tot,final_Res )


"""def compute_best_complexity_dissection_real(n,k, w, r=120, traitement_final=True,nb_points = 4): 

    M_tot = [];T_tot = [];col_tot = [];final_ls = []
    steps =  int((int(0.037*n) - int(0.025*n))/nb_points)
    ls = range(int(0.025*n), int(0.037*n),steps)
    for l in ls:
        print(l)
        (M,T,col) = dissection_real(r,n, k, l,w)
        M_tot = np.append(M_tot, M)
        T_tot = np.append(T_tot, T)
        col_tot = np.append(col_tot, col)
        final_ls = np.append(final_ls, [l for _ in range(0,len(M)) ])
        #plt.plot(M, T,markersize = 3, linewidth = 0, marker = '.', )     
        if (traitement_final):
            M_tot, T_tot, final_ls = traitement_array(M_tot,T_tot, final_ls)
    return(M_tot, T_tot,final_ls )    """
    