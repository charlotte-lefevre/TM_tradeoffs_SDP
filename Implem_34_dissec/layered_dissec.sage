## prints expected lists cardinals for these values of n,k,l
# When adding a new n parameters, one new case shall be considered
def print_expected_list_cardinals(n,k,l,level,slices_merges):
    if n == 560  and k == 379 and l == 34:
        cardinal_list_6 = 2**6
        cardinal_list_7 = 108
        cardinal_1_3 = (cardinal_list_6**3*cardinal_list_7)/(3**(slices_merges[0]+slices_merges[1]))
        cardinal_1_13 = (cardinal_list_6**2*cardinal_list_7**2)/(3**(slices_merges[0]+slices_merges[1]))
        cardinal_2_1 = (cardinal_1_3**3*cardinal_1_13)/(3**(slices_merges[2]+slices_merges[3]))
        cardinal_2_3 = (cardinal_1_13**4)/(3**(slices_merges[2]+slices_merges[3]))
        nb_sol = (cardinal_2_1*cardinal_2_3**3)/(3**(slices_merges[4]+slices_merges[5]))
        if level == 1:
            print("Expected list cardinal for the 3 first lists: " +str(round(cardinal_1_3,1)) + ", for the 13 others: " + str(round(cardinal_1_13,1)))
        elif level == 2:
            print("Expected list cardinal for the first list: " +str(round(cardinal_2_1,1) )+ ", for the 3 others: " + str(round(cardinal_2_3,1)))
        elif level == 3:
            print("Expected number of returned solutions " + str(round(nb_sol,1)))
        else:
            print("Error in print_expected_list_cardinals: 'level' parameter incorrect")
    elif n == 700 and k ==473 and l ==40:
        cardinal_list = int(3**(l/8))
        cardinal_1 = cardinal_list**4/(3**(slices_merges[0]+slices_merges[1]))
        cardinal_2 = (cardinal_1**4)/(3**(slices_merges[0]+slices_merges[1]))
        nb_sol = cardinal_2**4/(3**(slices_merges[2]+slices_merges[3]))
        if level == 1:
            print("Expected list cardinal: "+str(round(cardinal_1,1)))
        elif level == 2:
            print("Expected list cardinal: "+str(round(cardinal_2,1)))
        elif level == 3:
            print("Expected number of returned solutions " + str(round(nb_sol,1)))
        else:
            print("Error in print_expected_list_cardinals: 'level' parameter incorrect")
    elif n==875 and k==591 and l==48:
        cardinal_list = int(3**(l/8))
        cardinal_list_outlier = 2**9
        cardinal_1_1 = (cardinal_list**3*cardinal_list_outlier)/(3**(slices_merges[0]+slices_merges[1]))
        cardinal_1_15 = (cardinal_list**4)/(3**(slices_merges[0]+slices_merges[1]))
        cardinal_2_1 = (cardinal_1_1*cardinal_1_15**3)/(3**(slices_merges[2]+slices_merges[3]))
        cardinal_2_3 = (cardinal_1_15**4)/(3**(slices_merges[2]+slices_merges[3]))
        nb_sol = (cardinal_2_1*cardinal_2_3**3)/(3**(slices_merges[4]+slices_merges[5]))
        if level== 1:
            print("Expected list cardinal for the first list: " +str(round(cardinal_1_1,1)) + ", for the 15 others: " + str(round(cardinal_1_15,1)))
        elif level == 2:
            print("Expected list cardinal for the first list: " +str(round(cardinal_2_1,1)) + ", for the 3 others: " + str(round(cardinal_2_3,1)))
        elif level== 3:
            print("Expected number of returned solutions " + str(round(nb_sol,1)))
        else:
            print("Error in print_expected_list_cardinals: 'level' parameter incorrect")
    else:
        print("I don't know the theoretical expected cardinal for these particular n,k,l parameters")

## Solves the Subset Sum step  \f$H_s \times e  = sp\f$ using a four-dissection
# @param n,k,l: usual code parameters
# @param H22, H12, sp sub-matrices and syndrome
# @param target_w: target weight of the completion of the solution
# @param times a dict which keeps a trace of times spend in each part of the 3,4 dissec
# @param verbose prints lists cardinal at each step
# @returns e solution of SDP or 0 if no candidate has the desired hamming weight
def layered_dissec(n,k,l, H22, H12, sp, target_w,times,verbose = True):
    sp = matrix(sp, ncols=1)
    sub_s = matrix(sp[n-k-l:n-k], ncols =1)

    
    # tells what is the size of intermediate targets 
    # for example slices_merges[0] + slices_merges[1] is the target size of the first 4-dissection
    # slices_merges[0] is the size of the intermediat etarget within the 4-dissection
    # slices_merges[1] + slices_merges[2] is the target size of the second  4-dissection and so on 
    # when n is changed this shall be changed 
    
    # if n = 560
    if n == 560:
        slices_merges = [4,8,4,8,4,6]
    elif n == 700:
        slices_merges = [5,10,5,10,5,5]
    elif n == 875:
        slices_merges = [6,12,6,12,6,6]
    else:
        print("Error: this n value has not been implemented. You shall modify 'slices_merges' and 'build_lists_division' method, because of rounding problems we do not it automatically. \n I exit.")
        exit()

    # if n = 700 
    #slices_merges = [5,10,5,10,5,5]

    # if n = 875
    #slices_merges = [6,12,6,12,6,6]
    
    len_already_merged = 0

    Lists_level_2 = []
    Lists_level_3 = []

    # keeps track of which columns of H have already been explored when building the leaves lists
    index_H = 0

    array_ai = []


    ####################################
    ############# LEVEL 1 ##############
    ####################################
    if (verbose):
        print("\n\n LEVEL 1\n\n")
        print_expected_list_cardinals(n,k,l,1,slices_merges)

    # the ais represent the intermediate targets. Indeed, insead of searching for collisions of form L_{i} - L_{i+1} = 0, 
    # we search for collisions of form L_{i} - L{i+1} = a_i
    # and we need that a_1+a_2+..+a_16 = 0
    # by doing this we avoid to re-do a partial gaussian elimination when one iteration of the 3,4 dissection is not enough
    null_vect = matrix( F3,nrows=l, ncols = 1)
    sum_ais = matrix(F3, ncols = 1, nrows = l) 
    for i in range(0,4**3,4):
    # Build 4 lists and do a 4-dissection
    # The H[indexes] represent the columns of H upon which we should do linear combinations for each list
    #If n is changed, we need to change this step
        L1,L2,L3,L4, index_H,sum_ais =  build_lists_division(i, index_H,sum_ais,array_ai,H22, sub_s, slices_merges, n,k,l,times)
        # oracle_test_integrity enable to test if the method build_lists_division was well done, uncomment this part to use the test 
        """if i == 0:
            oracle_test_integrity(L1, H22, "L1 level 0", True, sub_s-array_ai[int(i/4)], include_s = True)
        else:
            oracle_test_integrity(L1, H22, "L1 level 0 i: "+str(i), True, -array_ai[int(i/4)], include_s = True)
        oracle_test_integrity(L2, H22, "L2 level 0", False)
        oracle_test_integrity(L3, H22, "L3 level 0", True)
        oracle_test_integrity(L4, H22, "L4 level 0", False)"""
        beg = time.time()
        L1234 = four_dissection(L1,L2,L3,L4,l, len_already_merged, slices_merges[0], slices_merges[1], slices_merges[2], (i/4 % 2 == 0), False,times) 
        end = time.time()
        times["Four dissection (partionned)"]  = times["Four dissection (partionned)"]  + end - beg
        if i == 0:
            oracle_test_integrity(L1234, H22, "L1234 level 1", (i/4 % 2 == 0), sub_s-array_ai[int(i/4)], include_s = True)
        else:
            oracle_test_integrity(L1234, H22, "L1234 level 1", (i/4 % 2 == 0), -array_ai[int(i/4)], include_s = True)

        Lists_level_2.append(L1234)
        if verbose:
            print("len list " + str(len(L1234)))

    len_already_merged = len_already_merged + slices_merges[0] + slices_merges[1]
    ####################################.
    ############# LEVEL 2 ##############
    ####################################
    if verbose:
        print("\n\n LEVEL 2\n\n")
        print_expected_list_cardinals(n,k,l,2,slices_merges)
    # doing 4-dissections

    for i in range(0,4**2,4):
        # Build 4 lists and do a 4-dissection
        beg = time.time()
        L1234 = four_dissection(Lists_level_2[i],Lists_level_2[i+1],Lists_level_2[i+2],Lists_level_2[i+3],l,len_already_merged,slices_merges[2], slices_merges[3], slices_merges[4], (i/4 % 2 == 0), False,times)
        end = time.time()
        times["Four dissection (partionned)"]  = times["Four dissection (partionned)"]  + end - beg

        if verbose:
            print("len list " + str(len(L1234)))

        Lists_level_3.append(L1234)

    len_already_merged = len_already_merged + slices_merges[2] + slices_merges[3]
    ####################################
    ############# LEVEL 3 ##############
    ####################################
    if verbose:
        print("\n\n LEVEL 3\n\n")
        print_expected_list_cardinals(n,k,l,3,slices_merges)

    beg = time.time()
    e = four_dissection( Lists_level_3[0], Lists_level_3[1], Lists_level_3[2], Lists_level_3[3],l, len_already_merged, slices_merges[4], slices_merges[5], 0, True, True,times,sp[0:n-k-l,0], H12, target_w) 
    end = time.time()
    times["Four dissection (partionned)"]  = times["Four dissection (partionned)"]  + end - beg
    return e