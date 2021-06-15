## @package syndrome_decoding
# Implements ISD algorithm.
# Calls for the SubSetSum step a 3,4 dissection


## From a given partity check matrix H with associated syndrome s, 
# It transforms H and s to have a partial reduced form
# @returns Hp, sp, perm and S such that
# Hp = S * H and has needed form in PGE+SS framework
# sp = S * s , pi permutation such that
# Hp x e = sp <=> H x pi(e) = s
# Algorithm is randomized : it returns a different solution at each step (though the number of solutions is limited)
def build_H(H, s,n,k,l):
    repeat = True
    H_original = copy(H)
    while (repeat):
        H = copy(H_original)

        # Seek for a good permutation pi such that the sub matrix H11 of  pi * His invertible  
        perm  = Permutations(n).random_element()
        H.permute_columns(perm)
        H11 = H[0:n-k-l, 0:n-k-l]
        while(not(H11.is_invertible())):
            H = copy(H_original)
            perm = Permutations(n).random_element()
            H.permute_columns(perm)
            H11 = H[0:n-k-l, 0:n-k-l]

        # From now, we will work on  H * Id.permute_columns(perm), call it   H *pi
        # Thus if we have e tel such that H * pi * e = s, then the solutions of H * v = s is pi * e
        # But anyway we will continue to transofrm H

        # Now we seek S such that S * H has desired form
        # S must look like that : 
        # H11-1 | 0
        # H12^t | -H12^t H11 H21^t(H21H21^t)-1  noted A et B in the code
        # the H11-1 and 0 are linked to the fact that we want identity on n-k-l
        # the 2 other values are explained by the fact that we want 0 below identity
        # For the two last values we have a 2 equations - 1 unknown : maybe there was a better choice ?

        H21 = H[n-k-l:n-k,0:n-k-l]
        H12 = H[0:n-k-l,n-k-l:n]
        H22 = H[n-k-l:n-k,n-k-l:n]
        H11_invert = ~H11
        B = identity_matrix(F3,l)
        A = -H21*H11_invert


        #A.augment(B) # partie du bas

        top_part = (H11_invert).augment(matrix(F3, nrows = n-k-l, ncols = l))
        bottom_part = A.augment(B)
        S = top_part.stack(bottom_part)

        H = S*H
        


        # Now
        # H_new e = s_new
        # <=> S H pi e = s_new
        # <=> H pi e = S-1 s_new
        # so that we want on veut S-1 s_new = s
        # So  s_new = S s
        s = S * s

        repeat = not(rank(H22) == l) # need full rank matrix

    return(H, s, perm, S)


## Solves the SSS H_ss * e  = sub_s using a four-dissection
# verbose option prints infromations about expected list cardinal & can display warning
#  The algorithm is deterministic, even thought we are supposed to pick in the lists in a random way
def four_tree_list_tree(l,sub_s, H22, verbose = False):

    # Because we are dealing with concrete instances, we can not have same list cardinal everywhere 
    # so we needed to taylor the list caridnal to have exactly what we want
    cardinal_list_6 = 2**6
    cardinal_list_7 = 108

    target_11 = 4
    target_12 = 8


    Lists_level_2 = []
    Lists_level_3 = []
    index_H = 0


    ####################################
    ############# LEVEL 1 ##############
    ####################################
    if (verbose):
        print("\n\n LEVELLL 1\n\n")
        print("Expected length 53.3 for three first ones otherwise / 90 ")
    for i in range(0,4**3,4):
    # Build 4 lists and do a 4-dissection
    # The H[ugly indexes] represent the columns of H upon which we should do linear combinations for each list

    #this is where we need to change list initialisation if the n paraùmeter is changed
        if (i < 4*3): 
            L1 = build_list_matrix(H22[:,index_H : index_H+6], i, sub_s,cardinal_list_6)
            L2 = build_list_matrix(H22[:,index_H+6 : index_H+12], i+1, sub_s,cardinal_list_6)
            L3 = build_list_matrix(H22[:,index_H+12 : index_H+18], i+3, sub_s,cardinal_list_6)
            L4 = build_list_matrix(H22[:,index_H+18 : index_H+18+7], i+4, sub_s,cardinal_list_7)
            index_H = index_H + 18+7
        else:
            L1 = build_list_matrix(H22[:,index_H : index_H+6], i, sub_s,cardinal_list_6)
            L2 = build_list_matrix(H22[:,index_H+6 : index_H+12], i+1, sub_s,cardinal_list_6)
            L3 = build_list_matrix(H22[:,index_H+12 : index_H+12+7], i+3, sub_s,cardinal_list_7)
            L4 = build_list_matrix(H22[:,index_H+12+7 : index_H+12+14], i+4, sub_s,cardinal_list_7)
            index_H = index_H + 12+14


        L1234 = four_dissection(L1,L2,L3,L4,l, 0,target_11,target_12, False, verbose)


        Lists_level_2.append(L1234)
        if verbose:
            print("len list " + str(len(L1234)))

    # Just for debug, to remove
    assert(len(Lists_level_2) == 4**2)

    ####################################.
    ############# LEVEL 2 ##############
    ####################################
    if verbose:
        print("\n\n LEVELLL 2\n\n")
        print("Expected length 25.6 for first one or 122.9")

    # Building the lists and doing 4-dissection
    for i in range(0,4**2,4):
        # Build 4 lists and do a 4-dissection
        L1 = Lists_level_2[i]
        L2 = Lists_level_2[i+1]
        L3 = Lists_level_2[i+2]
        L4 = Lists_level_2[i+3]

        L1234 = four_dissection(L1,L2,L3,L4,l,target_12+target_11, target_11,target_12, False, verbose)

        #check_syndrome = (i==0)
        #oracle_test_four_list_tree(Htmp,L1234,s, 0,(target_11+target_12)*2, check_syndrome)

        if verbose:
            print("len list " + str(len(L1234)))

        Lists_level_3.append(L1234)

    # Just for debug, to remove  
    assert(len(Lists_level_3) == 4)

    ####################################
    ############# LEVEL 3 ##############
    ####################################
    if verbose:
        print("\n\n LEVELLL 3\n\n")
        print("Expected length 804 ")

    L1 = Lists_level_3[0]
    L2 = Lists_level_3[1]
    L3 = Lists_level_3[2]
    L4 = Lists_level_3[3]

    L = four_dissection(L1,L2,L3,L4,l, (target_12+target_11)*2, target_11, l - (target_12+target_11)*2 +- target_11 , True, verbose)
    List_L = L[matrix(F3, nrows=l,ncols=1,immutable = True)]
    if verbose:
        print("Number of solutions " + str(len(List_L)))
        print("end")
        print("Theoretical expected number of solutions: 804 ")
    return List_L




## Main implementation of SDP
# @param H partity check matrix
# @param s syndrome target
# @param weight weight taregt
# @param n,k,l usual code parameters
def SDP(H,s,n,k,l,weight):

    nb_iter = 0
    perm = Permutations(n).identity()


    target_weight =  weight-k-l #for probabilistic step
    print("target weight " + str(target_weight))
    print("")
    condition = True
    while(condition): #check about rank of instanced
        print("One iteration !")
        nb_iter = nb_iter + 1

        # Transform H such that it has desired form
        (Hp,sp,perm,S) = build_H(H, s,n,k,l)

        H22 = Hp[n-k-l:n-k,n-k-l:n]
        H12 = Hp[0:n-k-l,n-k-l:n]
        sub_s = matrix(sp[n-k-l:n-k], ncols = 1) 
        
        # Do the 3,4-dissection
        List_L = four_tree_list_tree(l,sub_s, H22, verbose = False)
                
        # check if target hamming weight is reached somewhere 
        hamming_weights = [ (sp[0:n-k-l]- vector(H12*e_bottom)).hamming_weight() for e_bottom in List_L ]
        print("Here is max hamming weight I found :  " +str(max(hamming_weights)))
        if(target_weight in hamming_weights): # we found a solution
            condition = False
            hamming_weights = np.array(hamming_weights)
            bool_hamming_weights = (hamming_weights == target_weight)

            e_bottom = matrix(np.array(List_L)[bool_hamming_weights][0], ncols = 1)
            e_top = matrix(sp[0:n-k-l], ncols = 1) - H12*e_bottom
            e = e_top.stack(e_bottom)
            PI = identity_matrix(F3,n)
            PI.permute_columns(perm)

            """ print("We have a solution ! Let's just test if the solutions are well formed")


            for e_bottom in List_L:        
                e_top = matrix(sp[0:n-k-l], ncols = 1) - H12*e_bottom
                e = e_top.stack(e_bottom)
                if (H*PI*e != matrix(s, ncols = 1)):
                    print("Error !!!!!!!!!  merging list does not work")
                    print(vector(s))
                    print("")
                    print(vector(H*PI*e))
                    print("")
                    print(vector(H*e))
                    break"""
            


            print(" \n\n nb iteration "+ str(nb_iter) + "\n\n")
            return PI*e
