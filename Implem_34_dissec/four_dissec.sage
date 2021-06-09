## @package four_dissec
# Impelmentation of 4 - dissection


# A few comments about matrix version :
# We no longer deal with vectors, but with matrices with ncol = 1
# Why ? because a priori there is no way to concatenate vectors unless if we go through lists
# But according to some experiments, lists are really expensive


## Performs a 4-dissection 
# The target is assumed to be zero, so we must have syndrome already included in the lists
# @param start_target_point  tells in which index do we start the lookup
# @param L_i a dict with key->value as linear_combin(vectors to merge) -> underlying vector  : hash values : (Valeur CL target a s, e1)
# @param target1, target2 are the size of the targets for level 1 and level 2. Normally, target1+target2 = 3*x, target1 = x target2 = 2x
# @param return_as_list : if we know that the size of the linear combinations is very big compared to the expected number of solutions, then set it as false, so that a dictionnary will be returned
# Else, if for example this is the last step of the dissection, set is as True, so that a list will be returned
# @attention vectors must not be given under the form of matrice sfor practical reasons
def four_dissection(L1,L2,L3,L4, l, start_target_point, target1, target2, return_as_list=False, verbose=True):
    L_final = {}

    nb_Ls = []; i = 0
    nb_L12 = []
    nb_L34 =[]
    """ncol_H = H.ncols()
    H1 = H[:,0:int(ncol_H/2)]
    H2 = H[:,int(ncol_H/2):]"""
    
    # We will nullify bits starting form start_target_point up to max(l, 3*(level)*size_targets)
    end_target_point = min(l, start_target_point+target1+target2)

    ctarget_vector = matrix(F3, nrows = target1, ncols = 1)
    ctarget_vector_top_pad = matrix(F3, nrows = start_target_point, ncols = 1)
    ctarget_vector_bottom_pad = matrix(F3, nrows = l-target1-start_target_point, ncols = 1)
    
    #for ctarget in return_all_vectors(size_targets):
    for i in range(0, 3**target1):
        nb_Ls = nb_Ls + [0]
        nb_L12 = nb_L12 + [0]
        nb_L34 = nb_L34 + [0]
        # We need to rebuild the lists
        L12 = {}
        # We want L_1 + L_2 = -c, L_3 + L_4 = c
        # So -L_1 - c = L2 / c - L_3 = L4 :: we will create pools that enumerate indexes on the left part
        # Pool Lij is of form : keys = key on size_targets values/ value : key
        Pool_L12 = {}
        Pool_L34 = {}
        
        L1_new={}
        L3_new ={}
        
        
        full_vector_ctarget = ctarget_vector_top_pad.stack(ctarget_vector.stack(ctarget_vector_bottom_pad))

        #Play with mutability
        for key_lc in L1: #key_lc : key_linearcombination
            vect = L1[key_lc]
            new_key_L1 = -key_lc  - full_vector_ctarget
            key_hashed = new_key_L1[start_target_point:start_target_point+target1]
            key_hashed.set_immutable()
            new_key_L1.set_immutable()
            if key_hashed in Pool_L12:
                Pool_L12[key_hashed] = Pool_L12[key_hashed] + [new_key_L1]
            else:
                Pool_L12[key_hashed] = [new_key_L1]
                
            L1_new[new_key_L1] = -vect
            
        for key_lc in L3:
            vect = L3[key_lc]
            new_key_L3 = -key_lc  + full_vector_ctarget
            key_hashed = new_key_L3[start_target_point:start_target_point+target1]
            key_hashed.set_immutable()
            new_key_L3.set_immutable()
            
            if key_hashed in Pool_L34:
                Pool_L34[key_hashed] = Pool_L34[key_hashed] + [new_key_L3]
            else:
                Pool_L34[key_hashed] = [new_key_L3]
                
                
            L3_new[new_key_L3] = -vect
            
               
        #
        # Layer 1
        #
        
        Pool_L1234 = {}
        
        # For each element in L2, hash its key, go in the pool to see if the key appears
        # if key_hashed appears in the pool : 
        # 1) retrieve a = -L1 - c which is in the pool and corresponding vector
        # 2) retrieve b = L2 which is in L2
        # 3) Store in L12 b-a = L1 + L2 + c associated with vector concatenated
        
        for key_lc in L2:
            key_hashed = key_lc[start_target_point:start_target_point+target1]
            key_hashed.set_immutable()
            if key_hashed in Pool_L12:
                list_as = Pool_L12[key_hashed]
                for a in list_as:
                    new_key = (key_lc - a)
                    new_key.set_immutable()
                    #L12[new_key] = vector((-L1_new[a]).list() + L2[key_lc].list())
                    L12[new_key] = (-L1_new[a]).stack(L2[key_lc])
                    nb_L12[i] = nb_L12[i] + 1
                    key_hashed2 = new_key[target1+start_target_point:end_target_point]
                    key_hashed2.set_immutable()
                    if key_hashed2 in Pool_L1234:
                        Pool_L1234[key_hashed2] = Pool_L1234[key_hashed2] + [new_key]
                    else:
                        Pool_L1234[key_hashed2] = [new_key]
                   
        # For each element in L4, hash its key, go in the pool to see if the key appears
        # if key_hashed appears in the pool : 
        # 1) retrieve a = -L3 + c which is in the poola nd corresponding vector
        # 2) retrieve b = L4 which is in L4
        # 3) Store in L34 a-b = -L3 - L4 +c associated with vector concatenated
        # 4) Take advantage of the loop to constitude Pool_L1234 
        
                
        #    
        # Layer 2 : Pool_L1234
        #
        
        # Do the same as previously 
        # L12 contains elements of the form L1 + L2 + c
        # Pool_L1234 of the form -L3 - L4 + c
        # So if there is a collision :
        # 1) retrieve a = -L3 - L4 - c which is in the pool and and corresponding vector
        # 2) retrive b = L1 + L2 + c
        # 3) Store in the final list b - a as a key along with new vector     
        for key_lc in L4:
            key_hashed = key_lc[start_target_point:start_target_point+target1]
            key_hashed.set_immutable()
            if key_hashed in Pool_L34:
                list_as = Pool_L34[key_hashed]
                for a in list_as:
                    #new_vector = vector(L3_new[a].list() + (-L4[key_lc]).list())
                    new_vector = (L3_new[a]).stack(-L4[key_lc])
                    new_key = a - key_lc
                    new_key.set_immutable()
                    #L34[new_key] =new_vector
                    nb_L34[i] = nb_L34[i] + 1
                    key_hashed = new_key[target1+start_target_point:end_target_point]
                    key_hashed.set_immutable()
                    if key_hashed in Pool_L1234:
                        list_as = Pool_L1234[key_hashed]
                        for a in list_as:
                            #new_vector_final  = vector(L12[a].list() + (-new_vector).list() )
                            new_vector_final  = L12[a].stack(-new_vector)
                            new_key_final = -new_key + a
                            new_key_final.set_immutable()
                            if return_as_list ==False:
                                if new_key_final in L_final and verbose:
                                    print("WARNING : You choose to return dissection result as a dict that maps ONE key to ONE vector, but because of that we detected that one solution will be simply throwed away. Either this is really really bad luck, either you shall consider setting return_as_list as True, but it will make processing much more complicated. This message will be displayed once for each throwed away solution")

                                else:
                                    L_final[new_key_final] = new_vector_final
                            else:
                                if new_key_final in L_final:
                                    L_final[new_key_final] = [new_vector_final] + L_final[new_key_final]
                                else:
                                    L_final[new_key_final] = [new_vector_final]

                            nb_Ls[i] = nb_Ls[i]+1
                 

        # update the target vector
        if i+1 < 3**target1:
            condition = True
            index = 0
            while(condition):
                if(ctarget_vector[index,0] == F3(2)):
                    ctarget_vector[index,0] = 0
                    index = index + 1
                                   
                else:
                    condition = False
                    ctarget_vector[index,0] = ctarget_vector[index,0] + F3(1)

    return L_final