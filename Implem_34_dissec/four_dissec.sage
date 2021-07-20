## @package four_dissec
# Implementation of 4-dissection


# A few comments about matrix version :
# We no longer deal with vectors, but with matrices with ncol = 1
# Why ? because a priori there is no way to concatenate vectors unless if we go through lists
# But according to some experiments, lists are really expensive




## Performs a 4-dissection 
# The target is assumed to be zero, so we must have syndrome already included in the lists, for example replace L1 by L1 - s
# @param L_i a dict with key:\f$H_ie_i in F_3^l \f$ -> value:\f$e_i\f$
# @param l the sub-syndrome size
# @param start_target_point, target1, target2 are such that:
# * The merge is done on on target1+taregt2 coordinates
# * The merge starts at index start_target_point
# * The first level of the dissection does a merge on target1 coordinates. In particular, it iterates over all intermediate targets (so on a space of search size \f$3^{target1]\f$)
# @param lastLevel : shall be True at the last dissection done. When set to False, the what is returned is a dict of form \f$H_ie_i in F_3^l \f$ -> value:\f$e_i\f$ (so it is useful for the first two levels of the 3,4 dissection)
# @param target_weight,H12,sp only set if lastLevel=True: it performs a hamming weight check on-the-fly
# In this case, if there are two distinct vectors e1 and e2 such that He1 = He2, then one of the two solutions is thrown away. With l = 34, this happens with probability 1/3**24 for the first level and 1/3**12 for the second level, so it is rare enough to be ignored here.
# However, at the last level, a check on-the-fly is performed so that we don't need to store the solutions in a list
# @attention vectors must not be given under the form of matrices 
def four_dissection(L1,L2,L3,L4, l, start_target_point, target1, target2, lastLevel=False, target_weight=0, H12=0, sp=0, verbose=True):
    if lastLevel ==False:
        L_final = {}
    else:
        print(target_weight)



    # We nullify elements starting form start_target_point up to min(start_target_point+target1+target2,l)
    end_target_point = min(l, start_target_point+target1+target2)

    ctarget_vector = matrix(F3, nrows = target1, ncols = 1)
    ctarget_vector_top_pad = matrix(F3, nrows = start_target_point, ncols = 1)
    ctarget_vector_bottom_pad = matrix(F3, nrows = l-target1-start_target_point, ncols = 1)
    
    #a 4-dissection is exhaustive, so we iterate over all intermediate targets
    for i in range(0, 3**target1):
        
        # We want L_1 + L_2 = -ctarget, L_3 + L_4 = ctarget
        # So -L_1 - ctarget = L2 / ctarget - L_3 = L4 
        # Therefore we create L1_new (resp. L3_new) which constains the left handside of the first (resp. second) equation
        # As we need to find partial collisions (and a priori a dict only allows for full collisions), we need to build a new dictionnary which maps the partial targets to the full taregts
        # This is what does the Pool_Li variables
        Pool_L1 = {}
        Pool_L3 = {}
        
        L1_new={}
        L3_new ={}


        full_vector_ctarget = ctarget_vector_top_pad.stack(ctarget_vector.stack(ctarget_vector_bottom_pad))

        # build L1_new and Pool_L1
        for key_L1 in L1:
            vect = L1[key_L1]
            new_key_L1 = -key_L1  - full_vector_ctarget
            key_hashed = new_key_L1[start_target_point:start_target_point+target1]
            key_hashed.set_immutable() # keys of dict must be immutable
            new_key_L1.set_immutable()
            if key_hashed in Pool_L1:
                Pool_L1[key_hashed] = Pool_L1[key_hashed] + [new_key_L1]
            else:
                Pool_L1[key_hashed] = [new_key_L1]
                
            L1_new[new_key_L1] = -vect
            
        # Do the same with L3
        for key_L3 in L3:
            vect = L3[key_L3]
            new_key_L3 = -key_L3  + full_vector_ctarget
            key_hashed = new_key_L3[start_target_point:start_target_point+target1]
            key_hashed.set_immutable()
            new_key_L3.set_immutable()
            
            if key_hashed in Pool_L3:
                Pool_L3[key_hashed] = Pool_L3[key_hashed] + [new_key_L3]
            else:
                Pool_L3[key_hashed] = [new_key_L3]
                
                
            L3_new[new_key_L3] = -vect
            
               
        

        # For each element in L2, hash its key, go in the pool to see if the key appears
        # if key_hashed appears in the pool : 
        # 1) retrieve a = -L1 - c which is in the pool and corresponding vector
        # 2) retrieve b = L2 which is in L2
        # 3) Store in L12 b-a = L1 + L2 + c associated with vector concatenated
        # Search for all element in L2 partial collisions with elements in Pool_L1, retrieve the full key "new_key", store in Pool_L12 the mapping between hashed key (a slice of the whole key - so that it is ready for the next merge) and the full key, 
        # And store in L12 the mapping full key (H_ie_i) -> vector (e_i)
        Pool_L12 = {}
        L12 = {}
        for key_L2 in L2:
            key_hashed = key_L2[start_target_point:start_target_point+target1]
            key_hashed.set_immutable()
            if key_hashed in Pool_L1:
                list_keys_L1 = Pool_L1[key_hashed]
                for key_L1 in list_keys_L1:
                    new_key = (key_L2 - key_L1)
                    new_key.set_immutable()
                    L12[new_key] = (-L1_new[key_L1]).stack(L2[key_L2])
                    key_hashed2 = new_key[target1+start_target_point:end_target_point]
                    key_hashed2.set_immutable()
                    if key_hashed2 in Pool_L12:
                        Pool_L12[key_hashed2] = Pool_L12[key_hashed2] + [new_key]
                    else:
                        Pool_L12[key_hashed2] = [new_key]
                   
  
                

        # There is no need to build Pool_L34 and L34 because collisions can be fully checked on-the-fly   
        # We do the following: 
        # for each l4 in L4:
        ## Check for collisions in L3, for each found collision:
        ### Retrieve the key new_key_L34
        ### Check for collisions in Pool_L12, for each found collision
        #### Retrieve the key and the vector: we have a candidate
        for key_L4 in L4:
            key_hashed = key_L4[start_target_point:start_target_point+target1]
            key_hashed.set_immutable()
            if key_hashed in Pool_L3:
                list_keys_L3 = Pool_L3[key_hashed]
                for key_L3 in list_keys_L3:
                    new_vector = (L3_new[key_L3]).stack(-L4[key_L4])
                    new_key_L34 = key_L3 - key_L4
                    new_key_L34.set_immutable()
                    key_hashed = new_key_L34[target1+start_target_point:end_target_point]
                    key_hashed.set_immutable()
                    if key_hashed in Pool_L12:
                        list_keys_L12 = Pool_L12[key_hashed]
                        for key_L12 in list_keys_L12:
                            #new_vector_final  = vector(L12[a].list() + (-new_vector).list() )
                            new_vector_final  = L12[key_L12].stack(-new_vector)
                            new_key_final = -new_key_L34 + key_L12
                            new_key_final.set_immutable()
                            if lastLevel ==False:
                                L_final[new_key_final] = new_vector_final
                                # if new_key_final in L_final means that we throw away a solution
                            else:
                                # check-on-the-fly
                                if (sp[0:n-k-l]- vector(H12*new_vector_final)).hamming_weight() == target_weight:
                                    return new_vector_final




        # update the target vector
        # see it at adding 1 to ctarget (addition with carry)
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

    if lastLevel ==False:
        return L_final
    else:
        return 0


