
## The division of the parity-check matrix depends on the values of n,k,l and each list contains theoretically (k+l)/64 columns of H22, which if oftenly not an integer
# Therefore a manual splitting of H22 must be done, this is the aim of this function
# for example with n==875 and k==591 and l=48, all lists contain 10 columns of H22 except the first one 
def build_lists_division(index_dissec, index_H,sum_ais,array_ai,H22, sub_s, slices_merges, n,k,l,times):
    # Remark : the a_is enable here to consider new runs of the 3,4 dissections when one iteration does not provide enough solutions
    null_vect = matrix( F3,nrows=l, ncols = 1)
    if n == 560  and k == 379 and l == 34:
        cardinal_list_6 = 2**6
        cardinal_list_7 = 108
        if (index_dissec < 4*3): 
            a_i =  random_matrix(F3, ncols = 1, nrows = l)
            array_ai.append(a_i)
            sum_ais = sum_ais - a_i
            if (index_dissec == 0):# need to include the sub-syndrome in the lists
                L1 = build_list_leaves(H22, slices_merges, -sub_s+a_i,cardinal_list_6,index_H,index_H+6,True,times)
                
            else:
                L1 = build_list_leaves(H22, slices_merges, a_i,cardinal_list_6,index_H,index_H+6,True,times)

            L2 = build_list_leaves(H22,slices_merges, null_vect,cardinal_list_6,index_H+6,index_H+12,False,times)
            L3 = build_list_leaves(H22, slices_merges, null_vect,cardinal_list_6,index_H+12,index_H+18,True,times)
            L4 = build_list_leaves(H22, slices_merges, null_vect,cardinal_list_7,index_H+18,index_H+18+7,False,times)
            index_H = index_H + 18+7
        else:

            if(index_dissec == 4**3-4): # means that we are dealing wit a_16, which is fixed because sumof a_i = 0
                L1 = build_list_leaves(H22, slices_merges, sum_ais,cardinal_list_6,index_H,index_H+6,True,times)
                array_ai.append(sum_ais)
            else:
                a_i = random_matrix(F3, ncols = 1, nrows = l)
                array_ai.append(a_i)
                sum_ais = sum_ais - a_i
                L1 = build_list_leaves(H22, slices_merges, a_i,cardinal_list_6,index_H,index_H+6,True,times)


            L2 = build_list_leaves(H22,slices_merges, null_vect,cardinal_list_6,index_H+6,index_H+12, False,times)
            L3 = build_list_leaves(H22, slices_merges, null_vect,cardinal_list_7,index_H+12,index_H+12+7, True,times)
            L4 = build_list_leaves(H22, slices_merges, null_vect,cardinal_list_7,index_H+12+7, index_H+12+14, False,times)
            index_H = index_H + 12+14
    elif n == 700 and k ==473 and l ==40:
        cardinal_list = int(3**(l/8))
        coeff_mult = 8
        coeff_mult_outlier = 9
        if index_dissec ==0:
            # Build 4 lists and do a 4-dissection
            # The H[indexes] represent the columns of H upon which we should do linear combinations for each list
            #If n is changed, we need to change this step
            a_i =  random_matrix(F3, ncols = 1, nrows = l)
            array_ai.append(a_i)
            sum_ais = sum_ais - a_i
            L1 = build_list_leaves(H22, slices_merges, -sub_s+a_i,cardinal_list,index_H,index_H+coeff_mult_outlier,True,times)
            L2 = build_list_leaves(H22,slices_merges, null_vect,cardinal_list,index_H+coeff_mult_outlier,index_H+coeff_mult_outlier+coeff_mult,False,times)
            L3 = build_list_leaves(H22, slices_merges, null_vect,cardinal_list,index_H+coeff_mult_outlier+coeff_mult,index_H+coeff_mult_outlier+2*coeff_mult,True,times)
            L4 = build_list_leaves(H22, slices_merges, null_vect,cardinal_list,index_H+coeff_mult_outlier+2*coeff_mult,index_H+coeff_mult_outlier+3*coeff_mult,False,times)
            index_H = index_H +coeff_mult*3+coeff_mult_outlier

        else:
            if(index_dissec == 4**3-4):
                L1 = build_list_leaves(H22, slices_merges, sum_ais,cardinal_list,index_H,index_H+coeff_mult,True,times)
                array_ai.append(sum_ais)
            else:
                a_i = random_matrix(F3, ncols = 1, nrows = l)
                array_ai.append(a_i)
                sum_ais = sum_ais - a_i
                L1 = build_list_leaves(H22, slices_merges, a_i,cardinal_list,index_H,index_H+coeff_mult,True,times)
            L2 = build_list_leaves(H22,slices_merges, null_vect,cardinal_list,index_H+coeff_mult,index_H+2*coeff_mult, False,times)
            L3 = build_list_leaves(H22, slices_merges, null_vect,cardinal_list,index_H+coeff_mult*2,index_H+coeff_mult*3, True,times)
            L4 = build_list_leaves(H22, slices_merges, null_vect,cardinal_list,index_H+coeff_mult*3, index_H+coeff_mult*4, False,times)
            index_H = index_H + 4*coeff_mult
    elif n==875 and k==591 and l==48:
        #print("debug before " + str(index_H))
        coeff_mult = 10
        coeff_mult_outlier = 9
        cardinal_list = int(3**(l/8))
        cardinal_list_outlier = 2**9

        if index_dissec == 0:
            a_i = random_matrix(F3, ncols = 1, nrows = l)
            array_ai.append(a_i)
            sum_ais = sum_ais - a_i
            L1 = build_list_leaves(H22, slices_merges,-sub_s+ a_i,cardinal_list_outlier,index_H,index_H+coeff_mult_outlier,True,times)
            L2 = build_list_leaves(H22,slices_merges, null_vect,cardinal_list,index_H+coeff_mult_outlier,index_H+coeff_mult+coeff_mult_outlier, False,times)
            L3 = build_list_leaves(H22, slices_merges, null_vect,cardinal_list,index_H+coeff_mult_outlier+coeff_mult,index_H+coeff_mult*2+coeff_mult_outlier, True,times)
            L4 = build_list_leaves(H22, slices_merges, null_vect,cardinal_list,index_H+coeff_mult*2+coeff_mult_outlier, index_H+coeff_mult*3+coeff_mult_outlier, False,times)
            index_H = index_H + 3*coeff_mult+coeff_mult_outlier

        else:
            if(index_dissec == 4**3-4):
                L1 = build_list_leaves(H22, slices_merges, sum_ais,cardinal_list,index_H,index_H+coeff_mult,True,times)
                array_ai.append(sum_ais)
            else:
                a_i = random_matrix(F3, ncols = 1, nrows = l)
                array_ai.append(a_i)
                sum_ais = sum_ais - a_i
                L1 = build_list_leaves(H22, slices_merges, a_i,cardinal_list,index_H,index_H+coeff_mult,True,times)
            L2 = build_list_leaves(H22,slices_merges, null_vect,cardinal_list,index_H+coeff_mult,index_H+2*coeff_mult, False,times)
            L3 = build_list_leaves(H22, slices_merges, null_vect,cardinal_list,index_H+coeff_mult*2,index_H+coeff_mult*3, True,times)
            L4 = build_list_leaves(H22, slices_merges, null_vect,cardinal_list,index_H+coeff_mult*3, index_H+coeff_mult*4, False,times)
            index_H = index_H + 4*coeff_mult
    


    else:
        print("This instanciation with this n is not implemented ! you should modify the function 'build_lists_division' and the variable 'slices_merge'")
        exit()

    #print("debug after " + str(index_H))
    return L1,L2,L3,L4, index_H,sum_ais


## From a matrix H22, it builds a list with <tt>cardinal</tt> elements with each element containing the values {+-( H22[:,beg_index : end_index] * v + s), v, TRUNCATED(+-( H22[:,beg_index : end_index] * v + s))} with v a full-weight vector 
# @param slices_merges the division of the partial merges 
# @param even if False we negate the linear combination
# @returns a sorted list according to <tt>linear_combin_trunc</tt>
def build_list_leaves(H22, slices_merges, s,cardinal, beg_index, end_index, even,times):
    beg_build = time.time()
    H_sub = H22[:,beg_index : end_index]
    Li = []
    s_target = slices_merges[0]
    cols = H_sub.ncols()
    if(cardinal > 2**(cols)):
        raise Exception("Cannot return the demanded number of linear combinations : M has " + str(H_sub.ncols()) + " columns and you ask " + str(cardinal) + " full weight linear combinations of it : impossible" )
    vector_iteration = matrix(F3, [1 for _ in range(0,cols)], ncols = 1)
    
    elmt = element_list()
    if even:
        elmt.linear_combin = H_sub*vector_iteration + s
    else:
        elmt.linear_combin = -(H_sub*vector_iteration + s)
        
    #elmt.linear_combin_trunc = elmt.linear_combin[0:s_target].list()
    elmt.linear_combin_trunc = 0
    four_pow = 1
    for i in range(0,s_target):
        elmt.linear_combin_trunc = elmt.linear_combin_trunc + int(elmt.linear_combin[i,0])*four_pow
        four_pow = 4*four_pow
    elmt.vector_value = [copy(vector_iteration)]
    #elmt.vector_value = copy(vector_iteration)
    elmt.begin_index = beg_index
    elmt.end_index = end_index
    Li.append(elmt)

    for _ in range(1,cardinal):
        # We update vector_iteration by adding ones (and taking care that we never have zeros)
        # until we have sucessfully incremented it
        # We cannot update the vector at the end of the loop and need to start the loop at one, because of case cardinal = 2**(cols)
        index = 0
        condition = True
        while(condition):
            if(vector_iteration[index,0] == F3(2)):
                vector_iteration[index,0] = 1
                index = index+1
                if(index == cols):
                    raise Exception("Problem when updating vector_iteration: out of range index")
    
            else:
                condition = False
                vector_iteration[index,0] = F3(2)
                
        elmt = element_list()
        if even:
            elmt.linear_combin = H_sub*vector_iteration + s
        else:
            elmt.linear_combin = -(H_sub*vector_iteration + s)
        #elmt.linear_combin_trunc = elmt.linear_combin[0:s_target].list()
        elmt.linear_combin_trunc = 0
        four_pow = 1
        for i in range(0,s_target):
            elmt.linear_combin_trunc = elmt.linear_combin_trunc + int(elmt.linear_combin[i,0])*four_pow
            four_pow = 4*four_pow
        elmt.vector_value = [copy(vector_iteration)]
        #elmt.vector_value = copy(vector_iteration)
        elmt.begin_index = beg_index
        elmt.end_index = end_index
        Li.append(elmt)

    beg = time.time()
    Li.sort()
    end = time.time()
    times["Sort"] = times["Sort"] + (end-beg)
    end_build = time.time()
    times["Build list leaves"]  = times["Build list leaves"]  + end_build - beg_build


    return Li


