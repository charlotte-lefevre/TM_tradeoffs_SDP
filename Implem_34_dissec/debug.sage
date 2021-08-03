## @package debug 

## Tests whether the elements_lists in L behave well : the vector_value multiplied by the matrix H must give the linear combination (of the opposite of the linear combination if even is False)
# @param s the syndrome if it needs to be included in the verification (so only if include_s = True)
# @param missing_slices when the linear combinations only keeps track of the coordinates starting from missing_slices
"""def oracle_test_integrity(L, H, text, even, s=0, missing_slices = 0, include_s = False):
    for element in L:
        lc = element.linear_combin
        l_v = element.vector_value
        H_sub = H[:,element.begin_index:element.end_index]
        list_v = []
        for sub_v in l_v:
            list_v = list_v + sub_v.list()
        vect = matrix(F3, list_v, ncols = 1)
        #print(vect)
        #print(H_sub)

        completed_vector = matrix(F3,[0 for _ in range(0,missing_slices)] + lc[missing_slices:len(lc)].list(), ncols = 1)
        if include_s:
            if even:
                if H_sub*vect -s != completed_vector :
                    print("Not working " + str(text) +" the vector * matrix gives us")
                    print(H_sub*vect -s)
                    print("versus obtained vector")
                    print(completed_vector)
            else:
                if H_sub*vect -s != -completed_vector :
                    print("Not working " + str(text)+" the vector * matrix gives us")
                    print(H_sub*vect -s)
                    print("versus obtained vector")
                    print(-completed_vector)                
        else:
            if even:
                if H_sub*vect != completed_vector:
                    print("Not working " + str(text)+" the vector * matrix gives us")
                    print(H_sub*vect)
                    print("versus obtained vector")
                    print(completed_vector)
            else:
                if H_sub*vect != -completed_vector:
                    print("Not working " + str(text)+" the vector * matrix gives us")
                    print(H_sub*vect)
                    print("versus obtained vector")
                    print(-completed_vector)"""

## Tests whether the elements_lists in L behave well : the vector_value multiplied by the matrix H must give the linear combination (of the opposite of the linear combination if even is False)
# @param s the syndrome if it needs to be included in the verification (so only if include_s = True)
# @param missing_slices when the linear combinations only keeps track of the coordinates starting from missing_slices
def oracle_test_integrity(L, H, text, even, s=0, include_s = False):
    for element in L:
        if even:
            lc = element.linear_combin
        else:
            lc = -element.linear_combin
        if include_s==False:
            syndrome = matrix(GF(3), ncols = 1, nrows = H.nrows())
        else:
            syndrome = s
        l_v = element.vector_value
        candidate = matrix(GF(3), nrows=0,ncols=1)
        for vi in l_v:
            candidate = candidate.stack(vi)

        H_sub = H[:,element.begin_index:element.end_index]

        #print(vect)
        #print(H_sub)
        if H_sub*candidate -syndrome != lc :
                print("Not working " + str(text) +" the vector * matrix gives us")
                print(vector(H_sub*candidate -syndrome))
                print("versus obtained vector")
                print(vector(lc))
                print("substract")
                print(vector(H_sub*candidate -syndrome -lc))
                exit()


## copied a list 
def copy_list(L_tocopy):
    L_return = []
    for element in L_tocopy:
        el = element_list()
        el.linear_combin = copy(element.linear_combin)
        el.linear_combin_trunc = copy(element.linear_combin_trunc)
        el.vector_value: copy(element.vector_value)
        L_return = L_return + [el]
    print("et")
    print(len(L_return))
    return L_return