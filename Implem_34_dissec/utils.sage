## @package utils util functions


## from M, creates cardinal linear combinations of matrix M
# Returns the vectors under a matrix form for practicla reasons
# @return L1 the linear combinations 
# @return L2 the associated vectors (in a mtrix form)
# if lc = L1[i], e = L2[i] then M*e = lc
def return_all_LC_matrix(M,cardinal):
    cols = M.ncols()
    if(cardinal > 2**(cols)):
        raise Exception("Cannot return the demanded number of linear combinations : M has " + str(M.ncols()) + " columns and you ask " + str(cardinal) + " full weight linear combinations of it : impossible" )
    
    vector_iteration = matrix(F3, [1 for _ in range(0,cols)],ncols = 1)
    
    vectors = [] ; linear_combin = []
    
    linear_combin.append(M*vector_iteration)
    vectors.append(copy(vector_iteration))

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
                    raise Exception("Badly managed errors in return_all_LC_matrix_no_recursivity : this case was supposed to be taken care before")
    
            else:
                condition = False
                vector_iteration[index,0] = F3(2)
                
        linear_combin.append(M*vector_iteration)
        vectors.append(copy(vector_iteration))

    return(linear_combin, vectors)


## From a submatrix H_sub, it builds a list containing cardinal_list linear combinations of vectors
# if index is equal to 0, it takes into account the syndrom (that is populates the list with "linear_combinations - s")
# The result is returned as a dict "Linear combination ---> vector"
# @attention for practical reasons everything is returned as a matrix object (appending is much easier)
def build_list_matrix(H_sub, index, s,cardinal_list):
    L = {}
    (Mtimesvs, vs) = return_all_LC_matrix(H_sub,cardinal_list) #on a la liste des CL des colones (target) avec les vecteurs associés
    for j in range(0,len(Mtimesvs)):
        M_times_v = 0
        if (index == 0): # We shall include the syndrome
            M_times_v = (Mtimesvs)[j] - s
        else:
            M_times_v = (Mtimesvs)[j]
        v = vs[j]
        M_times_v.set_immutable()
        L[M_times_v] = v
    return L
