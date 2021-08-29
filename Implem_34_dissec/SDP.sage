## @package SDP contains the function doing the PGE and the one solving the SDP

## From a given parity check matrix H with associated syndrome s, 
# It transforms H and s to have a partial reduced form
# @returns Hp, sp, P and S such that
# * \f$Hp = S  H  P \f$ and has the needed form in PGE+SS framework
# * \f$sp = S  \cdot s \f$
# * \f$ Hp \cdot e = sp \iff H \cdot P \cdot e = s\f$
# <br> The algorithm is randomized
def build_H(H:matrix, s:matrix,n:int,k:int,l:int):
    repeat = True
    H_original = copy(H)
    while (repeat):
        H = copy(H_original)

        # Seek for a good permutation pi such that the sub matrix H11 of  pi * His invertible  
        P  = Permutations(n).random_element()
        H.permute_columns(P)
        H11 = H[0:n-k-l, 0:n-k-l]
        while(not(H11.is_invertible())):
            H = copy(H_original)
            P = Permutations(n).random_element()
            H.permute_columns(P)
            H11 = H[0:n-k-l, 0:n-k-l]

        # From now, we will work on  H * Id.permute_columns(P)
        # Thus if we have e tel such that H * pi * e = s, then the solutions of H * v = s is pi * e


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



        top_part = (H11_invert).augment(matrix(F3, nrows = n-k-l, ncols = l))
        bottom_part = A.augment(B)
        S = top_part.stack(bottom_part)

        H = S*H
        


        # Now
        # H_new * e = s_new
        # <=> S H pi e = s_new
        # <=> H pi e = S-1 s_new
        # so that we want on veut S-1 s_new = s
        # So  s_new = S s
        s = S * s

        repeat = not(rank(H22) == l) # need full rank matrix

    return(H, s, P, S)

## Main implementation of SDP
# @param H parity check matrix
# @param s syndrome target
# @param n,k,l usual code parameters
# @param weight weight target
# @param verbose set it to True to have the intermediate lists cardinal printed
# @returns solution of SDP and numbe rof needed iterations
def SDP(H:matrix,s:matrix,n:int,k:int,l:int,weight:int, verbose=False):
    nb_iter = 0
    perm = Permutations(n).identity()
    times = {}

    target_weight =  weight-k-l #for the probabilistic step
    print("")

    # Transform H such that it has desired form
    (Hp,sp,perm,S) = build_H(H, s,n,k,l)

    sp = vector(sp)
    H22 = Hp[n-k-l:n-k,n-k-l:n]
    H12 = Hp[0:n-k-l,n-k-l:n]

    # The loop stops whenever a solution is found. In such cases it returns the solution of the SDP and the  numbe rof needed iterations
    while(True): #check about rank of instanced
        start = time.time()
        nb_iter = nb_iter + 1


        times["Sort"] = 0
        times["Build list leaves"] = 0
        times["Four dissection (partionned)"] = 0
        times["- Add intermediate target in L1 and L3 "] = 0
        times["- Update modulo"] = 0
        times["- Merging on-the-fly"] = 0
        times["---- Check on-the-fly"] = 0
        times["- Merging not on-the-fly"] = 0
        times["--Appending elements "] = 0
        times["--Creating new elements (partitionned) "] = 0
        times["---- New element constructor + update linear_combin"] = 0
        times["---- Update truncated 'vector' "] = 0
        times["---- Update candidate vector of SDP "] = 0


        
        # Do the 3,4-dissection
        e = layered_dissec(n,k,l,  H22, H12, sp, target_weight, times,verbose )
        end = time.time()
        print(str(end-start))
        print("Global : time spent sorting " + str(times["Sort"]))
        print("Time spent building the leaves lists " + str(times["Build list leaves"] ))
        print("Time spent in four dissec (partitioned) " + str(times["Four dissection (partionned)"] ))
        print("- Add intermediate target for L1 and L3 "+ str(times["- Add intermediate target in L1 and L3 "] ))
        print("- Update modulo "+ str(times["- Update modulo"]))
        print("- Time spent merging on-the-fly " + str(times["- Merging on-the-fly"]))
        print("---- check on-the-fly " + str(times["---- Check on-the-fly"]))
        print("- Time spent merging " + str(times["- Merging not on-the-fly"] ))
        print("--appending elements " + str(times["--Appending elements "] ))
        print("--creating new elements (partitionned) " + str(times["--Creating new elements (partitionned) "]))
        print("---- new element constructor + matrix operations (add/substract) " + str(times["---- New element constructor + update linear_combin"]))
        print("---- Fill-in the truncated 'vector' " + str(times["---- Update truncated 'vector' "] ))
        print("---- Update candidate vector of SDP  " + str(times["---- Update candidate vector of SDP "]))

        if e != 0:
            PI = identity_matrix(F3,n)
            PI.permute_columns(perm)
            return PI*e,nb_iter

