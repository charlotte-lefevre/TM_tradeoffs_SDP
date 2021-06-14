## @package main
# Main file which shall be used to test POC of ternary syndrome decoding using a 3,4-dissection
import numpy as np

load('four_dissec.sage')
load('utils.sage')
load('syndrome_decoding.sage')

F3 = Zmod(3)

# Used parameters
n = 560 ; k = 379 ; l = 34
W =0.948366

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


    PI = identity_matrix(F3,n)
    PI.permute_columns(perm)
    if ((H != S* H_original* PI)):
        print("problem !!!!!!!!")

        print(H11.is_invertible())
    print("")
    return(H, s, perm, S)



# Used parameters
n = 560 ; k = 379 ; l = 34
W =0.948366

for _ in range(0,50):
    H = random_matrix(F3, nrows = n-k, ncols = n)
    s = random_vector(F3, n-k)

    (Hp,sp,perm,S) = build_H(H, s,n,k,l)