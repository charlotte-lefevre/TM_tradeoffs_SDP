## @package main Quick starter to use 3,4-dissection

import numpy as np

load('four_dissec.sage')
load('utils.sage')
load('syndrome_decoding.sage')

F3 = Zmod(3)

## Code length
n = 560 
## Code dimension
k = 379
## sub-target size
l = 34
## Target weight size
W =0.948366

## Main function which shall be used to test POC of ternary syndrome decoding using a 3,4-dissection
# It samples uniformly a partity check matrix and syndrome and looks for solution to SDP
def main():
    H = random_matrix(F3, nrows = n-k, ncols = n)
    s = random_vector(F3, n-k)
    w = int(n*W)


    print("Run main algo ..")
    e = SDP(H,s,n,k,l,w)
    e = vector(e)
    print("Solution returned ! Let us check is it works ..")

    if (H*e - s == 0):
        print("It works :) ")
    else:
        print("\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!It does not work .. ! ")
        H = random_matrix(F3, nrows = n-k, ncols = n)
        s = random_vector(F3, n-k)
        e = SDP(H,s,n,k,l,w)
        e = vector(e)
        print((H*e - s == 0))

main()