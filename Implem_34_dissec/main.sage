

## @package main Quick starter to use 3,4-dissection
F3 = GF(3)

load("debug.sage")
load( "element_list.py", "merges.spyx", "SDP.sage", "build_list.sage", "layered_dissec.sage")

import numpy as np
import time 
import argparse
parser = argparse.ArgumentParser()

# available instanciations 
#n =700; k =473; l =40
#n = 560; k = 379; l = 34
#n=875; k=591; l=48

# the seed is fixed for benchmarks
#set_random_seed(20)
#set_random_seed(100)
## Main function which shall be used to test POC of ternary syndrome decoding using a 3,4-dissection
# It samples uniformly a parity check matrix and syndrome and looks for solution to the SDP
def main(verbose= True):
    
    H = random_matrix(F3, nrows = n-k, ncols = n)
    s = random_matrix(F3, nrows= n-k, ncols = 1)

    #print("Run main algo ..")
    e, nb_iter = SDP(H,s,n,k,l,w, verbose)
    e = vector(e)
    print("Solution returned ! Let us check is it works ..")

    if (H*e - vector(s) == 0 and e.hamming_weight() == w):
        print("It works :) ")
        return e,nb_iter
    else:#means that there is a bug somewhere
        print("\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!It does not work .. ! ")
        print(e.hamming_weight())
        print(vector(H*e - vector(s)))
        return
    
        


        
parser.add_argument("-it", "--nbiter", dest="nbiter", help="Number of iterations of SDP, default to one", type=int, default=1)
parser.add_argument("-n", "--codelength", dest="codelength", help="Code length (available:875,700,560)", required=True, type=int)
parser.add_argument("-v", "--verbose", dest="verbose", help="verbose output", type=bool, default=False)
args = parser.parse_args()



if args.codelength == 700:
    n =700; k =473; l =40
elif args.codelength == 560:
    n = 560; k = 379; l = 34
elif args.codelength == 875:
    n=875; k=591; l=48
else:
    print("The code length " +str(args.codelength) + " is not implemented. Available code lengths : 700,560,875. I exit")
    exit()
w = int(0.948366*n) 
# available instanciations 
#n =700; k =473; l =40
#n = 560; k = 379; l = 34
#n=875; k=591; l=48
avg_iter = 0
for _ in range(0,args.nbiter):
    # set verbose to true to have more information about length of lists
    _, nb_iter =  main(verbose = args.verbose)
    avg_iter = avg_iter + nb_iter

print("Average number of iterations: " + str(RR(avg_iter/args.nbiter)))
