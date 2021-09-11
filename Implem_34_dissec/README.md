# Usage
main.sage [-r R] -n N [-v V]; where 
- R denotes the desired number of runs of the full algorithm with new random parity-check matrix and  syndrome (useful for benchmarks)
-  N the length of the code (note that there are currently three possible n implemented)
-  V the verbosity

The parity-check matrix and the syndrome are sampled uniformly at random, but this can be easily changed in <tt>main.sage</tt>.
If you have problems of kind "not found packages" rename <tt>merges.spyx</tt> to <tt>merges.sage</tt> (and do the same renaming in the file <tt>main.sage</tt> when loading the files). 

# About the implemented values of n, k and l
There are three implemented values, the smallest one solves the SDP in a few seconds while the largest one yields to running times up to a few hours.
The reason why the code size can not be chosen arbitrarily is that in such case, the asymptotic quantities are not enough to guide us, and there are a few choices to make.
For example, if \f$n=560; k= 379 ,l = 34\f$, then \f$(k+l)/64 = 6.453125 \f$
(remember that we want to partition k+l columns of the matrix called "H22" into 64 sets). 
Therefore, to some lists will be allocated linear combinations of 6 columns of H22, and for others 7. Moreover, the (theoretical) required list cardinal to be in a constant amortized time can not be reached with lists which have only 6 columns of H22.
Additionally, when 8 does not divide l, we also need to adjust by hand how do we want to spread the intermediate targets (this is done through the array named <tt>slices_merges</tt>)

To add to this implementation a new value of n, one needs therefore to make the following changes:
- In file <tt>main.sage</tt>, add a case args.codelength==n
- In <tt>layered_dissec.sage</tt>, add one case for the value of \f$n\f$ and assign the corresponding value of <tt>slices_merges</tt> (in the i^eth merge, the merge will be done with a target size of <tt>slices_merges[i-1]</tt>). 
- In <tt>build_list.sage</tt>, add one case for the value of \f$n\f$ in the function <tt>build_lists_division</tt>. The key line is <tt>L_i = build_list_leaves(H22,., .,cardinal_list,beg_index,end_index,.,.)</tt>, which uses <tt>H22[beg_index:end_index]</tt> to build the list <tt>L_i</tt> of size <tt>cardinal_list</tt>


