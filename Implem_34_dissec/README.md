# Usage
main.sage [-r R] -n N [-v V]; where R denotes the desired number of runs of the full algorithm with new random parity-check matrix and  syndrome (useful for benchmarks), N the length of the code, V the verbosity
Note that there are only three possible n implemented for now.
The parity-check matrix and the syndrome are sampled uniformly at random, but this can be easily changed in <tt>main.sage</tt>.
If you have problems of not found packages, rename <tt>merges.spyx</tt> to <tt>merges.sage</tt> (don't forget to change the name when loading the corresponding file in the file <tt>main.sage</tt>

# How to change the value of n, k and l
Changing the value of \f$n\f$ is not a trivial task
Indeed, the theoretical number of columns of H allocated to each list on the leaves is given by \f$(k+l)/64\f$. However, in practise this number if not an integer, so some manual work needs to be done. For example, with \f$n=560; k= 379 ,l = 34, (k+l)/64 = 6.453125 \f$, so to some lists will be allocated 6 columns of H's submatrix H22, and for others 7. Moreover, the required list cardinal to be in a constant amortized time can not be reached with lists which have only 6 columns of H22.
Additionally, when 8 does not divide l, we also need to adjust by hand how do we want to spread the intermediate targets (this is done in the array <tt>slices_merges</tt>)
To change the value of n, we need to change:
- In file <tt>main.sage</tt>, add a case args.codelength==n
- In <tt>layered_dissec.sage</tt>, add one case for the value of \f$n\f$ and assign the corresponding value of <tt>slices_merges</tt>. 
- In <tt>build_list.sage</tt>, add one case for the value of \f$n\f$ in the function <tt>build_lists_division</tt>.
