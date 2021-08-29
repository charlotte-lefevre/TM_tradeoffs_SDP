## @package element_list 

## @brief Structure stored in a list 
# An element_list is defined by
# * <tt>vector_value</tt>: the candidate for the sub-SDP, call it \f$e_2\f$. It is implemented as a list of matrices because it is better for performances
# * <tt>begin_index</tt>, <tt>end_index</tt>: the indexes of the parity-check matrix's columns which has been used to compute the value of <tt>linear_combin</tt>. In other words \f$H_2\f$ := H[,begin_index:end_index]
# * <tt>linear_combin </tt>: the value of \f$H_2e_2\f$ (resp. \f$-H_2e_2\f$ ) if the element belongs to a list with even (resp. odd) index. 
#  We sometimes negate the linear combination to make mergings easy. Indeed, we want collisions of form \f$H_{i}e_{i} + H_{i+1}e_{i+1} = 0\f$, so \f$H_{i}e_{i} = -H_{i+1}e_{i+1} \f$.
# * <tt>linear_combin_trunc</tt>: it is equal to <tt>linear_combin</tt> on a slice which is the one on which we want to do a merging. For performance reasons it is implemented as an integer in base 4. For example, to the vector (2,0,1) is assigned the integer 1*4**0 + 0*4**1 + 2*4**2
class element_list:
    def __init__(self) -> None:
        self.vector_value: list[matrix] = None # because we concatenate vectors together we impelment this as a list of vectors
        # the vector is a slice of H's columns between begin index and end index 
        self.linear_combin:matrix = None
        self.linear_combin_trunc:int = None
        self.begin_index = 0
        self.end_index = 0

    ## The ">"operator is defined here in order to be able to use built-in sorting methods
    # As we sort according to  linear_combin_trunc, we overload the operator according to these fields    
    def __gt__(self, other):
        return (self.linear_combin_trunc>other.linear_combin_trunc)




                