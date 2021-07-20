# Computation of time Memory trade offs for ternary syndrome decoding with large weight 

# Aim
This project aims at computing time-memory trade-offs costs for ternary syndrome decoding in \f$F_3\f$ with large weight and Wave signature scheme parameters.
The most efficient state-to-the-art algorithms use a framework called PGE+SS framework. In the latter, the most expensive part is due to a step called Subset Sum step.
Here we aim at computing the Time/Memory cost of a range of birthday-based algorithms to solve the Subset Sum step.
Everything here has been created with Wave signature scheme in mind, so only the most interesting trade offs in the Wave signature scheme case are exhibited. <br />
**Warning** There is no implementation of any algorithm here, the aim here is just to compute the asymptotic cost of algorithms to help to select the prefered trade off.


# Common notation and terminologies

In the functions used here, the notation is quite redundant. We define here the most common terminologies.

## Two different regimes

The cost can be computed in two regimes, the *asymptotic* one and the so-called '*bit cost* 'one.
In the *asymptotc* regime, the code size \f$n\f$ is not fixed, and the aim is to derive an estimation of the cost when \f$n \rightarrow +\infty\f$.
In other words, multplicative constants are not taken into account, and only exponential coefficients are considered.
On the countrary, the '*bit cost*' regimes aims at providing a more accurate cost estimation. 
Here, the code size is fixed and multiplicative constants are considered. It aims at providing a cost in terms of number of bits stored / number of bit operations done.<br />

In both cases, the memory and time costs are given in the *logarithmic* scale, and relative to \f$n\f$ in asymptotic cost (*viz.*\f$memory = \log_2(M)/n, time = \log_2(T)/n\f$)

## Most common parameters

To avoid redundant comments, we define here the most used parameters in the project's functions:

### Code parameters
*  \f$n\f$ is the size of the code. 
*  \f$k\f$ is the code dimension. In asymptotic notation, \f$k := R \times n\f$
*  \f$l\f$ is size of the target (or size of sub syndrome decoding problem). In asymptotic notation, \f$l := Re \times n\f$
*  \f$w\f$ is the weight target in SDP. In asymptotic notation, \f$w := W \times n\f$
* \f$Skkl\f$ is the \f$\log_3\f$ os desired number of solutions for Probabilistic step (divided by \f$n\f$ in the  *asymptotic* regime)



### (Layered) dissection parameters
*  \f$r\f$ is such that base dissection is an \f$r\f$-dissection
* \f$h\f$ is the number of levels in the tree
* \f$c\f$ is the dissection work-factor. In other words, an exhaustive dissection with memory \f$M\f$ has time cost \f$M^c\f$


# Sample examples
Examples are contained in "Examples" directory 
*  <tt>Mainplots.ipynb</tt> contains our main results
*  <tt>DissectionPlot.ipynb</tt> plot the dissection algorithm
* <tt>plot_one_layered_dissection.ipynb</tt> plots the time memory graph of one dissection



