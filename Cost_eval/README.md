# Computation of time Memory trade offs for ternary syndrome decoding with large weight 

# Aim
This project aims at computing time-memory trade-offs costs for ternary syndrome decoding in \f$F_3\f$ with large weight.
The most efficient state-to-the-art algorithms use a framework called PGE+SS framework. In the latter, the most expansive part is due to a step called SubsetSum step.
Here we aim at computing the Time/Memory cost of a pallet of birthday-based algorithms to solve the Subset Sum step
We aim at providing a pallet of algorithms, each one presenting a different trade-off
Everything here has been created with Wave signature scheme in mind, so only the most interesting tarde offs in wave signature scheme case are exhibited <br />
**Warning** There is No implementation of any algorithm here, the aim here is just to compute asymptotic cost of algorithms to help the user to select his prefered trade off.


# Requirements
* Sagemath 9.2
* Python 3.9.2

# Common notation and terminologies

In the functions used here, the notation is quite redundant. We define here the most common terminologies.

## Two different regimes

The cost can be computed in two regimes, the *asymptotic* one and the so-called '*real* 'one.
In the *asymptotc* regime, the code size \f$n\f$ is not fixed, and the aim is to derive an estimation of the cost when \f$n \rightarrow +\infty\f$.
In other words, multplicative constants are not taken into account, and only eexponential coefficients are considered.
On the countrary, the '*real*' regimes aims at provideing more accurate cost estimation. Here, the code size is fixed and multilicartive constants are consdiered. <br />

In both cases, the memory and time costs are given in the *logarithmic* scale, and relative to \f$n\f$ in asymptotic cost (*viz.*\f$m = \log_2(M)/n\f$)

## Most common parameters

To avoid a heavy redundancy, we define here the most used parameters in the project's functions

### Code parameters
*  \f$n\f$ is the size of the code. 
*  \f$k\f$ is the code dimension. In asymptotic notation, \f$k := R \cdot n\f$
*  \f$l\f$ is size of the target (or size of sub syndrome decoding problem). In asymptotic notation, \f$l := Re \cdot n\f$
*  \f$w\f$ is the weight target in SDP. In asymptotic notation, \f$w := W \cdot n\f$
* \f$Skkl\f$ is the \f$\log_3\f$ os desired number of solutions for Probabilistic step (divided by \f$n\f$ in the  *asymptotic* regime)



### (Layered) dissection parameters
*  \f$r\f$ is such that base dissection is an \f$r\f$-dissection
* \f$h\f$ is the number of levels in the tree
* \f$c\f$ is the dissection work-factor. In other words, an exhaustive dissection with memory \f$M\f$ has time cost \f$M^c\f$


# Sample examples
Examples are contained in "Examples" directory 
*  <tt>Mainplots.ipynb</tt> contains our main resuts
*  <tt>DissectionPlot.ipynb</tt> calles the dissection algorithm
* <tt>plot_one_layered_dissection.ipynb</tt> plots the time memory graph of one dissection



