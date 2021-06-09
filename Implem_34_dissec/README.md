# Implementation of 4-dissection with 3 levels

# Aim
This project is a PoC which implements a 4-dissection with 3 levels embedded in an Information-Set Decoding framework.
It aims at solving a ternary syndrome decoding with large weight.

# Requirements
* Sagemath 9.2
* Python 3.9.2

# Usage

File <tt> main.sage </tt> provides a minimal example, with \f$n=560\f$ and Wave asymptotic parameters. <br />
Be careful, if you want to use one other \f$n\f$ value than the one provided in this example, you need to change leaves list initialisation in 
function  <tt>four_tree_list_tree</tt> : this is a quite tedious process.

