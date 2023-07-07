# Randomized and Inverse-Free Matrix Pencil Diagonalization via Pseudospectral Divide-and-Conquer

This repository contains a Matlab implementation of a recent randomized (and inverse-free) algorithm for diagonalizing any matrix pencil $(A,B)$. The main routine **RPD** takes as inputs the two matrices $A,B \in {\mathbb C}^{n \times n}$ and a desired backward error $\varepsilon < 1$.

Requiring only $||A||_2, ||B||_2 \leq 1$, **RPD** produces invertible matrices $S,T$ and diagonal $D$, such that 

$$ \max \left[||A - SDT^{-1}||_2, ||B - SIT^{-1}||_2 \right] \leq \varepsilon $$

with high probability. The algorithm is detailed extensively in a [forthcoming paper.](https://arxiv.org/abs/2306.03700)

## Usage
In addition to the source code we include the following scripts, each of which applies **RPD** to a different example pencil.
1. planted_spectrum.m
2. jordan_block.m
3. increasingly_singular.m
4. singular_pencil.m
   
Running each script reproduces the corresponding figures in [1].

## Contact
Thanks for taking the time to look through our work! If you have questions or feedback about using this repository, reach out to Ryan Schneider at ryschnei@ucsd.edu.

## References
<a id="1">[1]</a> 
J. Demmel, I. Dumitriu, and R. Schneider. *Generalized Pseudospectral Shattering and Inverse-Free Matrix Pencil Diagonalization.* arXiv preprint (2023).
