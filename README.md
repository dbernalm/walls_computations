This is a algorithm to find the lambda walls for complexes in $D^b(X)$ ($X$ is a 3fold of Picard 1) with Chern character $(-R,0,D,0)$. Those are cases of: instanton sheaves, shifts of ideal sheaves of curves (twisted by a certain number), and structure sheaves of one-dimensional curves.

So far, the case $X=P^3$ and $X$ abelian have been completed successfully. A part of the program is concerned with Fano 3folds of Picard 1, but there is still some work missing.

More specifically, if $A$ is a destabilizing object of $E$ with Chern character $ch^\beta = (-R, 0, D, 0)$ with $\beta= \frac{1}{k}$ or $0$ and $ch^\beta A = (r,c,d,e)$, then the algorithm computes all possibilities of such Chern character using the integral conditions (N $\beta$ 1), (N$\beta$2) and (N$\beta$3) and the numerical conditions (I$\beta$1), (I$\beta$2) and (I$\beta$3) which can be found in [arXiv:2511.13930](https://arxiv.org/abs/2511.13930). In order to improve the speed for computing the cases $D=3$ with $\beta=\frac{1}{3}$ and $D=4$ with $\beta = \frac{1}{4}$, the algorithm relies on the multiprocessing package from Python, and thus the program would work better on a physical computer with at least 4 cores. 

In this case, we define an object by calling its class

    a=Sheaf(R,D,k)
where $\beta=\frac{1}{k}$. If $k=1$ then the program assumes $\beta=0$. After it you use the method 

    a.num_dest(d)
where the argument counts the number for which you want the process to be divided. This is used for the multiprocessing aspect, each process analyzes simultaneously $\frac{1}{d}$ of the list. The output is then a list of numbers of the form $\left([r_i]_{i=1}^l,c,d,e\right)$ where the $r_i$ are the possible ranks, and $c,d,e$ are the Chern characters of the destabilizing object generating the numerical wall.
