*********************************************************************
********            GillespieSSA     (release 1)             ********
*********************************************************************

March 8, 2010
Tomas Vejchodsky
www.math.cas.cz/vejchod/gillespiessa.html

Implementation of the Gillespie Stochastic Simulation Algorithm for
Matlab. The code performs simulation of the dynamics of a system
of q chemical reactions of Nchs chemical species.
The i-th, i=1,2,...,q, has the form
  c(1,i) A(1) + c(2,i) A(2) + ... + c(Nchs,i) A(Nchs)
  --->
  p(1,i) A(1) + p(2,i) A(2) + ... + p(Nchs,i) A(Nchs)
where
  c(j,i) is the coefficient at the j-th chemical species at the
         reactant side of the i-th reaction
  p(j,i) is the coefficient at the j-th chemical species at the
         product side of the i-th reaction
  A(j) stands for the j-th chemical species, j=1,2,...,Nchs
The dynamics of this reaction is described by the corresponding rate
constant knu(i).

Quick Start
===========
- Unpack the archive  gillespiessa.tar.gz  with source files.
  Keep everything in the same subdirectory (say gillespiessa).
- Run Matlab.
- In Matlab, change current directory to say gillespiessa.
- Run the test computations.
  >> test1
  >> test2
  >> test5   [this test produces two output files]
  >> test5_gssa_anim    [example of the animation]

User's Guide
============
The Gillespie Stochastic Simulation Algorithm (SSA) follows these
steps:
 (a) Generate two randon numbers r1 and r2 uniformly distributed
     in (0,1).
 (b) Compute alpha0 = \sum_{i=1}^q alpha_i(t), where alpha_i(t)
     stands for the propensity function of i-th reaction.
 (c) The next reaction take place at time t+tau, where
     tau = 1/alpha0 * log(1/r1)   (log is the natural logarithm).
 (d) Determine which reaction occurs. Find j such that
     r2 >= 1/alpha0 \sum_{i=1}^{j-1} alpha_i(t)
     r2 <  1/alpha0 \sum_{i=1}^j     alpha_i(t)
     Update the number of reactants and products of the j-th reaction.
 (e) Go to (a) with t = t + tau


The main function is "gillespiessa.m". This function performs the
stochastic simulation and returns its complete history. It has five
input parameters:
knu ... 1-by-q row vector; rate constants of the chemical reactions;
      knu = k/nu^(order-1), where nu stands for the volume and order
      for the order of the chemical reaction
c ... Nchs-by-q matrix; c(j,i) is the coefficient at the j-th chemical
      species at the reactant side of the i-th reaction
p ... Nchs-by-q matrix; p(j,i) is the coefficient at the j-th chemical
      species at the product side of the i-th reaction
A0 ... Nchs-by-1 column vector; A0(j) is the number of molecules of
      j-th chemical species at the initial time t=0
tfin ... the final time; the simulation stops if it reaches the time
      tfin

Remark: Nchs is determined as size(A0,1) and q as size(knu,2).

There are two outputs of "gillespiessa.m":
A ... Nchs-by-N matrix; A(j,m) is the number of j-molecules in time
      interval [t(m), t(m+1)]. The initial value is A(:,1) == A0(:).
t ... 1-by-N vector; A number of time values in which the number
      of molecules changes.

The function "gillespiessa_anim.m" does the same as "gillespiessa.m",
but it shows in addition a real time animation of the progress of the
simulation.

The function "mkplotdata.m" modifies the output data from
"gillespiessa.m" such that the modified data (plotted by the standard
Matlab commands) correspond to the proper interpretation of the
results of the Gillespie SSA.

For further information see the description directly in the
particular m-files.

Enjoy computations with GilespieSSA!
