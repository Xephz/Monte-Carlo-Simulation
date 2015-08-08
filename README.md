#Monte Carlo Simulation

![The result of a simulation (a cluster crystal) displayed in Jmol](img.jpg?raw=true)

This program was written for the purpose of running Monte Carlo simulations with a generalised exponential model potential to explore soft matter systems. It can however easily be modified to simulate any potential by changing the function PairPotential in coreMCfunc.c

The xyz files outputted by this program are intended for use with Jmol to visualise the results of the simulation. The GRAPH and DIST files simply give lists of numbers representing radial distribution functions and cluster distributions respectively and be imported into your favourite data analysis tools.

The program is used by editing the GEM.c source file with the desired parameters, inputs and outputs.

A program written as part of our final year project by myself(stu.wood.94@googlemail.com) and my partner Jack Porter.
 
Note that due to licensing issues with redistribution I have subsituted the pseudo-random number generator used origionally with the funtions from stdlib. This results in much slower simulations that are relying on what are likely lower quality distributions.
