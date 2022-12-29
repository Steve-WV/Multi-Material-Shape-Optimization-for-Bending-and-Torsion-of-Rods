Guidelines how to use the code: 

REQUIREMENTS:
 - Cmake version no less than 3.7.2
 - A suitable C++ compiler
 - The linear algebra library Eigen
 - The Boost libraries


After compilation using cmake and make, a target called "s2d" is produced 
which can be executed by: s2d -f [path to]/Square_144K.obj 

By this the simulation is started on the reference domain Omega=[-1,1]^2.

The main code of the simulation is provided in the executable 
file "pde.cpp". Here the relevant parameters are described 
in lines 16-27, where e.g. the weigthing factors can be changed 
in order to simulate certain gradient flows. 

The default configuration of the weigthing factors is appropriate 
to the fourth experiment provided in the thesis Section 7.
