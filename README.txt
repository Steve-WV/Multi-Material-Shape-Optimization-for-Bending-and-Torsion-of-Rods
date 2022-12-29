Guidelines how to use the code: 

REQUIREMENTS:
 - Cmake version no less than 3.7.2
 - A suitable C++ compiler
 - The linear algebra library Eigen
 - The Boost libraries


After compilation using cmake and make, a target called "s2d" is produced 
which can be executed by: s2d -f [path to]/circle160K.obj -c [path to]/default1.cfg 

By this the simulation is started on a circular reference domain.

The main code of the simulation is provided in the executable 
file "pde.cpp". 

The relevant parameters are included via the config file "default1.cfg" where factors can be changed 
in order to simulate different shape optimization problems.

The default configuration of the weigthing factors is appropriate 
to the fourth experiment provided in [1, Section 5].
