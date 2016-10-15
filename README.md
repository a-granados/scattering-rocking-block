# Computation of the Scattering map for a non-smooth system: the rocking block
This software was to compute the numerical results presented in the paper "The scattering map in two coupled piecewise-smooth systems, with numerical application to rocking blocks" (see <a href="http://agranados.no-ip.org"> my webpage</a>). See the paper for details on the system and parameters definition.


**General description**

These routines execute the method described in section 6 of the paper. The main goal is to find, for given coordinates at the reference manifold, (tau,v,s), the intersection of the corresponding heteroclinic fiber with the section given by x=0.<br>
This intersection is given by the parameter zeta. At the highest level, the program performs a Bolzano method to find zeta^* that gives a zero of the distance between the intersection of the stable and unstable fiber with x=0.<br>
The intersections between the stable and unstable fibers with the switching manifold x=0 are found by means of a secant method. At each step, several initial conditions are flowed in parallel, to accelerate the process, and one decides whether trajectories scape (invariant manifold lies below) or come back (invariant manifold lies above). Each client integrates one initial condition using multiple precision libraries ( arprec), and requires login via ssh or rsh. In a network providing around 100 cores, the overall process may take around 8 hours when requiring a tolerance of 10^{-27}.<br>
Once the heteroclinic fiber has been found, another set of programs evaluate the scattering map and perform a 2d mapping varying the initial conditions tau, v or s.

**Main files and routines**

Each routine contains comments that may help to follow what it is being doing.
* system_params.dat contains system parameters: omega (frequency of the forcing) te (rescaling for the amplitude of the forcing) and tk (rescaling for the amplitude of the coupling)
* clients_sm and clints_um are the list of clients to compute the intersections between x=0 and the stable and unstable fibers, respectively 
* integrate_forwards.cpp, find_sm.cpp, integrate_bacwards.cpp and find_um.cpp integrate the system forwards and backwards 
* find_sm.sh and find_um.sh the secant method to find the intersection between x=0 and the stable and unstable fibers, respectively.
* find_hetN.sh launches the whole process and contains the Bolzano's loop. Last version of the script is 4.2, and is the one that I recommend to use.
