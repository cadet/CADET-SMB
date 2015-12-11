FOURTH VERSION

In this version, the regular Simulated Moving Bed (SMB) code was developed, which both integrated the four_column case and eight_column case. So it is a general one now. Two cases was replicated, while the case I is from the paper, Numerical Method for Accelaration Calculation of Cyclic Steady State of ModiCon-SMB-Processes. It contains both the four_column case and eight_column case; 

k_a & k_d    : [5.72 7.7];
d_{col}      : 0.02 m;
L            : 0.25 m;
Porosity     : 0.83;
d_{par}      : not given, in my case, 0.0005 m;
Switching time          : 180 s;
concentration of feed   : 0.55 g/l;
Recycle flow rate       : 9.62e-7 m^3/s;
Feed flow rate          : 0.98e-7 m^3/s;
Desorbent flow rate     : 1.96e-7 m^3/s;
Extract flow rate       : 1.54e-7 m^3/s;
Raffinate flow rate     : 1.40e-7 m^3/s;

while the case II is from the paper, Model-based Control of a Simulated Moving Bed Chromatographic Process for the Separation of Fructose and Glucose. It is a merely eight_column case.

k_a & k_d    : [0.54 0.28];
d_{col}      : 2.6 cm;
L            : 53.6 cm;
Porosity     : 0.38;
d_{par}      : 0.325 cm;
Switching time          : 1552 s;
concentration of feed   : 0.5 g/cm^3;
Recycle flow rate       : 0.1395 cm^3/s;
Feed flow rate          : 0.0200 cm^3/s;
Desorbent flow rate     : 0.0414 cm^3/s;
Extract flow rate       : 0.0348 cm^3/s;
Raffinate flow rate     : 0.0266 cm^3/s;


Both cases, actually three cases involved, were reproducible. In regard to the simulation time, 

	- four_column case
		Time: 159.37 sec
		Round: 23
		Switch: 92

	- eight_column case
		Time: 997.72 sec
		Round: 51
		Switch: 408


-----------------------------------------------------------------------------
THIRD VERSION 

All the bugs mentioned above were all fixed in the third version. The reproduced data were from the paper Numerical Method for Accelaration Calculation of Cyclic Steady State of ModiCon-SMB-Processes.

k_a & k_d    : [5.72 7.7];
d_{col}      : 0.02 m;
L            : 0.25 m;
Porosity     : 0.83;
d_{par}      : not given, in my case, 0.0005 m;
Switching time          : 180 s;
concentration of feed   : 0.55 g/l;
Recycle flow rate       : 9.62e-7;
Feed flow rate          : 0.98e-7;
Desorbent flow rate     : 1.96e-7;
Extract flow rate       : 1.54e-7;
Raffinate flow rate     : 1.40e-7;

The convergence to the Cyclic Steady State was achieved after 68 switching. This resulted in approx. 99 second on my desktop under the relative tolerence 1e-4. 

-----------------------------------------------------------------------------
SECOND VERSION
Since there is not so many cases of four-column SMB (at lease I did not find
so many), I did not try so many cases. 

The parameters for the numerical test for the second version is withdrawn from
the paper Numerical Method for Accelerated Calculation of Cyclic Steady State
of ModiCon-SMB-Processes. It seems that the paper, On Simplified Modeling Approaches
to SMB Process, use the same data for test. However, in both cases, the
parameter PARTICLE_RADIUS can not be found.

As for the bug, T_{column} less than the switch time T_{s}. It was alleviated
by a trick (assuming the tubes are employed between the adjacent columns).

And the another trick is that the injection from the FEED inlet were adopted
periodicly (four switches one injection). In regard to the reason, just for
implementation. Lastly, the so-called CSS was reached after approximately one
hour simulation. 

-----------------------------------------------------------------------------

This is the code of the four-column Simulated Moving Bed (SMB) chromatography for separating two components mixture. And it has also extended to the situation of eight-column. Although it is self-explanatory, here still comes a simple introduction.

In the CADET, the general rate model (GRM) is adopted, while in this case, GRM is simplified into the equilibrium dispersive model. In other words, the beads are treated as solid and no porous diffusion, and film diffusion are taken into consideration (please see the parameter configuration in getParameter.m). Additionally, the simplest Linear adsorption isotherm is implemented in this case. 

In each connectiong nodes between two columns, the mass conservation is used to calculated the inlet concentration of next column (see massConservation.m).

secColumn.m is the subroutine used for calculation of single column.

simulatedMovingBed is the main funciton. 

This code is based on CADET, and the code for simulated moving bed is put in the file "simulatedMovingBed". The binary file is for ubuntu 14.04. If you are using other OS, you need to build the binary file from the source file (http://github.com/modsim/CADET.git) yourself; if you are using Windows, you can download the file directly from https://github.com/modsim/cadet/releases. 

