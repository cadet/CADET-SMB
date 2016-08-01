![](https://github.com/modsim/CADET/blob/master/doc/logo/CADET-GitHub.png)

# CADET-SMB

CADET-SMB is a comprehensive simulator for analysis and design of simulated moving bed (SMB) chromatographic processes. It is developed at the Institute of Bio- and Geosciences 1 (IBG-1) of Forschungszentrum Jülich (FZJ) under supervision of Dr. Eric von Lieres. CADET-SMB uses the simulation engine of the CADET framework, which provides a fast and accurate solver for the general rate model (GRM) of packed bed liquid chromatography. 

# One-column analog approach

This variant, one-column analog approach, is an attractive alternative to the standard_SMB approach on the master branch of the repository, in terms of less time-consuming on the computational side and cost-saving on the operational side. A manifest advantage of this approach is that, if the cyclic steady state (CSS) is the state in SMB processes that you are concerned, it speeds up the convergence to the CSS comparing to the standard_SMB approach.  

In one-column analog approach, the conventionally periodical SMB chromatographic processes are reproduced by a single column setup with a recycle lag of (N-1)t_s time units, where N is the column numbers and t_s is the switching time. It has been proved that this ideal one-column analog is theoretically indistinguishable from the equivalent SMB, except for the discontinuous use of the inlets/outlets (we refer you to the papers listed for details). By the way, we are not going to discuss how those processes are able to be implemented in practice, since some external tanks/tubes for the discontinuous inlets/outlets are required. 
		
		[1]  Abunasser N, Wanket P. C, Kim Y. S, et al. One-column chromatography with recycle analogous to a four-zone simulated moving bed[J], Industrial & engineering chemistry research, 2003, 42(21): 5268-5279.
		[2] Mota J. P. B, Araujo J. M. M, Single-column simulated-moving-bed process with recycle lag[J], AIChE journal, 2005, 51(6): 1641-1653.

![](https://github.com/modsim/CADET-SMB/blob/One-column_analog/doc/scheme.JPG)

*The schematic of the one-column analog approach*

As observed from above figures, there is only one column (the red marked one) packed in each SMB unit, and the other columns are replaced with tanks (dot-line rectangles) to provide necessary recycle streams. By the way, the chosen column can be an arbitrary one in the SMB unit. 

First of all, there is no adjacent columns any more in this variant as there is only one column packed. So what is the inlet concentration profile of each column? And secondly, let's have an insight into the recycle pattern of the one-column analog approach. The simulation sequence in the standard_SMB approach is clockwise, however in this approach it is counter-clockwise. (Assuming that each column is sticked with a Arabic number permanently). The simulation sequence here is the column 3 in zone III, column 3 in zone II, column 3 in zone I, and column 3 in zone IV. 


![](https://github.com/modsim/CADET-SMB/blob/One-column_analog/doc/recycle_pattern.JPG)

The recycle pattern could be briefly illustrated by the following two figures.

![](https://github.com/modsim/CADET-SMB/blob/One-column_analog/doc/analog_1.JPG)
![](https://github.com/modsim/CADET-SMB/blob/One-column_analog/doc/analog_2.JPG)



# Features

* Binary separation is available using four-zone scheme;

* Ternary separation can be achieved by using integrated five-zone scheme or eight-zone scheme;

![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_binary.JPG)
![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_ternary_5.JPG)
![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_ternary_8.JPG)

* In both binary and ternary separations, arbitrary column configurations are available, in addition to basic column configurations such as 1-1-1-1, 2-2-2-2-2, 3-3-3-3, 4-4-4-4-4;

* It converges to the same cyclic steady state (CSS) with the naive standard_SMB approach; but the trajectories are distinct and it speeds up the convergence to the cyclic steady state.

![](https://github.com/modsim/CADET-SMB/blob/One-column_analog/doc/time_comparison.JPG)

* Continuous stirred tank reactor (CSTR) and dispersive plug flow reactor (DPFR) models can be placed before and after each column to account for dead volumes in pumps, tubing, and valves;

![](https://github.com/modsim/CADET-SMB/blob/master/doc/dead_volumes.JPG)

* MATLAB interface allows to monitor the dynamic characteristics of each column in the SMB unit;

* Optimization of decision variables for improving, say, productivity, purity, and operating costs;

* Column models include transport dispersive model, equilibrium dispersive model, and general rate model;

* Wide range of standard equilibrium/isotherm models allow to simulate either pure component or multi-component/competitive behaviour;

* Further features of the CADET framework can be found at https://github.com/modsim/CADET.


# Dependency and Platforms

* Matlab (R2010b or higher);
* CADET (version 2.3.2 or later);
* platforms, please see the Dependencies section in the CADET wiki.


# Tutorial and Instructions

First, download the CADET software from https://github.com/modsim/CADET-SMB/releases, as CADET-SMB is based on the CADET simulator.
Then, download the latest release of CADET-SMB from https://github.com/modsim/CADET-SMB/releases.

## Installation of CADET

* download the latest release for your platform;
* unzip the archive to your destination directory;
* start MATLAB;
* change to the unzipped CADET directory and run installCADET.m. (You can save the MATLAB path to avoid calling installCADET.m every time you restart MATLAB);
* Try one of the examples (e.g., examples/forward/loadWashElution.m) to check if everything works.

## Installation of CADET-SMB

* create a directory, simulatedMovingBed, in your unzipped CADET directory;
* unzip the CADET-SMB archive to the simulatedMovingBed directory;
* Change the working directory to the simulatedMovingBed directory and run isSMBupdateAvailable.m script (Along side checking the existence of the newest version, it also attach the current path to the MATLAB path); 
* To test a forward simulation, copy any getParameter routine from the examples/Forward folder to the simulatedMovingBed folder and change the name to getParameters.m; then run simulatedMovingBed.m.
* To test an optimization, copy any getParameter routine from the examples/Optimization folder to the simulatedMovingBed folder and also change the name to getParameters.m; then run SMBOptimization.m.

# Demonstration 

Several examples are provided in the repository. 

* The four-column case is taken from the paper http://dx.doi.org/10.1016/j.compchemeng.2006.06.013;
* the eight-column case is taken from the paper http://dx.doi.org/10.1016/S0959-1524(01)00005-1; 
* the parameters of the ternary component in the binary separation example is made up from the eight-column case;
* the five-column case for ternary separation is taken from the paper http://dx.doi.org/10.1016/j.chroma.2011.09.015; 
* the ten-column case for ternary separation is taken from the paper http://dx.doi.org/10.1016/j.ces.2004.10.007.

By the way, the demonstration cases can directly be run coping them from the examples directory to the simulatedMovingBed directory, then changing the file name to getParameters.m. Apparently, the examples can also be modified by replacing your own models, operating conditions, and optimization routines. 

* How to write your own getParameter routine?
* How to adopt another type of equilibrium isotherm models?

# Documentation 

For more details of the CADET-SMB software, see the file doc.pdf in the repository and the doc.pdf in the master branch.

# Further Development 

SMB is actively developed. Hence, breaking changes and extensive restructuring may occur in any commit and release. For non-developers it is recommended to upgrade from release to release instead of always working with the most recent commit.
