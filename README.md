![](https://github.com/modsim/CADET/blob/master/doc/logo/CADET-GitHub.png)

# CADET-SMB

CADET-SMB is a comprehensive simulator for analysis and design of simulated moving bed (SMB) chromatographic processes. It has been developed at the Institute of Bio- and Geosciences 1 (IBG-1) of Forschungszentrum Jülich (FZJ) under supervision of Dr. Eric von Lieres. CADET-SMB uses the simulation engine of the CADET framework, which provides a fast and accurate solver for the general rate model (GRM) of packed bed liquid chromatography. Further details can be found in the following publications:

* He, Q.-L.; Leweke, S.; von Lieres, E.: [Efficient numerical simulation of simulated moving bed chromatography with a single-column solver](http://doi.org/10.1016/j.compchemeng.2017.12.022), Computers and Chemical Engineering 111 (2018), 183–198.
* Leweke, S.; von Lieres, E.: [Chromatography Analysis and Design Toolkit (CADET)](http://doi.org/10.1016/j.compchemeng.2018.02.025), Computers and Chemical Engineering 113 (2018), 274–294.

# Important Note

At the time when CADET-SMB was developed, the CADET engine was limited to a single chromatography column. The CADET engine now supports strongly coupled networks of unit operations, while CADET-SMB is based on weak coupling. In addition, the model family of the CADET engine has been enlarged, including CSTR and DPFR units (and many other). Hence, **we stongly recommend using the CADET engine for most applications** (https://github.com/modsim/CADET). CADET-SMB can still be beneficial (faster) when only the cyclic steady state (CSC) is of interest and the one-column analog can be applied. CADET-SMB has been used for verifying the network capabilities of the CADET engine. CADET-SMB also demonstrates the feasibility of the lag-aware operator splitting approach, which we are proud to have found.

# Introduction

There are various practical modes of preparative chromatography. Cyclic batch elution chromatography is most frequently applied, and an efficient simulator is provided in the CADET framework (https://github.com/modsim/CADET). In counter-current chromatography, the fluid and solid phases are moved through the column in opposite directions. Since the true moving bed (TMB) process is technically hard to implement, the simulated moving bed (SMB) process is usually applied. In this repository, we offer an extension of the CADET framework, CADET-SMB, for simulating SMB chromatographic processes. CADET-SMB is generally capable of simulating arbitrary column networks with and withoput closed loops. In addition to chromatography columns, the CSTR and DPFR models can be included in the network, e.g. to account for hold-up volumes.

# Features

Code features are organized into network setup, numerical methods, and inverse problems.

## Network topology

SMB chromatography has originally been developed for binary (two components) separations. This is typically achieved by using four distinct zones with one or more columns each. Later, SMB variants have been developed for ternary (three components) separations. Two major strategies can be distinguished, both of which have advantages and disadvantages: a) sequential cascade of two conventional SMB units with eight zones in total, and b) integrated SMB units with eight or down to five zones. Moreover, CADET-SMB can be set-up with arbitrary column configurations, e.g., for simulating multicolumn counter-current solvent gradient purification (MCSGP).

## Numerical methods

CADET-SMB provides two classes of numerical solution approaches: a) fixed point iteration (FPI) for computing the cyclic steady state (CSS) of an SMB unit, and b) operator splitting (OSP) for computing the dynamic trajectory (DTR) from any initial system state into the CSS. Two variants are implemented for each approach, standard versions (STD-FPI, STD-OSP) and alternatives with significantly improved numerical efficiency, namely fixed point iteration for the one-column analog (OCA-FPI) and lag-aware operator splitting (LAW-OSP). The improved performance of these numerical methods can be particularly useful in optimization settings. Details on all four approaches can be found in the documentation, and a scientific journal publication is currently in preparation.

## Inverse problems

In SMB chromatography, both the operating conditions (column dimensions, flow rates, switch times) and the column configurations (network topology) can be optimized, leading to a mixed-integer nonlinear programming problem. However, optimization of the (discrete) column configuration is not yet implemented in CADET-SMB. For any given network topology, the operating conditions can be optimized with respect to user-specified objectives, e.g., purity, yield, cost. As these objectives are to be optimized in CSS, only the FPI approach is supported. Available search strategies include standard MATLAB functionality, particle swarm optimization (PSO), differential evolution (DE), Markov Chain Monte Carlo (DRAM version) and Metropolis adjusted differential evolution (MADE).

## Detailed feature list

* Binary separation is available using the four-zone scheme;

![](https://github.com/modsim/CADET-SMB/blob/master/doc/scheme_binary.JPG)
![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_binary.JPG)

*Four zone scheme for binary separations and the chromatogram of the four-zone SMB*

* Ternary separation is available using the cascade scheme, the integrated eight-zone or five-zone schemes;

![](https://github.com/modsim/CADET-SMB/blob/master/doc/scheme_cascade.JPG)
![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_cascade.JPG)

*Cascade scheme for ternary separations and the respective chromatograms of the cascade system*

![](https://github.com/modsim/CADET-SMB/blob/master/doc/scheme_ternary_8.JPG)
![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_ternary_8.JPG)

*Eight-zone scheme for ternary separations and the chromatogram of the eight-zone scheme*

![](https://github.com/modsim/CADET-SMB/blob/master/doc/scheme_ternary_5.JPG)
![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_ternary_5.JPG)

*Five zone scheme for ternary separations and the chromatogram of the five-zone scheme*

* In both binary and ternary separations, arbitrary column configurations are available, in addition to basic column configurations such as 1-1-1-1, 2-2-2-2-2, 3-3-3-3, and 4-4-4-4;

* The dynamics of the concentration profiles (system trajectories) can be computed by operator splitting (OSP);

* The convergence of the fixed point iteration (FPI) to the cyclic steady state (CSS) can be accelerated by the one-column analog (OCA) approach;

* Continuous stirred tank reactor (CSTR) and dispersive plug flow reactor (DPFR) models can be placed before and after each column to account for dead volumes in pumps, tubing, and valves;

![](https://github.com/modsim/CADET-SMB/blob/master/doc/dead_volumes.JPG)

*The schematic of the dead volumes consideration*

* In binary separation, the ModiCon process is also available;

* The MATLAB interface allows to monitor the dynamic characteristics of each column in the SMB unit;

* Optimization of decision variables for improving, e.g., productivity, purity, operating costs;

* Columns are modeled using the general rate model (GRM);

* A wide range of binding models is avaliable for simulating both single-component and competitive multi-component binding;

* Further features of the CADET framework can be found at https://github.com/modsim/CADET.

# Dependencies and Platforms

* Matlab (R2010b or higher);
* CADET (version 3.0.0);
* Linux, Windows, Mac OS.

# Tutorial and Instructions

* First, download CADET from https://github.com/modsim/CADET-SMB/releases, as CADET-SMB is based on the CADET simulator.
* Then, download CADET-SMB from https://github.com/modsim/CADET-SMB/releases.

## Installation of CADET

* Unzip the archive to your destination directory;
* start MATLAB;
* change directory to the unzipped CADET directory and run installCADET.m (you can save the MATLAB path to avoid calling installCADET.m every time you restart MATLAB);
* Try one of the examples (e.g., examples/forward/loadWashElution.m) to check if everything works.
 
## Installation of CADET-SMB

* Create a directory, simulatedMovingBed, in your unzipped CADET directory;
* Unzip CADET-SMB archive to the simulatedMovingBed directory;
* Change the working directory to the simulatedMovingBed directory and run isSMBupdateAvailable.m script; 
* Test a forward simulation by copying any getParameter_something script from the examples/forward folder to the simulatedMovingBed folder and change the name to getParameters.m and then run simulatedMovingBed.m;
* Test an optimization by copying any getParameter_something script from the examples/Optimization folder to the simulatedMovingBed folder and change the name to getParameters.m and then run SMBOptimization.m.

# Demonstration 

Several examples are provided in the repository. 

* Four-column case for binary separation is taken from the paper http://dx.doi.org/10.1016/j.compchemeng.2006.06.013;
* Eight-column case for binary separation is taken from the paper http://dx.doi.org/10.1016/S0959-1524(01)00005-1; 
* Five-column case for ternary separation is taken from the paper http://dx.doi.org/10.1016/j.chroma.2011.09.015; 
* Ten-column case for ternary separation is taken from the paper http://dx.doi.org/10.1016/j.ces.2004.10.007.

The demonstration cases can directly be run by coping them from the examples directory to the simulatedMovingBed directory, then changing the file name to getParameters.m. The examples can be adapted to other use cases by changing the column models, operating conditions, and optimization routines. 

# Documentation 

The methods are described in the following publications:

* He, Q.-L.; Leweke, S.; von Lieres, E.: [Efficient numerical simulation of simulated moving bed chromatography with a single-column solver](http://doi.org/10.1016/j.compchemeng.2017.12.022), Computers and Chemical Engineering 111 (2018), 183–198.
* Leweke, S.; von Lieres, E.: [Chromatography Analysis and Design Toolkit (CADET)](http://doi.org/10.1016/j.compchemeng.2018.02.025), Computers and Chemical Engineering 113 (2018), 274–294.

For more details of the CADET-SMB software, see the file doc.pdf in the repository.

# Further Development 

CADET-SMB is currently not actively developed.
