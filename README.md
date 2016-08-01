![](https://github.com/modsim/CADET/blob/master/doc/logo/CADET-GitHub.png)

# CADET-SMB

CADET-SMB is a comprehensive simulator for analysis and design of simulated moving bed (SMB) chromatographic processes. It is developed at the Institute of Bio- and Geosciences 1 (IBG-1) of Forschungszentrum Jülich (FZJ) under supervision of Dr. Eric von Lieres. CADET-SMB uses the simulation engine of the CADET framework, which provides a fast and accurate solver for the general rate model (GRM) of packed bed liquid chromatography.

# Introduction

There are various practical modes of carrying out industrial purification by preparative chromatography. The most straightforward and frequently used one is the cyclic batch elution chromatography, simulations of that have been provided in the CADET framework (https://github.com/modsim/CADET). Another important mode is the counter-current chromatography, in which the fluid phase and the solid phase flow through the column in opposite directions. Since the continuous true moving bed (TMB) process is technically hard to implement, the periodically simulated moving bed (SMB) process is usually applied. In this repository, we offer an extension of the CADET framework, CADET-SMB, for simulating SMB chromatographic processes.

# Features

The introduction of features is divided into three parts, network setup, numerical computing, and inverse problem respectively. Correspondingly comprehensive diagrams will be provided to help you pinpoint the simulation approaches you intended, as there are several different methods that are integrated into the CADET-SMB for the sake of distinct simulation objects.

## Network setup

![](https://github.com/modsim/CADET-SMB/blob/master/doc/Network_setup.JPG)

Binary separation scenarios have been widely concerned in SMB chromatographic processes. In general, the typically four zone network setup (hence four-zone scheme) is, and merely, used for binary separations. However, the necessity of performing ternary separation has been raised. In ternary separation scenarios, two strategies can be opted. The first strategy is the cascade strategy, in which two separate SMB units are connected sequentially. The second strategy is an integrated SMB unit with five zones (hence five-zone scheme). Since both of them have advantages and disadvantages. So the main selecting criterion between those two network setups is the varying extent of affinity degree of the components to the solid beads. Moreover, arbitrary column configurations could be coped with rather than, say, 1-1-1-1 or 2-2-2-2-2. By the way, the code of the cascade scheme can be separately found in the Cascade branch of the GitHub.

## Numerical computing

![](https://github.com/modsim/CADET-SMB/blob/master/doc/Numerical_computing.JPG)

Specifically, there are two major types, overall four methods that are involved in the numerical computing part. They are standard SMB approach, one-column analog approach, operator-splitting approach, advanced operator-splitting approach respectively. Additionally, the corresponding code of the four methods can be found in the homonymous GitHub branch.

It is trivial to brief the standard SMB approach here. The details are shown in the documentation file. The one-column analog approach is derived from the works of Nadia Abunasser and Phillip Wankat, Jose Mota and Joao Araujo, which is aimed to speed up convergence to the cyclic steady state (CSS) from the point of the view of the simulation. It is also regarded as an attractive alternative to the standard SMB approach, notably in the optimization situations. It is necessary to emphasize that in this branch (CSS), more attentions are paid to the convergence speed. In the following branch (trajectory), concentration dynamics will be concerned.

Although in standard SMB approach and one-column analog approach the same correct CSS is converged, both of them cannot recur the actual dynamics of the concentration profiles that monitored at raffinate and extract ports (afterwards, I will just call these dynamics trajectories). It results from the concentration tracking manner, as in chromatography one can only detect the concentrations at the outlet end of the columns rather than immediate tracking. So the simulation sequence of a SMB unit is sequential, in reality it is actually simultaneous. Via operator-splitting approach in which the simulations of the columns are split into several time section intervals, the goal of obtaining correct trajectories is achieved. However in the operator-splitting approach, a big amount of time sections are required to eliminate the errors that are introduced into the system. In the advanced operator-splitting approach, the column at the very beginning of the simulations is simulated twice to compose both "true" column state and column outlet concentration profile. Hence we can approach the true trajectories by using operator-splitting and advanced operator-splitting approaches. Apparently, operator-splitting is more time-consuming than the advanced operator-splitting approach. There is a piece of evidence to back the statement that the true trajectories are obtained, by starting from different simulation starting points (e.g. feed, raffinate, extract, desorbent) the obtained trajectories coincide well with each other. In other words, the trajectory is unique. 


## Inverse problems

![](https://github.com/modsim/CADET-SMB/blob/master/doc/Inverse_problems.JPG)

Optimization of the CADET-SMB includes two parts, the optimization of decision variables and the optimization of the column configurations. The optimization of the column configuration part is still in development, so it is not presented yet. Regarding the optimization part of the decision variables, once the column configuration (it can be an arbitrary one, e.g. 1-2-3-1 or 2-1-3-2-4) is confirmed whenever by the structural optimization algorithm or mentally decision, the processing parameters (e.g. flow rates in different zones, switch time, column length) can be optimized with the user-defined objective function, which might aims to improve the productivity, or purity, or operation costs. As for the optimization algorithm, there are four options, a MATLAB build-in deterministic algorithm, the particle swarm optimization (PSO), the differential evolution (DE), and the Metropolis adjusted differential evolution (MADE). 

So far, the standard SMB approach and the one-column analog approach are combined with the optimization functionality, whereas the operator-splitting approach and the advanced operator-splitting approach are not. As the major goal of these two approaches is to approach the true trajectory.

## Feature list

![](https://github.com/modsim/CADET-SMB/blob/master/doc/diagram.JPG)

* Binary separation is available using four-zone scheme; ternary components in binary separation are also possible;

* Ternary separation can be achieved by using integrated five-zone scheme, eight-zone scheme, or cascade scheme; quaternary components in ternary separation are possible;

![](https://github.com/modsim/CADET-SMB/blob/master/doc/scheme.JPG)
*Four-zone scheme and five-zone scheme*

![](https://github.com/modsim/CADET-SMB/blob/master/doc/scheme_8.JPG)
*Eight-zone scheme with internal raffinate and internal extract*

![](https://github.com/modsim/CADET-SMB/blob/master/doc/cascade.JPG)
*Cascade scheme*


![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_binary.JPG)
![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_ternary_5.JPG)
*Chromatogram with four-zone and five-zone scheme*

![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_ternary_8.JPG)
*Chromatogram with eight-zone scheme*

* In both binary and ternary separations, arbitrary column configurations are available, in addition to basic column configurations such as 1-1-1-1, 2-2-2-2-2, 3-3-3-3, and 4-4-4-4-4;

* By using operator-splitting, the true dynamics of the concentration profiles (also, trajectories) can be reproduced. Also using one-column analog, the convergence speed to the cyclic steady state (CSS) could be accelerated.

* Continuous stirred tank reactor (CSTR) and dispersive plug flow reactor (DPFR) models can be placed before and after each column to account for dead volumes in pumps, tubing, and valves;
![](https://github.com/modsim/CADET-SMB/blob/master/doc/dead_volumes.JPG)
*The schematic of the dead volumes consideration*


* In binary separation, the ModiCon process is also available;

* MATLAB interface allows to monitor the dynamic characteristics of each column in the SMB unit;

* Optimization of decision variables for improving, e.g., productivity, purity, operating costs;

* Parameter estimation from experimental data will be implemented in future versions;

* Column models include transport dispersive model, equilibrium dispersive model, and general rate model;

* Wide range of standard equilibrium/isotherm models allow to simulate either pure component or multi-component/competitive behaviour;

* Further features of the CADET framework can be found at https://github.com/modsim/CADET.


# Dependency and Platforms

* Matlab(R2010b or higher);
* CADET (version 2.3.2 or later);
* platforms, please see the Dependencies section in the CADET wiki.

# Tutorial and Instructions

First, download the CADET software from https://github.com/modsim/CADET-SMB/releases, as CADET-SMB is based on the CADET simulator.
Then, download the latest release of CADET-SMB from https://github.com/modsim/CADET-SMB/releases.

## Installation of the CADET

* download the latest release for your platform;
* unzip the archive to your destination directory;
* start MATLAB;
* change directory to the unzipped CADET directory and run installCADET.m (you can save the MATLAB path to avoid calling installCADET.m every time you restart MATLAB);
* Try one of the examples (e.g., examples/forward/loadWashElution.m) to check if everything works.
 
## Installation of the CADET-SMB,

* create a directory, simulatedMovingBed, in your unzipped CADET directory;
* unzip the CADET-SMB archive to the simulatedMovingBed directory;
* Change the working directory to the simulatedMovingBed directory and run isSMBupdateAvailable.m script (Along side checking the existence of the newest version, it also attach the current path the MATLAB path); 
* To test a forward simulation, copy any getParameter_something script from the examples/Forward folder to the simulatedMovingBed folder and change the name to getParameters.m; then run simulatedMovingBed.m.
* To test an optimization, copy any getParameter_something script from the examples/Optimization folder to the simulatedMovingBed folder and also change the name to getParameters.m; then run SMBOptimization.m.

# Demonstration 

Several examples are provided in the repository. 

* The four-column case for binary separation is taken from the paper http://dx.doi.org/10.1016/j.compchemeng.2006.06.013;
* the eight-column case for binary separation is taken from the paper http://dx.doi.org/10.1016/S0959-1524(01)00005-1; 
* the parameters of the ternary component in the binary separation example is made up from the eight-column case artificially;
* the five-column case for ternary separation is taken from the paper http://dx.doi.org/10.1016/j.chroma.2011.09.015; 
* the ten-column case for ternary separation is taken from the paper http://dx.doi.org/10.1016/j.ces.2004.10.007.

By the way, the demonstration cases can directly be run by coping them from the examples directory to the simulatedMovingBed directory, then changing the file name to getParameters.m. Apparently, the examples can also be modified by replacing your own models, operating conditions, and optimization routines. 

* How to write your own getParameters routine?
* How to adopt other isotherm models?

# Documentation 

For more details of the CADET-SMB software, see the file doc.pdf in the repository.

# Further Development 

CADET-SMB is actively developed. Hence, breaking changes and extensive restructuring may occur in any commit and release. For non-developers it is recommended to upgrade from release to release instead of always working with the most recent commit.
