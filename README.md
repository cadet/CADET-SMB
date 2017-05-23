![](https://github.com/modsim/CADET/blob/master/doc/logo/CADET-GitHub.png)

# CADET-SMB
CADET-SMB is a comprehensive simulator for analysis and design of simulated moving bed (SMB) chromatographic processes. It is developed at the Institute of Bio- and Geosciences 1 (IBG-1) of Forschungszentrum Jülich (FZJ) under supervision of Dr. Eric von Lieres. CADET-SMB uses the simulation engine of the CADET framework, which provides a fast and accurate solver for the general rate model (GRM) of packed bed liquid chromatography.

# Introduction
There are various practical modes of preparative chromatography. Cyclic batch elution chromatography is most frequently applied. In counter-current chromatography, the fluid and solid phase are moved through the column in opposite directions. Since the true moving be (TMB) process is technically hard to implement, the simulated moving bed (SMB) process is usually applied. In this repository, we offer an extension of the CADET framework, CADET-SMB, for simulating SMB chromatographic processes.

# Branch: cascade of STD-FPI
For more general introduction, please see the main branch, in which we introduce the standard version of fixed point iteration (STD-FPI) method.

In this branch, we introduce the cascade of two four-zone SMB in a row using STD-FPI method to achieve the goal of ternary separations. 

![](https://github.com/modsim/CADET-SMB/blob/Cascade/doc/scheme.JPG)

*Cascade scheme of two four-zone SMB in a row*

# Detailed feature list
* Ternary separation is available by using the cascade scheme;

![](https://github.com/modsim/CADET-SMB/blob/Cascade/doc/profile_cascade.JPG)

*The respective chromatograms of the cascade system*

* Arbitrary column configurations are available, in addition to basic column configurations such as 1-1-1-1, 2-2-2-2, 3-3-3-3, 4-4-4-4;

* Continuous stirred tank reactor (CSTR) and dispersive plug flow reactor (DPFR) models can be placed before and after each column to account for dead volumes in pumps, tubing, and valves;

![](https://github.com/modsim/CADET-SMB/blob/master/doc/dead_volumes.JPG)

* MATLAB interface allows to monitor the dynamic characteristics of each column in the SMB unit;

* Column models include transport dispersive model, equilibrium dispersive model, and general rate model;

* Wide range of standard equilibrium/isotherm models allow to simulate either pure component or multi-component/competitive behaviour;

* Further features of the CADET framework can be found at https://github.com/modsim/CADET.

# Dependency and Platforms

* Matlab(R2010b or higher);
* CADET (version 2.3.2);
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
