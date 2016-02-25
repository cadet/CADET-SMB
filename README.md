![](https://github.com/modsim/CADET/blob/master/doc/logo/CADET-GitHub.png)

# CADET-SMB_Dynamic-Analog

The CADET-SMB is a comprehensive simulator for analysis and design of simulated moving bed (SMB) chromatograpy. And the variant, Dynamic Analog, is an improved version SMB which converges to the same cyclic steady state (CSS) with Full-simulation. However, it is little more time-consuming to compensate the benifits of obtaining the dynamic trajectory. This SMB simulation variant is regarded as the operator splitting.

# Introduction

There are various practical modes of carrying out industrial purification by preparative chromatography. The most straightforward and frequently used is cyclic batch elution chromatography, simulations of which are provided in the CADET framework (https://github.com/modsim/CADET). Another importent mode is counter-current chromatography, in which the fluid and the solid phases flow through the column in opposite directions. Since the continuous true moving bed (TMB) process is technically hard to implement, the cyclic simulated moving bed (SMB) is usually applied. In this repository, we offer an extension of the CADET framework for simulating SMB processes.
![The SMB scheme of binary separation](https://github.com/modsim/CADET-SMB/blob/master/doc/scheme_ternary.JPG)

CADET-SMB is developed at the Institute of Bio- and Geosciences 1 (IBG-1) of Forschungszentrum Jülich (FZJ) under supervision of Dr. Eric von Lieres. CADET-SMB uses the simulation engine of the CADET framework, which provides a fast and accurate solver for the general rate model (GRM) of packed bed liquid chromatography. Specific optimizers are also provided.

# Features

* It possesses the same features with the full-simulation, ternary separation and binary separation with arbitrary column configuration;
* It not only converges to the same statiionary state with the naive full-simulation, but also obtains the dynamic characteristics.
* Rearding the full-simulation, I refer you to https://github.com/modsim/CADET-SMB;

# Tutorial and Instructions

First, download the CADET software from https://github.com/modsim/CADET-SMB/releases, as CADET-SMB is based on the CADET simulator.

Regarding the installation of CADET,

* download the latest release for your platform;
* unzip the archive to your destination directory;
* start MATLAB;
* change to the unzipped CADET directory and run installCADET.m. (You can save the MATLAB path to avoid calling installCADET.m every time you restart MATLAB);
* Try one of the examples (e.g., examples/forward/loadWashElution.m) to check if everything works.
* Then, download the latest release of CADET-SMB from https://github.com/modsim/CADET-SMB/releases.

Regarding the installation of CADET-SMB,

* create a directory, simulatedMovingBed, in the CADET directory;
* unzip the CADET-SMB archive to that directory;
* Change the working directory in MATLAB to that directory; 
* To test a forward simulation, copy any getParameter routine from the examples/Forward folder to the simulatedMovingBed folder and change the name to getParameters.m; then run simulatedMovingBed.m.
* To test an optimization, copy any getParameter routine from the examples/Optimization folder to the simulatedMovingBed folder and also change the name to getParameters.m; then run SMBOptimization.m.

# Demenstration 

Several examples are provided in the repository. 

* The four-column case is taken from the paper http://dx.doi.org/10.1016/j.compchemeng.2006.06.013;
* the eight-column case is taken from the paper http://dx.doi.org/10.1016/S0959-1524(01)00005-1; 
* the parameters of the ternary component in the binary separation example is made up from the eight-column case;
* the five-column case for ternary separation is taken from the paper http://dx.doi.org/10.1016/j.chroma.2011.09.015; 
* the ten-column case for ternary separation is taken from the paper http://dx.doi.org/10.1016/j.ces.2004.10.007.

The demonstration cases can directly be run coping them from the examples directory to the simulatedMovingBed directory, then changing the file name to getParameters.m. 

Of course, the examples can be modified of replaced by own models, operating conditions, and optimization routines. 

* How to write your own getParameter routine?
* How to adopt another type of equilibrium isother models?

# Documentation 

For more details of the CADET-SMB software, see the file doc.pdf in the repository.

# Further Development 

SMB is actively developed. Hence, breaking changes and extensive restructuring may occur in any commit and release. For non-developers it is recommended to upgrade from release to release instead of always working with the most recent commit.
