![](https://github.com/modsim/CADET/blob/master/doc/logo/CADET-GitHub.png)

# CADET-SMB
The CADET-SMB is a comprehensive simulator for analysis and design of simulated moving bed (SMB) chromatograpy. 

# Introduction

There are various practical modes of carrying out industrial purification by preparative chromatography. The most straightforward and frequently used is cyclic batch elution chromatography, simulations of which are provided in the CADET framework (https://github.com/modsim/CADET). another importent mode is counter-current chromatography, in which the fluid and the solid phases flow through the column in opposite directions. Since the continuous true moving bed (TMB) process is technically hard to implement, the cyclic simulated moving bed (SMB) is usually applied. In this repository, we offer an extension of the CADET framework for simulating SMB processes.
![The SMB scheme of binary separation](https://github.com/modsim/CADET-SMB/blob/master/doc/scheme_binary.JPG)

CADET-SMB is developed at the Institute of Bio- and Geosciences 1 (IBG-1) of Forschungszentrum Jülich (FZJ) under supervision of Dr. Eric von Lieres. CADET-SMB uses the simulation engine of the CADET framework, which provides a fast and accurate solver for the general rate model (GRM) of packed bed liquid chromatography. Specific optimizers are also provided.

# Features

* Ternary separation is available, using five-zone technique; quaternary components in ternary separation are possible;
![The SMB scheme of ternary separation](https://github.com/modsim/CADET-SMB/blob/master/doc/scheme_ternary.JPG)
![Axial concentration profile in the ternary separation](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_ternary.JPG)
* In ternary separation, 1-1-1-1-1, 2-2-2-2-2, 3-3-3-3-3, and 4-4-4-4-4 configurations are currently provided;
* Binary separation is available; ternary components in binary separation are also possible;
* In binary separation, 1-1-1-1, 2-2-2-2, 3-3-3-3, 4-4-4-4 column configurations are currently provided;
* In binary separation, the ModiCon process is available;
![Axial concentration profile in the binary separation](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_binary.JPG)
* MATLAB interface allows to monitor the dynamic characteristics of each column;
* Optimization of decision variables for improving, e.g., productivity, purity, operating costs;
* Parameter estimation form experimental data will be implemented in future versions;
* Column models include transport dispersive, equilibrium dispersive and general rate model (GRM);
* Wide range of standard equilibrium/isotherm models that allow to simulate either pure component or multi-component/competitive behaviour;
* Further features of the CADET framework can be found at https://github.com/modsim/CADET.

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

* The four-column case is taken from the paper http://www.sciencedirect.com/science/article/pii/S009813540600192X;
* the eight-column case is taken from the paper http://www.sciencedirect.com/science/article/pii/S0959152401000051; 
* the parameters of the ternary component in the binary separation example is made up from the eight-column case;
* the five-column case for ternary separation is taken from the paper http://www.sciencedirect.com/science/article/pii/S002196731101363X; 
* the ten-column case for ternary separation is taken from the paper http://www.sciencedirect.com/science/article/pii/S0009250904007742.

The demonstration cases can directly be run coping them from the examples directory to the simulatedMovingBed directory, then changing the file name to getParameters.m. 

Of course, the examples can be modified of replaced by own models, operating conditions, and optimization routines. 

* How to write your own getParameter routine?
* How to adopt another type of equilibrium isother models?

# Documentation 

For more details of the CADET-SMB software, see the file doc.pdf in the repository.

# Further Devekopment 

SMB is actively developed. Hence, breaking changes and extensive restructuring may occur in any commit and release. For non-developers it is recommended to upgrade from release to release instead of always working with the most recent commit.
