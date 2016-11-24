![](https://github.com/modsim/CADET/blob/master/doc/logo/CADET-GitHub.png)

# CADET-SMB
CADET-SMB is a comprehensive simulator for analysis and design of simulated moving bed (SMB) chromatographic processes. It is developed at the Institute of Bio- and Geosciences 1 (IBG-1) of Forschungszentrum Jülich (FZJ) under supervision of Dr. Eric von Lieres. CADET-SMB uses the simulation engine of the CADET framework, which provides a fast and accurate solver for the general rate model (GRM) of packed bed liquid chromatography. 

# Introductionc
There are various practical modes of preparative chromatography. Cyclic batch elution chromatography is most frequently applied. In counter-current chromatography, the fluid and solid phases are moved through the column in opposite directions. Since the true moving bed (TMB) process is technically hard to implement, the simulated moving bed (SMB) process is usually applied. In this repository, we offer an extension of the CADET framework, CADET-SMB, for simulating SMB chromatographic processes.

# Branch: LAW-OPS
For more general introduction, please see the main branch, in which we introduce the standard version of fixed point iteration (STD-FPI) method. 

In this branch, we introduce the lag-aware operator-splitting (LAW-OPS) method. The underlying idea to have this method is that only the chromatograms under cyclic steady state (CSS) are provided in previous literatures, rather than the whole trajectories from starting time to the CSS. This is resulted from the "lag" problem in conventionally SMB simulations. Thus LAW-OPS is proposed to overcome this issue. 

![](https://github.com/modsim/CADET-SMB/blob/Operator-splitting/doc/flow_pattern.JPG)

*The outlet of the cell 1 in previous columns should be transferred immediately to the inlet of the cell 1 in latter columns*

# Detailed feature list
* Binary separation is available using four-zone scheme; 

![](https://github.com/modsim/CADET-SMB/blob/master/doc/scheme_binary.JPG)
![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_binary.JPG)

*Four zone scheme for binary separations and the chromatogram of the four-zone SMB*

* Ternary separation is available by using the cascade scheme, the integrated five-zone or eight-zone schemes;

![](https://github.com/modsim/CADET-SMB/blob/master/doc/scheme_cascade.JPG)
![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_cascade.JPG)

*Cascade scheme for ternary separations and the respective chromatograms of the cascade system*

![](https://github.com/modsim/CADET-SMB/blob/master/doc/scheme_ternary_8.JPG)
![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_ternary_8.JPG)

*Eight-zone scheme for ternary separations and the chromatogram of the eight-zone scheme*

![](https://github.com/modsim/CADET-SMB/blob/master/doc/scheme_ternary_5.JPG)
![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_ternary_5.JPG)

*Five-zone scheme for ternary separations and the chromatogram of the five-zone scheme*

* In both binary and ternary separations, arbitrary column configurations are available, in addition to basic column configurations such as 1-1-1-1, 2-2-2-2-2, 3-3-3-3, 4-4-4-4-4;

* We provide not only the CSS information like what have shown above, as well as the trajectory information by using LAW-OPS method.

![](https://github.com/modsim/CADET-SMB/blob/Operator-splitting/doc/trajectory_extract.JPG)
![](https://github.com/modsim/CADET-SMB/blob/Operator-splitting/doc/trajectory_raffinate.JPG)

*The trajectories from LAW-OPS. Left side is from extract port, right side is from raffinate port*

* Continuous stirred tank reactor (CSTR) and dispersive plug flow reactor (DPFR) models can be placed before and after each column to account for dead volumes in pumps, tubing, and valves;

![](https://github.com/modsim/CADET-SMB/blob/master/doc/dead_volumes.JPG)

* MATLAB interface allows to monitor the dynamic characteristics of each column in the SMB unit;

* Optimization of decision variables for improving, e.g., productivity, purity, operating costs;

* Parameter estimation from experimental data will be implemented in future versions;

* Column models include transport dispersive model, equilibrium dispersive model, and general rate model;

* Wide range of standard equilibrium/isotherm models allow to simulate either pure component or multi-component/competitive behaviour;

* Further features of the CADET framework can be found at https://github.com/modsim/CADET.


# Dependency and Platforms
* Matlab (R2010b or higher);
* CADET (version 2.3.2);
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
