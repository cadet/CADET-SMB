![](https://github.com/modsim/CADET/blob/master/doc/logo/CADET-GitHub.png)

# CADET-SMB
The CADET-SMB is a comprehensive simulator used for design and simulation of chromatographic Simulated Moving Bed processes. 

# Introduction

There are two important, practical modes of carrying out industrial purification by preparative chromatography. The most straightforward and the most frequently used of these is the cyclic batch elution chromatography. The simulation of the first mode is available in CADET (https://github.com/modsim/CADET). The second mode available is counter-current chromatography, in which the fluid and the solid phase flow through the column in opposite directions. In principle, there are two possibilities of this mode, the continuous true moving bed and the cyclic simulated moving bed method (SMB). In practical, only the latter method is widespread. In this repository, we offer a open-source software to carry out the simulation of Simulated Moving Bed (SMB) processes.
![The SMB scheme of binary separation](https://github.com/modsim/CADET-SMB/blob/master/doc/scheme_binary.JPG)

The CADET-SMB is developed at the Institute of Bio- and Geosciences 1 (IBG-1) of Forschungszentrum Jülich (FZJ) under supervision of Dr. Eric von Lieres. CADET-SMB comprises tow parts, CADET which is a fast and accurate solver for the General Rate Model (GRM) of packed bed liquid chromatography, and SMB which is the simulator for simulated moving bed. SMB also involves optimizer that covers a wide range of variants.

# Features

* Ternary separation is available, using Five-Zone technique; Quaternary components in ternary separation is possible;
![The SMB scheme of ternary separation](https://github.com/modsim/CADET-SMB/blob/master/doc/scheme_ternary.JPG)
![Axial concentration profile in the ternary separation](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_ternary.JPG)
* In ternary separation, 1-1-1-1-1, 2-2-2-2-2, 3-3-3-3-3, and 4-4-4-4-4 configurations are accessible;
* Binary separation is available; Ternary components in binary separation is also possible;
* In binary separation, 1-1-1-1, 2-2-2-2, 3-3-3-3, 4-4-4-4 column configurations are available;
* In binary separation, the using of ModiCon tech is possible;
![Axial concentration profile in the binary separation](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_binary.JPG)
* MATLAB interface, you are allowed to monitor the dynamic characteristics of each column;
* Optimization of decision variables to gain benefits in productivity, purity, operating costs;
* Fit model against experimental data will come later;
* Flexible model that permits the study of other SMB variants, such as One-Column Analog, Dynamic Analog, which are focused on the different objects. However, they are not published in GitHub now;
* Kinetic models of CADET including transport dispersive, equilibrium dispersive and general rate model (GRM);
* Wide range of standard equilibrium/isotherm models that allow either pure component or multi-component/competitive behaviour;
* As for the more features of CADET, I refer you to the website, https://github.com/modsim/CADET.

# Tutorial and Instructions

First of all ,please download the CADET software form https://github.com/modsim/CADET-SMB/releases, as CADET-SMB is based on the CADET software.

Regarding the installation of CADET,

* Download the latest release for your platform;
* Unzip the archive to your destination directory;
* Fire up MATLAB;
* Change to the unzipped cadet directory and run installCADET.m. Additionally you can save the MATLAB path to avoid calling installCADET.m every time you restart MATLAB;
* Try one of the examples (say, examples/forward/loadWashElution.m) to check if everything works.
* Then, download the latest release of SMB from https://github.com/modsim/CADET-SMB/releases.

Regarding the installation of SMB,

* Make a directory, simulatedMovingBed, under the directory of the CADET;
* Unzip the archive to the new directory, simulatedMovingBed;
* Change the directory to simulatedMovingBed 
* For the forward simulation, copy one getParameter routine from the examples/Forward to the directory of simulatedMovingBed, and change the name to getParameters.m; Then run the routine, simulatedMovingBed.m.
* For the optimization case, copy one getParameter routine from the examples/Optimization to the directory of simulatedMovingBed, and also change the name to getParameters.m; Then run the routine, SMBOptimization.m.

# Demenstration 

As for the several demonstrated cases in the examples repository are both laboratory cases. 

* While the four-column case is from the paper http://www.sciencedirect.com/science/article/pii/S009813540600192X;
* and the eight-column case is from the paper http://www.sciencedirect.com/science/article/pii/S0959152401000051; 
* But the data of case for the ternary components in binary separation is totally fake, which is extended from the eight-column case;
* The five-column case for the ternary separation is from the paper, http://www.sciencedirect.com/science/article/pii/S002196731101363X; 
* The ten-column case for the ternary separation is from the paper, http://www.sciencedirect.com/science/article/pii/S0009250904007742.

You can directly use the demonstration cases by coping it from the examples directory into the simulatedMovingBed directory, then changing the file name to getParameters.m. 

You can also write your own parameter routines by refering the getParameters.m. 

* How to write your own getParameter routine?
* How to adopt another type of equilibrium isother models?

# Documentation 

For details of the software, I refer you to the file, doc.pdf, in the repository.

# Keeping Devekopment 

Since SMB is actively developed. So breaking changes and extensive restructuring may occur in any commit and release. For non-developers it is recommended to upgrade from release to release instead of always working with the most recent commit.
