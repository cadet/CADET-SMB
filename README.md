

# CADET-SMB
The Simulated Moving Bed (SMB) chromatography for separation of target compounds, such as Fructose and Glucose.

# Introduction

There are two important, practical modes of carrying out industrial purification by preparative chromatography. The most straightforward and the most frequently used of these is the cyclic batch elution chromatography. The simulation of the first mode is available in CADET (https://github.com/modsim/CADET). The second mode available is counter-current chromatography, in which the fluid and the solid phase flow through the column in opposite directions. In principle, there are two possibilities of this mode, the continuous true moving bed and the cyclic simulated moving bed method (SMB). In practical, only the latter method is widespread. In this repository, we offer a open-source code to carry out the simulation of SMB.

The CADET-SMB is developed at the Institute of Bio- and Geosciences 1 (IBG-1) of Forschungszentrum Jülich (FZJ) under supervision of Dr. Eric von Lieres. CADET-SMB comprises tow parts, CADET which is a fast and accurate solver for the General Rate Model (GRM) of packed bed liquid chromatography, and SMB which is the simulator for simulated moving bed. SMB also involves optimizer that covers a wide range of variants.

# Features

* Binary separation is available now, the ternary separation will come later;
* Ternary components in binary separation is also possible;
* 1-1-1-1, 2-2-2-2, 3-3-3-3, 4-4-4-4 column configurations are available;
* MATLAB interface, you are allowed to monitor the dynamic characteristic of each column.
* Optimization of decision variables.
* As for the features of CADET, I refer you to the website, https://github.com/modsim/CADET;

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

* Make a directory, say, simulatedMovingBed, in the directory of the CADET;
* Unzip the archive to the new directory, in my case simulatedMovingBed;
* Change the directory to simulatedMovingBed and run simulatedMovingBed.m.

# Documentation 

I refer you to the file, doc.pdf, in the repository.

# Keeping Devekopment 

Since SMB is actively developed. So breaking changes and extensive restructuring may occur in any commit and release. For non-developers it is recommended to upgrade from release to release instead of always working with the most recent commit.
