![](https://github.com/modsim/CADET/blob/master/doc/logo/CADET-GitHub.png)

# CADET-SMB

CADET-SMB is a comprehensive simulator for analysis and design of simulated moving bed (SMB) chromatographic processes. It is developed at the Institute of Bio- and Geosciences 1 (IBG-1) of Forschungszentrum Jülich (FZJ) under supervision of Dr. Eric von Lieres. CADET-SMB uses the simulation engine of the CADET framework, which provides a fast and accurate solver for the general rate model (GRM) of packed bed liquid chromatography. 

# operator-splitting approach

This variant, operator-splitting approach, is a referee of the approach, advanced operator-splitting, which is on the branch of the operator-splitting of the GitHub repository. In both those two approaches, we intend to approach the real flow pattern of the SMB unit such that the actual convergence trajectories are able to be obtain. And Due to technical problems, we can not track the internal composition profile (the immediate concentration) in each column during processing, instead we only detect the concentration profile at the end of each column. Thus the simulation sequence in both standard_SMB and one-column analog approach is sequential, although clockwise in standard_SMB while counter-clockwise in one-column analog. However, the actually simulation sequence in reality should be simultaneous, like which is presented in following figure.

![](https://github.com/modsim/CADET-SMB/blob/Dynamic_SMB/doc/sequence.JPG)

It does not matter what the simulation sequence is, if only the cyclic steady state (CSS) is concerned in the simulations. Since both standard_SMB approach and one-column analog approach converge to the same CSS, though, via the different convergence dynamics (hereafter, I will just call the convergence dynamics the trajectory). However, if the actual trajectories, rather than the CSS, are concerned, neither standard_SMB approach or one-column analog approach can afford any more extra information. We need to resort to new approaches. The operator-splitting technique is used to overcome the sequential simulation order in both standard_SMB approach and the one-column analog approach.


Operator splitting is a kind of mathematical terminology. By utilizing splitting of the simulations of columns in a SMB unit into several time sections, we could approach the real flow pattern as what has presented in the above figure. The definition of the time section is t_s/n, where n is the amount of time sections. 

![](https://github.com/modsim/CADET-SMB/blob/Dynamic_SMB/doc/operator_splitting.JPG)

*Brief demonstration*

For simplicity, the amount of three is used. In considering simulations of one switch time interval, the simulations of all red time sections are implemented sequentially rather than the whole columns. Subsequently, the green time sections and blue time sections. Although the columns are still sequentially calculated, the really simultaneous flow pattern is able to be reproduced if the amount of time sections is big enough.

Specifically, assuming the simulations start from the desorbent node, the outlet of the column in zone I within first time section interval (t_s/3) can be calculated given the inlet profile. In other words, the red piece of the column in zone I is pushed out to the first interval chromatogram. Then we switch to simulate the column in zone II with the given inlet profile from the simulation the column in zone I. Similarly, the red piece in zone II is pushed out to the first interval chromatogram. After all the red pieces are pushed out, the simulations of first time section have been accomplished. 

Now there is a critical point that in the pushing-out of the green pieces of the column in zone I, what is the inlet profile? If the outlet profile stored in the zone IV chromatogram is adopted as the inlet profile, some errors rather than noise are introduced into the simulations (this will be illustrated later). Again, the green and blue pieces could be pushed out with given inlet profiles, and being stored in respectively chromatograms. 

The introduced errors can be eliminated by increasing the amount of the time sections in this variant. In the following figures (a) there is 20 time sections in each column, while 50 time sections in (b), 100 in (c). The introduced errors could be eliminated by increasing the amount drastically, correspondingly the price of the computational time what is the major disadvantage.

![](https://github.com/modsim/CADET-SMB/blob/Dynamic_SMB/doc/interval_20.JPG)
![](https://github.com/modsim/CADET-SMB/blob/Dynamic_SMB/doc/interval_50.JPG)
![](https://github.com/modsim/CADET-SMB/blob/Dynamic_SMB/doc/interval_100.JPG)

*The effects of the introduced errors which could be eliminated by the increasing of the amount of time sections*



# Features

* Binary separation is available using four-zone scheme; Ternary separation can be achieved by using either integrated five-zone scheme or cascade scheme;

![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_binary.JPG)
![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_ternary.JPG)

* In separations, arbitrary column configurations are available, in addition to basic column configurations such as 1-1-1-1, 2-2-2-2-2, 3-3-3-3, 4-4-4-4-4;

* It not only converges to the same cyclic steady state (CSS) with the naive standard_SMB approach and one-column analog approach, but also the actual convergence dynamics. Correspondingly, more time-consuming;

![](https://github.com/modsim/CADET-SMB/blob/Dynamic_SMB/doc/time_comparison.JPG)

* Continuous stirred tank reactor (CSTR) and dispersive plug flow reactor (DPFR) models can be placed before and after each column to account for dead volumes in pumps, tubing, and valves;

![](https://github.com/modsim/CADET-SMB/blob/master/doc/dead_volumes.JPG)

* MATLAB interface allows to monitor the dynamic characteristics of each column in the SMB unit;

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
