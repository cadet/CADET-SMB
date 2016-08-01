![](https://github.com/modsim/CADET/blob/master/doc/logo/CADET-GitHub.png)

# CADET-SMB

CADET-SMB is a comprehensive simulator for analysis and design of simulated moving bed (SMB) chromatographic processes. It is developed at the Institute of Bio- and Geosciences 1 (IBG-1) of Forschungszentrum Jülich (FZJ) under supervision of Dr. Eric von Lieres. CADET-SMB uses the simulation engine of the CADET framework, which provides a fast and accurate solver for the general rate model (GRM) of packed bed liquid chromatography. 

# Advanced operator-splitting approach

This variant, advanced operator-splitting approach, is a developed version of the operator-splitting approach which is on the branch of the Dynamic_SMB of the GitHub repository. In those two operator-splitting approaches, we intend to approach the real flow pattern of the SMB unit such that the actual convergence trajectories are able to be obtain. Due to technical problems, we can not track the internal composition profile (the immediate concentration) in each column during processing, instead we only detect the concentration profile at the end of each column. Thus the simulation sequence in both standard_SMB and one-column analog approach is sequential, although clockwise in standard_SMB while counter-clockwise in one-column analog. However, the actually simulation sequence in reality should be simultaneous, like which is presented in following figure.

![](https://github.com/modsim/CADET-SMB/blob/Operator-splitting/doc/sequence.JPG)

It does not matter what the simulation sequence is, if only the cyclic steady state (CSS) is concerned in the simulations. Since both standard_SMB approach and one-column analog approach converge to the same CSS, though, via the different convergence dynamics (hereafter, I will just call the convergence dynamics the trajectory). However, if the actual trajectories, rather than the CSS, are concerned, neither standard_SMB approach or one-column analog approach can afford any more extra information. 

We need to resort to new approaches. The operator-splitting technique is used to overcome the sequential simulation order in both standard_SMB approach and the one-column analog approach. Operator splitting is a kind of mathematical terminology. By utilizing splitting of the simulations of columns in a SMB unit into several time sections, we could approach the real flow pattern as what has presented in the above figure. The definition of the time section is t_s/n, where n is the amount of time sections. 

![](https://github.com/modsim/CADET-SMB/blob/Operator-splitting/doc/operator_splitting.JPG)

*Brief demonstration*

And in the Dynamic_SMB branch of the GitHub, one of the variants of operator-splitting approaches has been proposed. The documentation of the fundamentals have been made in its doc.pdf, so there is no redundant illustrations here any more. However it possesses a drawback that some errors are introduced into the simulations. A new idea is proposed to neglect this drawback. The underlying idea here is that the very first column (or the starting simulation point) is always be simulated twice, and with some tricks to obtain both the true column state and the true outlet profile. 

![](https://github.com/modsim/CADET-SMB/blob/Operator-splitting/doc/advanced.JPG)

*schematic_1*

The attention is paid to the very first column, say, column in the zone I, during one switch time processing, see Fig.(schematic_1). At the very beginning, as there is no knowledge regarding the inlet profile of the column in one time section. Without losing the generality, the empty concentration profile is assumed. After one time section simulations, the actual inlet profile of the very first column is known which is gained from the outlet of the column in zone IV.

![](https://github.com/modsim/CADET-SMB/blob/Operator-splitting/doc/schematic.JPG)

*schematic_2*

It is worth mentioning that the column state after the simulation of the very first column within first time section is not reserved, instead only the outlet profile from such simulation is collected, see the part of Fig.(schematic_2) with assumed inlet profile. Because with the assumed inlet profile (also incorrect), the column is contaminated with the incorrect inlets. As seen from the Fig.(schematic_2), the original column state consists of both true s1 and s2 parts. With injecting the assumed inlet profile of t_s/n time units, the approximate s1 state is pushed out and replaced with s2 and wrong parts. However, the concentration pushed out is correct. So only the outlet concentration of the simulation within one time section is stored rather the column state. 

After one loop time section simulations, the true outlet profile of the column in zone IV is known, as shown in Fig.(schematic_1). If the cyclic steady state is arrived, this is the actual inlet profile of the column in zone I in next t_s switch interval. Then using the non-preserved column state and the "actual" column inlet profile, the column state that has not been reserved is updated, see the part of Fig.(schematic_2) with true inlet profile. The originally true column state, s1 and s2, is updated with also true s2 and s3 column state. It is the so-called saying that simulating the very first column twice in each time section. 

Similarly, in the simulations of the second time section, it still starts with the assumed empty inlet profile and correct column state (only after twice computing of the very first column), and at the end of the second time section, the very first column is updated with the true column state as mentioned before.



# Features

* Binary separation is available using four-zone scheme; Ternary separation can be achieved by using either integrated five-zone scheme, eight-zone, or cascade scheme;

![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_binary.JPG)
![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_ternary_5.JPG)
![](https://github.com/modsim/CADET-SMB/blob/master/doc/profile_ternary_8.JPG)

* In separations, arbitrary column configurations are available, in addition to basic column configurations such as 1-1-1-1, 2-2-2-2-2, 3-3-3-3, 4-4-4-4-4;

* It not only converges to the same cyclic steady state (CSS) with the naive standard_SMB approach and one-column analog approach, but also the actual convergence dynamics;

* It is also less time-consuming than the operator-splitting approach that is shown in the Dynamic_SMB branch.

![](https://github.com/modsim/CADET-SMB/blob/Operator-splitting/doc/time_comparison.JPG)

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
