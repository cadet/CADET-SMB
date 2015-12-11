The fourth version of the SMB

In this version, the four_column case and eight_column case were integrated into one code. So it is general. In comparison to the regular SMB (which was presented in the master branch of the Github), a novel single-column setup for experimentally reproducing the steady periodic behaviour of simulated countercurrent multicolumn chromatography was presented in this branch.

The data of case I is from the paper, Numerical Method for Accelaration Calculation of Cyclic Steady State of ModiCon-SMB-Processes. It involves both four_column case and eight_column case.

k_a & k_d : [5.72 7.7]  (Linear isotherm);
d_{col}   : 0.02 m;
L         : 0.25 m;
porosity  : 0.83;
d_{par}   : not given, in my case, 0.0005 m;
Switching time        : 180 s
concentration of feed : 0.55 g/l
Recycle flow rate     : 9.62e-7 m^3/s
Feed flow rate        : 0.98e-7 m^3/s
Desorbent flow rate   : 1.96e-7 m^3/s
Extract flow rate     : 1.54e-7 m^3/s
Raffinate flow rate   : 1.40e-7 m^3/s

While the data of case II is from the paper, Model-based Control of a Simulated Moving Bed Chromatographic Process for the Separation of Fructose and Glucose. It is merely a eight_column case.

k_a & k_d : [0.54 0.28]  (Linear isotherm);
d_{col}   : 2.6cm;
L         : 53.6cm;
porosity  : 0.38;
d_{par}   : 0.325cm;
Switching time        : 1552 s
concentration of feed : 0.5 g/cm^3
Recycle flow rate     : 0.1395 cm^3/s
Feed flow rate        : 0.02 cm^3/s
Desorbent flow rate   : 0.0414 cm^3/s
Extract flow rate     : 0.0348 cm^3/s
Raffinate flow rate   : 0.0266 cm^3/s

The developed code replicated both cases successfully. In regard to the simulation time (under the tolerance 1e-4), 

	- four_column case:
		Time: 96 sec
		Round: 56
		Switch: 224

	- eight_column case
		Time: 105 sec
		Round: 34
		Switch: 272

As seen from above data, the simulation time does't increase exponentially according to the number of the column. And we also acquired the same cyclic steady state with regular SMB. This is the advantage.

------------------------------------------------------------------------------------------------
The third version of the eight_column SMB.

All the bugs mentioned aboved were all fixed in the third version. The reproduced data were from the paper Model-based Control of a Simulated Moving Bed Chromatographic Process for the Separation of Fructose and Glucose.

k_a & k_d : [0.54 0.28]  (Linear isotherm);
d_{col}   : 2.6cm;
L         : 53.6cm;
porosity  : 0.38;
d_{par}   : 0.325cm;
Switching time        : 1552 s
concentration of feed : 0.5 g/cm^3
Recycle flow rate     : 0.1395 cm^3/s
Feed flow rate        : 0.02 cm^3/s
Desorbent flow rate   : 0.0414 cm^3/s
Extract flow rate     : 0.0348 cm^3/s
Raffinate flow rate   : 0.0266 cm^3/s

The convergence to the Cyclic Steady State was achieved after 104 switching. This resulted in approx. 296 second on my desktop under the relative tolerence 1e-4.
--------------------------------------------------------------------------------------------------
I extended the four-column Simulated Moving Bed (SMB) chromatography to eight-column situation.

Two bugs remained in the version of four-column were fixed. 1) the constraint of Length/interstitial_velocity is larger than switch time, as the other papers adopted. 2) All the parameter data are from the paper "Model-based control of a simulated moving bed chromatographic process for the separation of fructose and glucose", rather than artificial as in the four-column case. 


The bug left now is the initial setup of feed concentration and the injection during the process. Also the code is not compatible with the four-column situation. This will be improved in the later version.
  
