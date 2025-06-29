## Empirical Acoustic Diffusion Equation Coefficient
In this repository, the optimization and diffusion equation model scripts for optimizing a spatially varying diffusion coefficient are at hand. This repository is the complementary material connected to the submitted paper _Fichera, I., Van hoorickx, C., Hornikx, M., An empirical diffusion coefficient function for the acoustic diffusion equation model in long rooms, Submitted to Applied Acoustics (2025)_.

The software implementation of the Acoustic Diffusion Equation is part of an ongoing research in the Building Acoustics Group at the Built Environment Department of Eindhoven University of Technology.
The main software can be found in the [GitHub Repository](https://github.com/Building-acoustics-TU-Eindhoven/Diffusion) and is currectly **UNDER DEVELOPMENT**. It is being implemented by Ilaria Fichera in Python programming language. The diffusion equation implementation of this repository is based on the numerical Finite Different Method (FDM) by Du Fort&Frankel (Navarro et al., 2012) from the main project [GitHub Repository](https://github.com/Building-acoustics-TU-Eindhoven/Diffusion).

## Setup and Usage instructions
In order to run the scripts in this repository, the following first steps will need to be done.
1. Download and install Anaconda or download and install any Python software IDE you prefer;
2. Clone/Fork this repository to a folder of your preference;
3. You will be able to open all the files through the preferred IDE and test the software.

To properly run the software in Python, the following libraries are needed:
- Python version 3.10.9 or above
- math
- matplotlib
- numpy
- scipy
- sys
- time
- os

### Repository structure
The repository consists on the following structure:
+ `Results` folder: 
This folder contains the results file for each simulation for each room analysed. 
Within this folder, there are folders with the following name code, e.g. `Results/NPR39x3x3-S1.5,1.5,1.5-Rx,1.5,1.5-alpha0.1` folder, where _NRP_ stands for non-proportional room, _39x3x3_ is the dimension of the room, _S_ stands for source, _1.5,1.5,1.5_ are the three coordinates of the source position, _R_ stands for recevier, _x,1.5,1.5_ are the coordinates of the receiver position and _alpha_ is the uniform absorption coefficient of the surfaces of the room _0.1_. It is included an _x_ in the receiver position because the objective is to find the Sound Pressure Level at each x position, keeping constant the y and z position. The coordinate _x_ of the receiver can be chosen as preferred, however the _x_ coordinate of the source and the _x_ coordinate of the receiver cannot be the same. Within each room folder (e.g. `Results/NPR39x3x3-S1.5,1.5,1.5-Rx,1.5,1.5-alpha0.1` folder), there is one simulation results folder, `results_rad_imp`, which cointain the results from the simulation with the radiosity method. This is the reference method. These have been simulated using the Radiosity software provided by the Technical University of Denmark (DTU) and Dr. George Koutsouris (this software is not included in this repository). 
The results from the simulation with the diffusion equation with constant diffusion coefficient are saved in the folder `results_diff_imp`, while the results from the simulation with the diffusion equation with spatially dependent diffusion coefficient are saved in the folder `results_diff_opt`.
+ `results_diff_emp` local folder to store results of the simulation done on the moment by the diffusion equation code with empirical diffusion coefficient
+ `results_diff_imp` local folder to store results of the simulation done on the moment by the diffusion equation code with constant diffusion coefficient
+ `results_diff_opt` local folder to store results of the simulation done on the moment by the diffusion equation with spatially dependent diffusion coefficient
+ `results_rad_imp` local folder to store results of the simulation done on the moment by the radiosity method
+ Python scripts

### Workflow for correctly running the software

1. Decide which room to simulate. For example let's consider for this tutorial the room 39x3x3 $m^3$ with Source S1.5m,1.5m,1.5m, Receiver R8m,1.5m,1.5m and alpha 0.1.
2. If you want to simulate the room 39x3x3 $m^3$ with Source S1.5m,1.5m,1.5m, Receiver R8m,1.5m,1.5m and alpha 0.1, you need to make sure that the local folder `results_rad_imp` has the Radiosity method results file of the room you are simulating. The Radiosity method results files are pre-prepared as the reference results and they are stored in each room result folder.
3. Go to the `Results/NPR39x3x3-S1.5,1.5,1.5-Rx,1.5,1.5-alpha0.1` and copy the `results_rad_imp` folder and paste it in the `C:\....\Empirical_Acoustic_Diffusion_Equation_Coefficient` folder.
4. Open the first python script `1-DiffEq3D.py`.
5. In the "INPUT VARIABLE" section of the script, enter the length (x-axis), width (y-axis), heigth (z-axis) of the room you want to simulate, the source position coordinates (x_source, y_source, z_source) and the receiver position coordinates (y_rec, z_rec) in meters and the absorption coefficient (alpha_1,alpha_2 etc...). For these specific simulations, the x-axis coordinate of the receiver is not important since the objective is to calculate all the acoustic variables on a line on the x-axis. the _x_ coordinate of the source and the _x_ coordinate of the receiver cannot be the same.
6. Press run. While running you will see printed the percentage of completion of calculation. Once arrived at 100%, some Figures will show and all the relevant results will be saved in the local folder `results_diff_imp`.
7. Open the second python script `2-DiffEq3D-WithLoop.py`. If you have set up correctly the first script, you need to only press run. This script will calculate the reverberation time for each position x of the room and it will save the results in the local folder `results_diff_imp`.
8. If you have done all the steps above, open the third python script `3-OptimizationDx.py` and press run. The results will be saved in the local folder `results_diff_opt`.
The time length of the optimization simulation is difficult to quantify as it depends on how much different are the SPL curve over the distance between the radiosity method and the diffusion equation method. It will run until 100% as many times as the quantity C2 will need to be optimised to minimise the difference between the two curves. The calculation could be in between 15 to 30 hours, depending on the room that you are simulating and on the optimization variables.
9. Open the fourth python script `4-DiffEq3DWithNewDx.py` and press run. The results will be saved in the local folder `results_diff_opt`.
10. Open the fourth python script `5-DiffEq3DWithEmpDx.py` and press run. The results will be saved in the local folder `results_diff_emp`.
11. Once you have all the results, you can create figures and compare the different curves.

If some error occurs it might mean that the results on the local folders `results_diff_imp`, `results_diff_opt` and `results_rad_imp` are of different rooms and therefore there are problems with length of vectors and arrays.

## Explanation of what each script calculates

+ File `1-DiffEq3D.py` calculates the Sound Pressure Level (SPL) over the space and over time with a constant diffusion coefficient, but also other room acoustics parameters such as early decay curve, reverberation time at the receiver position, clarity etc... It creates also some figures of interesting results for the room simulated.

+ File `2-DiffEq3D-WithLoop.py` calculates specifically the Reverberation Time (RT) over the distance x with the constant diffusion coefficient. 

+ File `3-OptimizationDx.py` optimize the x-axis spatially dependent diffusion coefficient to match the SPL results over distance to the reference radiosity method results using the least-square optimization method. In particular, the script will use some initial values of the constants C2  and optimises these constants until the difference between the diffusion equation SPL over the x-axis and the radiosity SPL over the x-axis is minimised.  

+ File `4-DiffEq3DWithNewDx.py` can calculate all the acoustics parameters of the room in question with the spatially dependend diffusion coefficient calculated by the optimization. The script incorporates the fact that the diffusion coefficient is now a vector compared to previously when it was only a constant variable.

+ File `5-DiffEq3DWithEmpDx.py` can calculate all the acoustics parameters of the room in question with the spatially dependend diffusion coefficient calculated by the empirical model. The script incorporates the fact that the diffusion coefficient is calculated not by optimization but with the empirical formula.

For all the scripts, there are some associated functions (these should not be modified):
+ `FunctionRT.py` calculates the reverberation time of the room in question
+ `FunctionEDT.py` calculates the early decay time of the room in question
+ `FunctionClarity.py` calculates the clarity $C_{80}$ of the room in question based on Barron's revised theory formula.
+ `FunctionDefinition.py` calculates the definition $D_{50}$ of the room in question based on Barron's revised theory formula.
+ `FunctionCentreTime.py` calculates the centre time $T_s$ of the room in question based on Barron's revised theory formula.
+ `DiffEq3DFunction.py` is the function used for the optimization.

## Authors
Software is being developed by Ilaria Fichera at Eindhoven University of Technology (TU/e).

## License
Diffusion Equation Code and Optimization Code are licensed under GNU General Public License v2.0. See LICENSE.md for more details.
