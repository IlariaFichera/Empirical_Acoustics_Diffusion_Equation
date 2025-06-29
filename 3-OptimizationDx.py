# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 09:49:08 2023

@author: 20225533
"""
import scipy
from scipy.optimize import least_squares
import numpy as np
from os.path import dirname, join as pjoin
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.linalg import norm
import os

def calculate_spl_rt(k):
    from DiffEq3DFunction import calculate_spl_rt_diffusion
    t30_x, Dx, spl_rec_x_t0 = calculate_spl_rt_diffusion(k,length,width,height,x_source,y_source,z_source,x_rec,y_rec,z_rec,alpha)
    return t30_x, Dx, spl_rec_x_t0

def obj_fun(k):
    print(k)
    D_list.append(k)
    t30_x, Dx, spl_rec_x_t0 = calculate_spl_rt(k)
    
    SPL_t0_R = spl_contents["SPL_t0_R"][0]
    t30_x_R = rt_contents["T30_x"][0]
    
    mean_free_path_rounded = round(mean_free_path*2)/2
    
    if x_source <= mean_free_path_rounded:
        pos_from = x_source + mean_free_path_rounded
        temp = abs(pos_from - x_axis)
        x_from = x_axis[np.argmin(temp)]
        idx_dist = np.where(x_axis == x_from)[0][0]
        SPL_t0_R_nosource = SPL_t0_R[idx_dist:-1] #radiosity
        spl_rec_x_t0_nosource = spl_rec_x_t0[idx_dist:-1] #diffusion
        t30_x_R_nosource = t30_x_R[idx_dist:-1] #radiosity
        t30_x_nosource = t30_x[idx_dist:-1] #diffusion
        
        spl_diff = np.sqrt(A)*((spl_rec_x_t0 - SPL_t0_R)/1)
        #spl_diff[0:idx_dist3] = 0
        rt_diff = np.sqrt(B)*((t30_x - t30_x_R)/(0.05*t30_x_R))
        rt_diff[0:idx_dist] = 0
        tot_length = len(spl_diff) + len(rt_diff)
        tot_diff = np.zeros((tot_length)) 
        tot_diff[0:len(spl_diff)] = spl_diff
        tot_diff[len(spl_diff):tot_length] = rt_diff
        
        
    else:
        pos_source = np.where(x_axis == x_source)
        pos_source_before = np.where(x_axis == x_source-mean_free_path_rounded)[0][0]
        pos_source_after = np.where(x_axis == x_source+mean_free_path_rounded)[0][0]
        SPL_t0_R[pos_source_before+1:pos_source_after-1] = 0 #radiosity
        
        SPL_t0_PY_after_nosource = spl_rec_x_t0 #diffusion
        SPL_t0_PY_after_nosource[pos_source_before+1:pos_source_after-1] = 0 #diffusion
        
        t30_x_R[pos_source_before+1:pos_source_after-1] = 0 #radiosity
        t30_x_after_nosource = t30_x #diffusion
        t30_x_after_nosource[pos_source_before+1:pos_source_after-1] = 0 #diffusion
        
        spl_diff = np.sqrt(A)*((spl_rec_x_t0 - SPL_t0_R)/1)
        
        rt_diff = np.sqrt(B)*((t30_x - t30_x_R)/(0.05*t30_x_R))
        rt_diff[pos_source_before+1:pos_source_after-1] = 0
        tot_length = len(spl_diff) + len(rt_diff)
        tot_diff = np.zeros((tot_length)) 
        tot_diff[0:len(spl_diff)] = spl_diff
        tot_diff[len(spl_diff):tot_length] = rt_diff
    

    
    cost_spl = 0.5*(sum(spl_diff**2))
    cost_rt = 0.5*(sum(rt_diff**2))
    cost = 0.5*(sum(spl_diff**2))+0.5*(sum(rt_diff**2))
    cost_list.append(cost)
    cost_spl_list.append(cost_spl)
    cost_rt_list.append(cost_rt)
    
    return tot_diff

#INPUT VARIABLES
current_path = os.getcwd()
results_diff_imp = os.path.join(current_path, 'results_diff_imp')
#results_rad_imp = os.path.join(current_path, 'results_rad_imp')

D_th = np.load(os.path.join(results_diff_imp, 'D_th.npy'))  
RT_Sabine = np.load(os.path.join(results_diff_imp, 'RT_Sabine.npy')) 
c0 = np.load(os.path.join(results_diff_imp, 'c0.npy')) 
alpha = np.load(os.path.join(results_diff_imp, 'alpha.npy')) 
length = np.load(os.path.join(results_diff_imp, 'length.npy')) 
width = np.load(os.path.join(results_diff_imp, 'width.npy'))  
height = np.load(os.path.join(results_diff_imp, 'height.npy')) 
x_source = np.load(os.path.join(results_diff_imp, 'x_source.npy')) 
y_source = np.load(os.path.join(results_diff_imp, 'y_source.npy')) 
z_source = np.load(os.path.join(results_diff_imp, 'z_source.npy')) 
x_rec = np.load(os.path.join(results_diff_imp, 'x_rec.npy')) 
y_rec = np.load(os.path.join(results_diff_imp, 'y_rec.npy'))
z_rec = np.load(os.path.join(results_diff_imp, 'z_rec.npy')) 
mean_free_path = np.load(os.path.join(results_diff_imp, 'mean_free_path.npy'))  
x_axis = np.load(os.path.join(results_diff_imp, 'x_axis.npy'))


C2 = 0
C1 = 0
C0 = D_th
A = 1 #weight for the spl part of the cost
B = 1 #weight for the rt part of the cost

k = [C2,C1,C0] #C2*r**2 + C1*r + C0

D_list = []
cost_list = []
cost_spl_list = []
cost_rt_list = []
data_dir = "results_rad_imp"
spl_fname = pjoin(data_dir, 'SPL_t0_R.mat')
spl_contents = sio.loadmat(spl_fname)
rt_fname = pjoin(data_dir, 'T30_x.mat')
rt_contents = sio.loadmat(rt_fname)

#MAIN OPTIMIZATION CALCULATION
result = least_squares(obj_fun, k, bounds=([0,0,D_th], [np.inf,0.001,D_th+0.001]))

optimal_D = result.x
np.save('results_diff_opt\\optimal_D',optimal_D)

c2_constants = []
for j in D_list:
    c2_constants.append(j[0])

plt.figure(20)
plt.title("Figure 20: D and cost")
plt.plot(c2_constants,cost_list)
plt.plot(c2_constants,cost_spl_list)
plt.plot(c2_constants,cost_rt_list)
plt.legend(['cost total','cost_spl','cost_rt'])
plt.xlabel("c2_constants")
plt.ylabel("cost")

#Calculation of error
n = spl_contents["SPL_t0_R"].shape[1] + rt_contents["T30_x"].shape[1]
RMSD_before = np.sqrt((2/n)*cost_list[0])
RMSD_after = np.sqrt((2/n)*cost_list[-1])

import pickle
import types

# Save all variables to a file
def save(filename):
    with open(filename, 'wb') as f:
        # Filter out modules, functions, and other unsupported types
        filtered_variables = {}
        for k, v in globals().items():
            try:
                # Check if the object can be pickled
                pickle.dumps(v)
                # Exclude some types explicitly known to cause issues
                if not k.startswith('__') and not isinstance(v, (types.ModuleType, types.FunctionType, types.BuiltinFunctionType, types.LambdaType, types.MethodType, types.MappingProxyType)):
                    filtered_variables[k] = v
            except Exception as e:
                print(f"Could not pickle {k}: {str(e)}")

        pickle.dump(filtered_variables, f)

# To save all current variables
save('Optimization.pkl')