#!/usr/bin/env python
# coding: utf-8

### Testing Model Predictive Control on Trait-Based gLV systems
# willsharpless@berkeley.edu or ucsd.edu

# In[1]:

from pandas.core.arrays.sparse import dtype
import numpy as np
import pandas as pd
import sys, os
from casadi import *
import matplotlib.pyplot as plt

# Add do_mpc to path. This is not necessary if it was installed via pip
sys.path.append('../../../')

# Import do_mpc package:
import do_mpc
import matplotlib.pyplot as plt
import itertools
import time

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def limiting_similarity(T, E):
    A = np.zeros([len(T), len(T)])
    for i in range(0, len(T)):
        for j in range(0, len(T)):
            if T[i] != T[j]: 
                A[i,j] = (-1 / (1 + E)) * abs(1 / (T[i] - T[j]))
            else:
                A[i,j] = -1
    return(A)

def competitive_hierarchy(T, E):
    A = np.zeros([len(T), len(T)])
    for i in range(0, len(T)):
        for j in range(0, len(T)):
            if i != j: 
                comp = 0
                if T[i] > T[j]:
                    comp = -(T[i] - T[j])
                A[i,j] = (1 / (1 + E)) * comp
            else:
                A[i,j] = -1
    return(A)

# def shannon(sb_N, n_sp):
#     sum_N = np.sum(sb_N,axis=1)
#     pi = np.divide(sb_N,np.repeat(sum_N.reshape(len(sum_N),1), n_sp, axis=1))
#     sh_N = -1*np.log(np.prod(numpy.nan_to_num(np.power(pi, pi), nan=1.0),axis=1))
#     return sum_N, sh_N

# %% Define Constants

save_path = 'NEED TO BE DEFINED'

# General
n_rand = 2; #10
E = 8; 
delta_t = 1e-1; 
n_step = 100; 
n_horizon = 3; 
model_type = 'continuous'
n_sp = 3; #10
target_cwms = [0.2, 0.4, 0.6, 0.8]

# Model
r = np.repeat(1, n_sp).reshape(n_sp,1)
B = np.repeat(0, n_sp).reshape(n_sp,1)
B[n_sp-1], B[0] = 1, 1 # species w max & min traits are controllable

# Scoring
success_radius = 0.1
columns = ["Final_Trait_Composition",
              "Target_Trait_Composition",
              "Numerical_Convergence",
              "Initial_State",
              "Final_State",
              "Initial_Richness",
            ]
alpha = 0.05 #error
z = 1 - alpha/2

# %% Run MPC Sims

for hi, hypothesis in enumerate([limiting_similarity, competitive_hierarchy]):

    # MPC scoreboard data arrays
    n_runs = len(target_cwms) * n_rand * (2**n_sp)
    sb_t_cwm = 100*np.ones((n_runs, 1)) # Final T_cwm
    sb_target_t_cwm = 100*np.ones((n_runs, 1)) # Target T_cwm
    sb_conv = np.zeros((n_runs, 1), dtype=bool) # Numerical convergence?
    sb_Ni = np.zeros((n_runs, n_sp)) + 1e-4 # Initial _N (w a numerical error patch)
    sb_Nf = 100*np.ones((n_runs, n_sp)) # Final _N
    sb_n_sub_ic = np.zeros((n_runs, 1)) # Number of nontrivial species in initial subcommunity == richness

    # auxiliary
    skip = np.zeros((n_runs), dtype=bool) # Skip N_j < 0 equil
    c_run = 0

    for rand in range(n_rand):

        start = time.time()
        print("Beginning Randomization ", rand, " for hypothesis ", hypothesis.__name__, "...")
        np.random.seed(10 + rand) # same background

        for target_cwm in target_cwms:

            print(" Beginning MPC runs with "+str(n_sp)+" sp and "+str(target_cwm)+" target...")

            #### Model & MPC Definition (sadly, do_mpc does not support parameter redefinition)
            model = do_mpc.model.Model(model_type)

            # States
            _N = model.set_variable(var_type='_x', var_name='_N', shape=(n_sp,1))
            _u = model.set_variable(var_type='_u', var_name='_u', shape=(n_sp,1))

            # draw trait values from uniform random distribution, sort in order
            T = np.sort(np.random.uniform(0,1,n_sp)).reshape(n_sp,1)
            A = hypothesis(T,E)

            # Make i/c at all equilibria of defined glv system
            c_run += 1 #trivial fp is trivial
            for eq_sz in range(1, n_sp+1):
                for sub_idxs in itertools.combinations(np.arange(n_sp), eq_sz):
                    sub_idxs = np.array(sub_idxs)
                    eq_sub = -1*np.linalg.inv(A[:,sub_idxs][sub_idxs, :])@r[sub_idxs]
                    sb_Ni[c_run, sub_idxs] = eq_sub.T
                    sb_n_sub_ic[c_run] = eq_sz
                    if any(eq_sub < 0): skip[c_run] = 1 #skip subcomms with not all-positive eq
                    c_run+=1
            c_run -= int(2**n_sp) #reset for sims (num of eq = 2^n)

            # Set the differential equation (GLV plus a perturbation to abundance)
            model.set_rhs('_N', (_N*(r+A@_N)+B*_u))

            # Make a CWM target
            T_cwm = model.set_expression(expr_name='T_cwm', expr = sum1(T*_N)/sum1(_N))
            N_tot = model.set_expression(expr_name='N_tot', expr = sum1(_N))
            richness = model.set_expression(expr_name='richness', expr = sum1(_N>0.001))
            _N_dot = model.set_expression(expr_name='_N_dot', expr = (_N*(r+A@_N)+B*_u))

            model.setup()
            mpc = do_mpc.controller.MPC(model)
            setup_mpc = {
                'n_horizon': n_horizon,
                'n_robust': 1,
                'open_loop': 0,
                't_step': delta_t,
                'state_discretization': 'collocation',
                'collocation_type': 'radau',
                'collocation_deg': 2,
                'collocation_ni': 2,
                'store_full_solution': True,
            }
            mpc.set_param(**setup_mpc)

            # Define MPC Objective
            mterm = sum1(model.aux['T_cwm'] - target_cwm)**2 * 1/(model.aux['richness']+0.1)
            lterm = sum1(model.aux['T_cwm'] - target_cwm)**2 * 1/(model.aux['richness']+0.1)

            mpc.set_objective(mterm=mterm, lterm=lterm)
            mpc.set_rterm(_u=1e-3) # input change penalty

            # State bounds
            min_N = np.repeat(1e-8, n_sp)
            max_N = np.repeat(5, n_sp)
            min_u = np.repeat(-10, n_sp)
            max_u = np.repeat(10, n_sp)
            mpc.bounds['lower','_x','_N'] = min_N
            mpc.bounds['upper','_x','_N'] = max_N
            mpc.bounds['lower','_u','_u'] = min_u
            mpc.bounds['upper','_u','_u'] =  max_u
            mpc.setup()

            # Define observer/simulator
            estimator = do_mpc.estimator.StateFeedback(model)
            simulator = do_mpc.simulator.Simulator(model)
            params_simulator = {
                'integration_tool': 'cvodes',
                'abstol': 1e-8,
                'reltol': 1e-8,
                't_step': delta_t
            }
            simulator.set_param(**params_simulator)
            simulator.setup()

            #### Simulation from all x0
            for _ in range(2**n_sp):
                x0 = sb_Ni[c_run,:]
                # print("Driving MPC from x0 =", np.round(x0,3))

                if skip[c_run]: c_run += 1; continue #print(" - Skipped");

                mpc.reset_history()
                x_curr = np.copy(x0)
                mpc.x0 = x_curr; simulator.x0 = x_curr; estimator.x_curr = x_curr
                mpc.set_initial_guess()

                # Simulate
                _N_ddot = np.zeros((n_step,n_sp))
                with HiddenPrints():

                    c_conv = 0 #convergence counter

                    for k in range(n_step):
                        c_conv = max(c_conv-1, 0)
                        u_curr = mpc.make_step(x_curr)
                        y_next = simulator.make_step(u_curr)
                        x_curr = estimator.make_step(y_next)

                        # Compute d^2/dt^2
                        if k > 0:
                            _N_dot_f = mpc.data['_aux'][-1,4:4+n_sp]
                            _N_f = mpc.data['_x'][-1, :]
                            _u_dot = np.diff(mpc.data['_u'],axis=0)/delta_t
                            _u_dot_f = _u_dot[-1,:]
                            _N_ddot_f = _N_dot_f*(r+A@_N_f)+_N_f*(r+A@_N_dot_f)+B*_u_dot_f
                            _N_ddot[k,:] = _N_ddot_f[0,:]
                        
                        # ~numerically converging := bounded L1(N_dot) and L1(N_ddot)
                        if k > 0 and sum(abs(_N_dot_f)) < 0.5 and sum(abs(_N_ddot_f[0,:])) < 7.5:
                            c_conv += 2 

                        # ~convergence := numerically converging for 25 of the last 50 steps
                        if c_conv > 49:
                            stat, statbool = " - Converged!", 1
                            break
                        elif k == n_step - 1:
                            stat, statbool = " - Failed to converge :(", 0
                            break
                               
                # Store data
                sb_t_cwm[c_run] = mpc.data['_aux'][-1,1] # Final T_cwm
                sb_target_t_cwm[c_run] = target_cwm  # Final T_cwm
                sb_conv[c_run] = statbool # Numerical convergence?
                sb_Nf[c_run,:] = x_curr.T # Final _N

                # Plot
                # fig, ax, graphics = do_mpc.graphics.default_plot(mpc.data)
                # graphics.plot_results()
                # graphics.reset_axes()
                # plt.suptitle("x0 = " + str(np.round(x0.T,2)))
                # plt.show()

                c_run += 1
                # break
            # break

        end = time.time()
        print(" completed in ", np.round(end - start,1), "seconds\n")
        # break
    
    # Clean scored data of skipped sims
    sb_t_cwm = sb_t_cwm[~skip]
    sb_target_t_cwm = sb_target_t_cwm[~skip]
    sb_conv = sb_conv[~skip] 
    sb_Ni = sb_Ni[~skip,:]
    sb_Nf = sb_Nf[~skip,:]
    sb_n_sub_ic = sb_n_sub_ic[~skip]
   
    # combine and add to hypothesis dataframe
    data = [sb_t_cwm, sb_target_t_cwm, sb_conv, sb_Ni, sb_Nf, sb_n_sub_ic]
    np.save(save_path + "trait_MPC_n10_rand" + str(rand) + "_data_" + hypothesis.__name__, data)

    ## %% Compute success, combine into dataframe and export to csv

    df = pd.DataFrame()
    for i, col in enumerate(columns):
        if i != 3 and i != 4:
            df[col] = data[i].flatten()
        else:
            df[col] = data[i].tolist()

    # df["Final_Richness"] = np.where(...)

    near_target = (np.abs(np.subtract(data[0], data[1])) <= success_radius)
    converged = (data[2] > 0)
    success_idx = np.where(near_target & converged)[0]
    success = np.zeros(len(data[0]), dtype=bool)
    success[success_idx] = True
    df["Controller_Success"] = success

    rand_sim_lengths = np.diff(np.append(np.where(data[5] == 0)[0][::4], len(data[5])))
    df["Randomization"] = np.concatenate([(i+1)*np.ones(ix, dtype=int) for i, ix in enumerate(rand_sim_lengths)])
    
    print("Writing " + ["limiting_similarity", "competitive_hierarchy"][hi] + " df to file!")
    df.to_pickle(["limiting_similarity", "competitive_hierarchy"][hi] + "_mpc_data_full.pkl")

    df_short = df[["Randomization", "Initial_Richness", "Target_Trait_Composition", "Controller_Success"]].copy()
    
    df_summary = df_short.groupby(['Initial_Richness','Target_Trait_Composition']).mean().reset_index()

    means = df_short.groupby(['Initial_Richness','Target_Trait_Composition']).mean()['Controller_Success'].values
    n_samples = df_short.groupby(['Initial_Richness','Target_Trait_Composition']).size().values

    ci_size = z * np.sqrt(np.divide(np.multiply(means, 1 - means), n_samples))
    df_summary["Confidence_Interval_95"] = ci_size

    df_summary.rename(columns = {"Controller_Success":"Mean_Controller_Success"})
    df_summary.drop(columns = ["Randomization"], inplace=True)
    
    df_summary.to_csv(["limiting_similarity", "competitive_hierarchy"][hi] + "_mpc_stats_10n.csv")
    print("Finished and scored simulations for hypothesis ", hypothesis.__name__, " ========================== \n")
