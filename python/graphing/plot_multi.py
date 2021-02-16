#!/bin/env/python

# Jonas Kaufman and Andres Cabrera
# August 2018

from __future__ import print_function
import numpy as np
import os,sys,json
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Argument parsing
import argparse

parser = argparse.ArgumentParser(description='Graph generator (multi)')
parser.add_argument('casm_path', type=str, help='input filename', default = '/Users/Jonas/Research/NaxCoO2/P3', nargs='?')
parser.add_argument('heating_path', type=str, help='input filename', default = 'monte/heating_01', nargs='?')
parser.add_argument('cooling_path', type=str, help='input filename', default = 'monte/cooling_01', nargs='?')
parser.add_argument('path_template', type=str, help='sub path template', default = "mu_{}", nargs='?')
parser.add_argument('temp_interest', type=str, help='temperature of interest', default = "", nargs='?')
parser.add_argument('title', type=str, help='title of the dataset', default = "", nargs='?')
parser.add_argument('output_dir', type=str, help='Output directory', default = 'cached_output', nargs='?')
parser.add_argument('output', type=str, help='output filename', default = 'mc_P3_01.png', nargs='?')

args = parser.parse_args()


casm_path = args.casm_path
# constT_path = 'monte/constT_01'
heating_path = args.heating_path
cooling_path = args.cooling_path
fig_save_path = args.output_dir + args.output
path_template = args.path_template
temp_interest = args.temp_interest

def getSlopes(XX,YY):
    slopes = np.zeros((len(XX),2))
    slopes[:,0] = XX
    slopes[:-1,1] = [(YY[i+1] - YY[i])/(XX[i+1] - XX[i]) for i in range(len(XX)-1)]
    slopes[-1,1] = slopes[-2,1]
    return slopes

# Limits for heating and cooling runs
T_limits = [5,1000]
# mu_limits = [-3.5,5.5] # param_chem_pot not mu itself
mu_limits = [-1.0, 1.0] # param_chem_pot not mu itself
diff_T = 5
diff_mu = 0.05

# num_T = int(np.rint((T_limits[1]-T_limits[0])/diff_T)+1)
# T_vals = np.linspace(T_limits[0],T_limits[1],num_T).astype(int) 
num_mu = int(np.rint((mu_limits[1]-mu_limits[0])/diff_mu)+1)
mu_vals = [round(i,3) for i in np.linspace(mu_limits[0],mu_limits[1],num_mu)]
mu_factor = 2.0 # factor between parametric mu and actual
mu_norm = [m/mu_factor for m in mu_vals]

save_fig = False
if len(temp_interest) >= 0:
    T_interest = int(float(temp_interest))	# temperature of interest for voltage curve
#    print(T_interest, sys.stderr)
    save_fig = True

x_heat = None
x_cool = None
# Extract results
for i,mu in enumerate(mu_vals):
    # Heating
    curr_path = os.path.join(casm_path,heating_path,path_template.format(np.round(mu,3)))
    heat_data = json.loads(open(os.path.join(curr_path,'results.json')).read())
    xh = np.array(heat_data["<comp_n(Na)>"])/np.array(heat_data["<comp_n(Co)>"])
    # Cooling
    curr_path = os.path.join(casm_path,cooling_path,path_template.format(np.round(mu,3)))
    cool_data = json.loads(open(os.path.join(curr_path,'results.json')).read())
    xc = np.array(cool_data["<comp_n(Na)>"])/np.array(cool_data["<comp_n(Co)>"])

    if x_heat is None:
        # Assumes "T" field is identical on all results. reasonable?
        T_vals = cool_data["T"]
        num_T = len(T_vals)

        print(T_vals, file=sys.stderr)
        x_heat = np.zeros(shape=(num_mu,num_T))
        x_cool = np.zeros(shape=(num_mu,num_T))

    for j in range(num_T):
        x_heat[i,j] = xh[j]
        x_cool[i,j] = xc[num_T-j-1]
        
# Plot
fig, ax = plt.subplots(3, sharex=True)
fig.subplots_adjust(hspace=0)
ax[0].set_ylabel('T (K) heat')
ax[1].set_ylabel('T (K) cool')
ax[2].set_ylabel('mu (eV/Na)')
plt.xlabel('$x$ in ' + args.title)
colors = cm.plasma(np.linspace(0,1,num_mu))
t = T_vals.index(T_interest)
# T vs x 
sz = 5
#print 'mu\tx_heat_{0}\tx_cool_{0}'.format(T_vals[t])
for i,mu in enumerate(mu_vals):
    ax[0].scatter(x_cool[i,:],T_vals,sz,edgecolors='none',color=colors[i])
    ax[1].scatter(x_heat[i,:],T_vals,sz,edgecolors='none',color=colors[i])
    #print '{0}\t{1}\t{2}'.format(mu_norm[i],x_heat[i,t],x_cool[i,t])
# mu vs x
ax[2].plot(x_cool[:,t],mu_norm,'-',label='{} K cool'.format(T_vals[t]))
ax[2].plot(x_heat[:,t],mu_norm,'-',label='{} K heat'.format(T_vals[t]))
ax[2].legend(loc=4)
if save_fig:
    plt.savefig(fig_save_path,dpi=300)
else:
    plt.show()
