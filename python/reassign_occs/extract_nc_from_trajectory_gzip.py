# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:23:17 2019

@author: andres
"""
import os
import numpy as np
import json
import glob
import gzip
import netCDF4

def convert_gzip(file_path, result_path, verbose):
    if not os.path.exists(result_path):
        if os.path.exists(file_path):

            with gzip.GzipFile(file_path, 'rb') as f:
                results_json = json.loads(f.read().decode('UTF-8'))
                ncfile = netCDF4.Dataset(result_path, mode='w', format='NETCDF4') 
                occupation_dim = ncfile.createDimension('index', size=len(results_json['DoF'][0]['occupation']))

                min_len = min(len(results_json['DoF']), len(results_json['Step']))
                step_dim = ncfile.createDimension('step', size=min_len)

                steps = ncfile.createVariable('steps', np.uint32,('step'))
                occupation = ncfile.createVariable('occupation_dofs', np.uint8,('step', 'index'))
                for i, time_step in enumerate(results_json['DoF'][:min_len]):
                    occupation[i] = np.array(time_step['occupation'], dtype = np.uint8)
                    # print(occupation[i][:])
                steps[:] = results_json['Step'][:min_len]

                ncfile.close()
                if verbose:
                    print("Wrote: " + result_path)
        else:
            if verbose:
                print("Trajectory gzip file not found")
    else:
#         Output exists, no need to recreate
        pass

def process_dir(dir, in_name, out_name, verbose):
    dirs_with_traj = [x[0] for x in os.walk(dir) if in_name in x[2]]
    for d in dirs_with_traj:
        if verbose:
            print("Converting in " + d)
        convert_gzip(d + "/" + in_name, d + "/" + out_name, verbose)

if __name__ == '__main__':

    import sys, os
    sys.path.append(
        os.path.abspath(os.path.join(os.path.dirname(__file__), "../../external/tinc/tinc-python/tinc-python")))
    from process_args import *

    parser = TincArgumentParser(description='Kinetic MonteCarlo data converter')
    args = parser.get_args()

    root_dir = args['__input_dir']
    process_dir(root_dir, args['__input_names'][0], args['__output_names'][0], args["__verbose"])
    #result_path = args['__output_dir'] + args['__output_name']


