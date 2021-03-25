# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 12:21:31 2018

@author: sphere
"""


import netCDF4
import numpy as np
import glob
import json
import os

if __name__ == '__main__':
    
    # from process_args import TincArgumentParser
    # parser = TincArgumentParser(description='Analyze parameter space')
    # parser.add_argument('__output_names', type=str, default="parameter_space.nc", nargs='?')
    # args = parser.get_args()

    sub_dirs = glob.glob('*/')
    for sub_dir in sub_dirs:
        # Quick and dirty way to check if this is a CASM run. Perhaps should be improved?
        if os.path.exists(sub_dir + 'conditions.0/' and os.path.exists(sub_dir + "results.json")):

            files = glob.glob(sub_dir + "*.nc")
            shell_sites_files = []
            perco_files = []

            for filename in files:
                ncfile = netCDF4.Dataset(filename, mode='r', format='NETCDF4') 
                varnames = ncfile.variables.keys()
                if "shell_sites" in varnames and "occ_ref" in varnames:
                    shell_sites_files.append(filename[len(sub_dir):])
                elif "occupation_dof" in varnames:
                    perco_files.append(filename[len(sub_dir):])
                else:
                    print("Unrecognized file:" + filename)

            if len(shell_sites_files) > 0:
                names = {"files": shell_sites_files}
                print(f"Writing shell_site_files.nc in {sub_dir}")
                with open(sub_dir + "shell_site_files.json", 'w') as j:
                    json.dump( names, j)
            
            
            if len(perco_files) > 0:
                names = {"files": perco_files}
                print(f"Writing perco_files.nc in {sub_dir}")
                with open(sub_dir + "perco_files.json", 'w') as j:
                    json.dump( names, j)

        
    
    

