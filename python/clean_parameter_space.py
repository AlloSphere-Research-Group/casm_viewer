# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 09:19:21 2019

@author: andres
"""

import sys
import os
import glob

if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description='Analyze parameter space')
    
    parser.add_argument('directory', type=str, default="/2809share/vdv_data_temp/monte_with_traj_quick/__2809share_vdv_data_temp_monte_with_traj_quick__config.json", nargs='?')
    parser.add_argument('output_name', type=str, default="_parameter_space.json", nargs='?')
    
    args = parser.parse_args()
    root_path = args.directory
    
    output_name = args.output_name
    print("Analyzing " + root_path, file = sys.stderr)
    
    if not root_path[-1] == '\\' and not root_path[-1] == '/':
        root_path += '/'
        
    if os.path.exists(root_path + output_name):
        os.remove(root_path + output_name)
    if os.path.exists(root_path + output_name + ".meta"):
        os.remove(root_path + output_name + ".meta")
    sub_dirs = glob.glob(root_path + '*/')
    for sub_dir in sub_dirs:
        if os.path.exists(sub_dir + output_name):
            os.remove(sub_dir + output_name)
        if os.path.exists(sub_dir + output_name + ".meta"):
            os.remove(sub_dir + output_name + ".meta")
            