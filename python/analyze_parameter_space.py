# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 12:21:31 2018

@author: sphere
"""


import netCDF4
import numpy as np

from multiprocessing.pool import ThreadPool
import glob, re, os, sys

import gzip
import contextlib
import json

import sys
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), "reassign_occs")))
from template_creator import create_template_from_files

def check_equal(iterator: iter):
    iterator = iter(iterator)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == rest for rest in iterator)
    
def reduce_list(input_list):
    if check_equal(input_list):
        return [input_list[0]]
    else:
        return input_list

def process_cur_dir():
    if not os.path.exists('results.json'):
        return {}
    all_params = {}
    with open('results.json') as resuts_file:
        results = json.load(resuts_file)
        all_params['T'] = reduce_list(results['T'])
        all_params['param_chem_pot(a)'] = reduce_list(results['param_chem_pot(a)'])
        if 'param_chem_pot(b)' in results:
            all_params['param_chem_pot(b)'] = reduce_list(results['param_chem_pot(b)'])
        #print('Processing dir ' + directory)
    
    return all_params
    
def analyze_config_files(root_path):
    config_files = {}
    with pushd(root_path):
        files = glob.glob("*.json")
        for filename in files:
            print("Trying file: " + filename)
            try:
                j = json.load(open(filename))
                if 'driver' in j and 'supercell' in j:
                    config_files[filename] = {'supercell': j['supercell'], 'driver' : j['driver']}
            except:
                pass
    return config_files
    
def extract_folder_template(root_path):
    with pushd(root_path):
        dirs = glob.glob("*/")
    #    print(dirs)
        regex = '(-*[0-9]*\.[0-9]+)(?!\.)'
        
        # Identify folder template
        templates = {}
        for dir_name in dirs:
            this_dir_params = []
            matches = re.finditer(regex, dir_name, re.MULTILINE)
            cur_index = 0
            template = ''
            for m in matches:
                template += dir_name[cur_index:m.start()] + '$'
                cur_index += m.end()
                this_dir_params.append(m.group())
            if not template in templates:
                templates[template] = 0
            templates[template] += 1
        max_count = 0
        for template, count in templates.items():
            if count > max_count:
                folder_template = template
                max_count = count
    return folder_template
    
def dir_sort(v):
    try:
        k = int(v)
        return k
    except ValueError:
        return 1000

def dataset_directories():
    dataset_dirs = []
    # Quick and dirty way to check if this is a CASM run. Perhaps should be improved?
    # if os.path.exists('conditions.0/') and os.path.exists("results.json"):
    #     dataset_dirs.append('')
    #sub_dirs = glob.glob('*/')
    sub_dirs =[]
    p=os.listdir('./')
    for i in p:
        if os.path.isdir(i):
            sub_dirs.append(i)
    
    sub_dirs = sorted(sub_dirs, key=dir_sort)
    for sub_dir in sub_dirs:
        # Quick and dirty way to check if this is a CASM run. Perhaps should be improved?
        if os.path.exists(sub_dir + '/conditions.0/') or os.path.exists(sub_dir + "/results.json"):
            dataset_dirs.append(sub_dir + "/")
    return dataset_dirs
    
def get_dataset_params():
    dir_params = []
    all_params = process_cur_dir()
    if len(all_params) > 0:
        new_params = {}
        for param_name, param_values in all_params.items():
            
            if check_equal(param_values):
                if len(param_values) > 1:
                    print("condition is " + param_name)
                value = param_values[0]
            else:
                value = param_values
            new_params[param_name] = value
            
        dir_params.append(new_params)
        
    return dir_params

def write_transfmat(dataset_path, is_collection):
    if is_collection:
        with pushd(dataset_path):
            sub_dirs = glob.glob('*/')
            for sub_dir in sub_dirs:
                config_files = analyze_config_files(sub_dir)
                if len(config_files) > 0:
                    break
    else:
        with pushd(dataset_path):
            config_files = analyze_config_files("./")

    print(f'Writing transfmat from {config_files}')
    if len(config_files) > 0:
        for config in config_files.values():
            # Write supercell into transfmat file
            # This file is used to generate the POSCAR template file
            if 'supercell' in config:
                transfmat = ''
                for row in config['supercell']:
                    transfmat += ' '.join([str(value) for value in row])
                    transfmat += '\n'
                    
                with open(dataset_path + 'cached_output/transfmat', 'w') as f:
                    f.write(transfmat)
                    return True
    return False

def write_dir(output_full_path, dataset_params):
    
    empty_lists = True
    for name,space in dataset_param_meta['internal_params'].items():
        try:
            if len(space) > 0:
                empty_lists = False
        except:
            # Single value
            pass

    for name,space in dataset_param_meta['index_params'].items():
        try:
            if len(space) > 0:
                empty_lists = False
        except:
            # Single value
            pass

    for name,space in dataset_param_meta['mapped_params'].items():
        try:
            if len(space) > 0:
                empty_lists = False
        except:
            # Single value
            pass

    if empty_lists:
        print(f"All lists are empty. Not writing paramter space for {output_full_path}")
        return
    print("Writing output: " + output_full_path)
    ncfile = netCDF4.Dataset(output_full_path, mode='w', format='NETCDF4') 
    # occupation_dim = ncfile.createDimension('index', size=len(results_json['DoF'][0]['occupation']))
    # internal states are values that change inside a specific sample.
    # In our case, this parameter will be read later whn loading the trajecotry file.

    internal_states_group = ncfile.createGroup('internal_dimensions')
    for name,space in dataset_param_meta['internal_params'].items():
        datatype = np.int32
        new_internal_group = internal_states_group.createGroup(name)
        l = 1
        if type(space) == list:
            l = len(space)
        step_dim = new_internal_group.createDimension('dim', size=len(space))
        steps = new_internal_group.createVariable('values', datatype,('dim'), zlib=True)
        steps[:] = space

    # Conditions reflect a parameter that has produced an individual folder with results
    # The naming of these folders follow an unmutable rule. "conditions.0"
    conditions_group = ncfile.createGroup('index_dimensions')
    for name,space in dataset_param_meta['index_params'].items():
        datatype = np.float32
        new_condition_group = conditions_group.createGroup(name)

        l = 1
        if type(space) == list:
            l = len(space)
        condition_dim = new_condition_group.createDimension('dim', size=l)
        condition_var = new_condition_group.createVariable('values', datatype,('dim'), zlib=True)
        condition_var[:] = space

    # Parameters have generated root folders the naming of the folders is arbitrary.
    # Each value and directory represent a full run of the data production process
    parameters_group = ncfile.createGroup('mapped_dimensions')
    for name,space in dataset_param_meta['mapped_params'].items():
        datatype = np.float32
        l = 1
        if type(space["value"]) == list:
            l = len(space["value"])
        new_parameter_group = parameters_group.createGroup(name)
        parameter_dim = new_parameter_group.createDimension(name + '_dim', size=l)
        parameter_var = new_parameter_group.createVariable('values', datatype,(name+ '_dim'), zlib=True)
        parameter_dir_var = new_parameter_group.createVariable('ids', str ,(name+ '_dim'), zlib=True)
        
        parameter_var[:] = space["value"]
        for i,d in enumerate(space["dir"]):
            parameter_dir_var[i] = d
        # for i, sample in enumerate(space.items()):
            # str_out = netCDF4.stringtochar(sample['dir'], 'S4')
            # parameter_var[i] = sample['value']
            # parameter_dir_var[i] = sample['dir']

    ncfile.close()


@contextlib.contextmanager
def pushd(new_dir):
    previous_dir = os.getcwd()
    os.chdir(new_dir)
    yield
    os.chdir(previous_dir)
    
if __name__ == '__main__':

    import sys, os
    from tinc import *
    parser = TincArgumentParser(description='Analyze parameter space')
    parser.add_argument('__output_names', type=str, default="parameter_space.nc", nargs='?')
    args = parser.get_args()

    output_name = args["__output_names"][0]

    dataset_dirs = dataset_directories()
      
    prim_path = ''
    results = []
    
    custom_conditions = []
    
    for file in glob.glob("*.json"):
        try:
            with open(file) as f:
                    j = json.load(f)
                    if "driver" in j:
                        if j["driver"]["mode"] == "custom":
                            custom_conditions = j["driver"]["custom_conditions"][0].keys()
                            print(f"found config file for KMC. Custom conditions: {custom_conditions}")
        except:
            print("failed to load: " + file +". Invalid json? Symlink?")


    dataset_params = {}  

    dataset_param_common = {}  
    dataset_param_meta = {}  
    for dataset_path in dataset_dirs:
        with pushd(dataset_path): 
            results = get_dataset_params()

        if len(results) > 0:    
            dataset_params[dataset_path] = results
        else:
            print("No results for path ", dataset_path)

        # Gather results together
        for dataset_param_meta_new in results:
            if len(dataset_param_common) == 0 and len(dataset_param_meta_new) > 0:
                dataset_param_meta = dataset_param_meta_new
            else:
                for parameter_name in dataset_param_meta_new['parameters'].keys():
                    dataset_param_meta['parameters'][parameter_name] += dataset_param_meta_new['parameters'][parameter_name]

        
    print("Generating shell file and perco lists.")

    if len(dataset_params) == 0:
        print("Directory not a casm dataset")    
        sys.exit(255)


    param_cache = {}
    param_consistent = []
    for path, params in dataset_params.items():
        for param in params:
            for param_name, param_data in param.items():
                if not param_name in param_cache:
                    param_cache[param_name] = param_data
                    param_consistent.append(param_name)
                else:
                    if not param_cache[param_name] == param_data and param_consistent.count(param_name) > 0:
                        param_consistent.remove(param_name)
    
    consistent_params_data = {}
    for param_name in param_consistent:
        consistent_params_data[param_name] = param_cache[param_name]
        
    # Remove consistent parameters from internal values
    for path, params in dataset_params.items():
        for param in params:
            for consistent_param in consistent_params_data.keys():
                del param[consistent_param]
    
    # Now check if KMC dataset and write time dimension

    for path, params in dataset_params.items():
        sub_dirs = glob.glob(path+'/*/')
        for sub_dir in sub_dirs:
            # TODO copy params
            dataset_param_meta = {"internal_params": {}, "index_params": {}, "mapped_params": {}}
            if os.path.exists(sub_dir + "/trajectory.nc"):
                dataset_param_meta["internal_params"]["time"] = []
                break
            
        if len(params) > 0:
            write_dir(path + output_name, dataset_param_meta)

    # dataset_param_meta_new = {}
    # if len(dir_params) > 0:
    #     dataset_param_meta_new = {'conditions': consistent_params_data,
    #                         'internal_states': internal_run_conditions,
    #                         'conditions_map': conditions_map, "parameters" : {}} 
    #     for dir_name, params in dir_params:
    #         for param_name, param_value in params.items():
    #             if not param_name in dataset_param_meta_new['parameters']:    
    #                 dataset_param_meta_new['parameters'][param_name] = []
    #             dataset_param_meta_new['parameters'][param_name].append({'value': param_value, 'dir': dir_name})
    dataset_param_meta = {}

    dir_params = {}
    for dir_name, contents in dataset_params.items():
        for entry in contents:
            for name, data_val in entry.items():
                if not name in dir_params:
                    dir_params[name] = []
                dir_params[name].append((dir_name, data_val))
                if dataset_dirs.count(dir_name) > 0:
                    dataset_dirs.remove(dir_name)

    from operator import itemgetter
    for param_name, param_map in dir_params.items():
        param_map.sort(key=itemgetter(1))

    # if consistent parameter is size 1, we can assume that we are dealing with a single directory, 
    # so remove the parameter as it should be handled by "dir" already
    keys_to_remove = []
    for key in consistent_params_data.keys():
        if type(consistent_params_data[key]) != list:
            if key in custom_conditions:
                keys_to_remove.append(key)
    
    for key in keys_to_remove:
        del consistent_params_data[key] 

    if len(dataset_dirs) > 0:
        dataset_param_meta["mapped_params"] = {"dir": {"value": [i for i in range(len(dataset_dirs))], "dir" : dataset_dirs }}
    else:
        dataset_param_meta["mapped_params"] = {}
        for dir_param_name, dir_param_values in dir_params.items():
            dataset_param_meta["mapped_params"][dir_param_name] = {"value": [i[1] for i in dir_param_values], "dir" : [i[0] for i in dir_param_values] }
    dataset_param_meta["index_params"] = consistent_params_data
    dataset_param_meta["internal_params"] = {"time": []}

    #if len(dataset_dirs) == 1:
    #    dataset_param_meta['mapped_params']['dir'] = {}
        
    write_dir(output_name, dataset_param_meta)
        

            

            
    
