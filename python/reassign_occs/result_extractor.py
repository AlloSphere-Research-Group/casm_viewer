# -*- coding: utf-8 -*-
"""
Extracts POSCAR data from gzip files
for dynamic montecarlo data sets
"""

import gzip
import json
import os

from reassign_occs import parse_prim, determine_species, write_poscar
from template_creator import create_template_from_files, write_template


def num_time_steps(dataset_path = '/alloshare/vdv group/monte_with_traj_quick/',
            sub_dir = '',
            target_file = 'trajectory.json'):
    file_path = dataset_path + sub_dir + target_file + '.gz'
    with gzip.GzipFile(file_path, 'rb') as f:
        results_json = json.loads(f.read().decode('UTF-8'))
#        step_list = results_json['Step']
        
        return len(results_json['DoF'])
    return None

def extract(dataset_path = '/alloshare/vdv group/monte_with_traj_quick/',
            prim_path = '/alloshare/vdv group/monte_with_traj_quick/',
            sub_dir = '',
            time_step: int = 0,
            prim_name = 'prim_labels.json',
            target_file = 'trajectory.json'):
    file_path = dataset_path + sub_dir + target_file + '.gz'
    print("Opening: " + file_path)
    with gzip.GzipFile(file_path, 'rb') as f:
        results_json = json.loads(f.read().decode('UTF-8'))
#        step_list = results_json['Step']
        
        
        with open(prim_path + prim_name) as f:
            basis = parse_prim(f)
            print("Counting species")
            species_count, species_label = determine_species(basis,dof)
            print("Counting species done")
        return species_count, species_label
        
    
if __name__ == '__main__':

    import argparse
    
    parser = argparse.ArgumentParser(description='Graph generator')
    
    parser.add_argument('config_file', type=str, default="config.json", nargs='?')
    
    args = json.load(open(parser.parse_args().config_file))
    
    print(args)
     
    sub_dir = 'conditions.' + str(int(args['condition'])) + '/'
    
    print(args['dataset_path'] + 'prim.json')
    if os.path.exists(args['dataset_path'] + '/prim.json'):
        prim_name = 'prim.json'
    else:
        prim_name = 'prim_labels.json'
    
    # Create POSCAR template
    # template_name = 'template_POSCAR'
    # if not os.path.exists(args['dataset_path'] + template_name):
    #     lat, coords = create_template_from_files(args['prim_path'], args['dataset_path']  +'/cached_output/transfmat')
          
    #     write_template(lat, coords, args['dataset_path'] + template_name)
    
    import os
    if not os.path.exists(args['__output_dir']):
        os.makedirs(args['__output_dir'])
    
    # Write POSCAR
    species_count, species_label = extract(args['dataset_path'], args['dataset_path'], sub_dir, int(args['time_step']), prim_name)
    write_poscar(template_name, args['__output_dir'] + args['__output_names'][0], species_count, species_label)
    print('wrote ' + args['__output_dir'] + args['__output_names'][0])
