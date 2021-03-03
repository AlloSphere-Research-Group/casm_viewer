#!/bin/env/python

# Sanjeev Kolli and Jonas Kaufman
# May 17, 2018
# Script to label sites in Monte Carlo supercells

from __future__ import print_function

import numpy as np
import json
import os

# read prim
def parse_prim(file_obj):
    prim=json.load(file_obj)
    basis=prim["basis"]
    return [ x["occupant_dof"] for x in basis]

def determine_species(basis,dof):
    # go through final state occupations and label
    vol=len(dof)/len(basis)
    species_map = []

    for i in range(len(dof)):
        # species_map[i]=basis[int(i/vol)][dof[i]]+str(int(i/vol))
        species_map.append(basis[int(i/vol)][dof[i]])

    return species_map

def write_poscar(template_path, output_file, species_count, species_label):
    # write labeled POSCAR
    with open(template_path,'r') as f:
        poslines=f.readlines()
    poslines[5] = " ".join(species_label)+ '\n'
    poslines[6] = " ".join(map(str,species_count)) + '\n'
    with open(output_file,'w') as f:
        f.writelines(poslines)

if __name__ == '__main__':
    
    import sys, os
    sys.path.append(
        os.path.abspath(os.path.join(os.path.dirname(__file__), "../../external/tinc/tinc-python/tinc-python")))
    from process_args import *
    
    parser = TincArgumentParser(description='POSCAR Template generator')

    # Add arguments for regular command line usage 
    parser.add_argument('prim_path', type=str, default="../prim.json", nargs='?')
    parser.add_argument('transfmat', type=str, default="transfmat", nargs='?')
    parser.add_argument('final_state_path', type=str, default="", nargs='?')
    parser.add_argument('__input_names', type=str, default="final_state.json", nargs='?')
    parser.add_argument('__output_names', type=str, default="template_POSCAR", nargs='?')
    parser.add_argument('__output_path', type=str, default="", nargs='?')
    
    args = parser.get_args()
    
    # input files
    
    output_name = args['__output_names'][0]
    
    if not os.path.exists(args['__output_dir']):
        os.makedirs(args['__output_dir'])
        
    with open(args['prim_path'],'r') as f:
        basis=parse_prim(f)

    # read final_state
    with open(args['__input_names'][0],'r') as f:
        dof=json.load(f)["occupation"]

    # sort POSCAR by species
    if False: # Make true to generate sorted
    #    print("Output:" + output_path + output_name)
        import poscar   # poscar.py (should be supplies)
        altered_pos=poscar.Poscar(args['template_path'] + 'sorted_' + output_name )
        altered_pos.write(args['template_path'] + output_name)

    # NetCDF ######
    import netCDF4

    out_name = args['__output_dir'] + output_name
    print("Using template " + args['template_path'])

    templatencfile = netCDF4.Dataset(args['template_path'], mode='r') 
    atoms_var = templatencfile.variables['atoms_var']
    template_atoms_var = atoms_var[:] # get numpy array for speed

    import shutil
    shutil.copyfile(args['template_path'], out_name)

    ncfile = netCDF4.Dataset(out_name, mode='w')
    atom = np.dtype([('basis_index',np.uint8),('occupant_dof',np.uint8)])
    atom_t = ncfile.createCompoundType(atom,'atom')
    index = ncfile.createDimension('index')
    output_atoms_var = ncfile.createVariable('occupancy',atom_t,('index'))

    vol=len(dof)/len(basis)
    
    output_atoms_var[len(template_atoms_var) - 1] = (0,0)
    for i,d in enumerate(dof):
        output_atoms_var[i] = (int(i/vol), dof[i])

    ncfile.close()
    templatencfile.close()
