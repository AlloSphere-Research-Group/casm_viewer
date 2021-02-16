# Not needed if tinc-python is installed
import sys
sys.path.append('../external/tinc/tinc-python/tinc-python')
from tinc_client import *
tclient = TincClient()
tclient.wait_for_server_available()

from parameter import *

eci1_param=tclient.create_parameter(Parameter,\
                                    "tet_oct_eci","casm",-0.375,0.375,\
                                    [-0.375,-0.125,0.125,0.375],-0.375)

eci2_param=tclient.create_parameter(Parameter,\
                                    "oct_tet_NN","casm",0.0,6.0,[0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,6.0],6.0)

eci3_param=tclient.create_parameter(Parameter,\
                                    "oct_oct_NN","casm",0.0,1.0,\
                                    [0.0,1.0],0.0)

eci4_param=tclient.create_parameter(Parameter,\
                                    "tet_tet_NN","casm",0.0,0.5,[0.0,0.5],0.0)

root_folder = "C:/Users/Andres/source/repos/vdv_data/visualization/"
# root_folder = "D:/data_for_visualization/visualization/"

import shutil
def create_dir_string_from_eci_param(eci_value):
    prefix="AMX2_spinel_diffusion_0.0_0.0"+\
    "_"+str(eci1_param.value)+\
    "_"+str(eci2_param.value)+\
    "_"+str(eci3_param.value)+\
    "_"+str(eci4_param.value)
    datasetname=root_folder+prefix+\
    "/kinetic_mc"
    tclient.get_parameter("dataset").value=datasetname
    voltage_buffer=tclient.get_disk_buffer('graph1')
    voltage_buffer.allow_cache(5)
    fname=voltage_buffer.get_filename_for_writing()
    print(fname)
    # shutil.copyfile("C:/Users/askje/OneDrive/Documents/flow_chart_pngs/"+prefix+"flow_chart_voltage.png",fname)
    # voltage_buffer.done_writing_file(fname)
    # diffusion_buffer=tclient.get_disk_buffer('graph2')
    # diffusion_buffer.allow_cache(5)
    # fname=diffusion_buffer.get_filename_for_writing()
    # shutil.copyfile("C:/Users/askje/OneDrive/Documents/flow_chart_pngs/"+prefix+"flow_chart_diffusion.png",fname)
    # diffusion_buffer.done_writing_file(fname)
    # corr_fact_buffer=tclient.get_disk_buffer('graph3')
    # corr_fact_buffer.allow_cache(5)
    # fname=corr_fact_buffer.get_filename_for_writing()
    # shutil.copyfile("C:/Users/askje/OneDrive/Documents/flow_chart_pngs/"+prefix+"flow_chart_corrfact.png",fname)
    # corr_fact_buffer.done_writing_file(fname)
    return

eci1_param.register_callback(create_dir_string_from_eci_param)

eci2_param.register_callback(create_dir_string_from_eci_param)

eci3_param.register_callback(create_dir_string_from_eci_param)

eci4_param.register_callback(create_dir_string_from_eci_param)

eci1_param.value=0.125
eci2_param.value=6.0
eci3_param.value=0.0
eci4_param.value=0.0
tclient.wait_for_server_available()
import time
time.sleep(1);

import random
#import matplotlib.pyplot as plt
import threading
import shutil
imageBuffer = tclient.get_disk_buffer('graph')
imageBuffer.allow_cache(5)
image_names=["oct_doub_va_cleared_path",
            "oct_doub_va_double_cleared_path",
             "oct_doub_va_queued_path",
             "oct_no_va",
             "oct_sing_va",
             "oct_trip_va_forced_path",
             "oct_trip_va_free_range_path",
             "oct_trip_va_semiforced_path",
             "oct_trip_va_trapped_path",
             "tet_no_va",
             "tet_sing_va"]
def update_shellsite_pic(shellSite):
    #print("Shell bitstring " + str(shellSite))
    b = int(shellSite)
    current_image_name=""
    for iname in image_names:
        if b & 1 == 1:
            current_image_name=iname
            break
        b = b >> 1
    with threading.Lock():
        fname = imageBuffer.get_filename_for_writing()
        print(fname)
        shutil.copyfile("C:/Users/askje/OneDrive/Documents/"+current_image_name+".png",fname)
        imageBuffer.done_writing_file(fname)

        #print(current_image_name)

tclient.get_disk_buffer("trajectory_buffer").get_filename_for_writing()

shellsites = tclient.get_parameter("ShellSiteTypes")
shellsites.register_callback(update_shellsite_pic, False)



import random
#import matplotlib.pyplot as plt
import threading
import shutil
imageBuffer = tclient.get_disk_buffer('graph')
imageBuffer.allow_cache(5)
image_names=["oct_doub_va_cleared_path",
            "oct_doub_va_double_cleared_path",
             "oct_doub_va_queued_path",
             "oct_no_va",
             "oct_sing_va",
             "oct_trip_va_forced_path",
             "oct_trip_va_free_range_path",
             "oct_trip_va_semiforced_path",
             "oct_trip_va_trapped_path",
             "tet_no_va",
             "tet_sing_va"]
def update_shellsite_pic(shellSite):
    #print("Shell bitstring " + str(shellSite))
    b = int(shellSite)
    current_image_name=""
    for iname in image_names:
        if b & 1 == 1:
            current_image_name=iname
            break
        b = b >> 1
    with threading.Lock():
        fname = imageBuffer.get_filename_for_writing()
        print(fname)
        shutil.copyfile("C:/Users/askje/OneDrive/Documents/"+current_image_name+".png",fname)
        imageBuffer.done_writing_file(fname)

        #print(current_image_name)

tclient.get_disk_buffer("trajectory_buffer").get_filename_for_writing()

shellsites = tclient.get_parameter("ShellSiteTypes")
shellsites.register_callback(update_shellsite_pic, False)

shellsites.value = 20



param=tclient.parameter_spaces[0]
root_path=param.get_root_path()
trajectories_path=root_path+param.get_current_relative_path()+"\\trajectory.nc"
print(trajectories_path)

import xarray as xr
import numpy as np

def step_generator():
    root_path=param.get_root_path()
    trajectories_path=root_path+param.get_current_relative_path()+"\\trajectory.nc"
    print(f'processing: {trajectories_path}')
    try:
        chunkedby10ds=xr.open_dataset(trajectories_path,chunks={'sample_dim':100})
    except FileNotFoundError:
        return []
    #choose an atom to track
    #find the first step that site index changes, record new position
    #track that position from then on
    # repeat if changes. 
    previous_dof=np.array(chunkedby10ds['occupation_dofs'][0][:].values)
    atom_locs=np.where(previous_dof==1)[0].tolist()
    atom_to_track=atom_locs
    steps_changed=[[] for i in range(len(atom_locs))]
    swap_pairs=[[] for i in range(3001)]
    site_labels=[[atom_locs[i]] for i in range(len(atom_locs))]
    previous_dof=np.array(chunkedby10ds['occupation_dofs'][0][:].values)
    for step_count in range(1,3001):
        new_dof=np.array(chunkedby10ds['occupation_dofs'][step_count][:].values)
        diff_dof=previous_dof-new_dof
        #where diff_dof is 1 this is the site index where atom 
        #disappeared 1->0 
        disap_ix=np.where(diff_dof==1)[0][0]
        track_ind=atom_to_track.index(disap_ix)
        #time value is 1 more than step count
        steps_changed[track_ind].append(step_count+1)
        #255 means -1 this is the site index where atom 
        #appeared
        ap_ix=np.where(diff_dof==255)[0][0]
        atom_to_track[track_ind]=ap_ix
        swap_pairs[step_count]=[int(disap_ix),int(ap_ix)]
        site_labels[track_ind].append(atom_to_track[track_ind])
        previous_dof=new_dof
    return (steps_changed,swap_pairs)

#tclient.synchronize()
ps = tclient.get_parameter_space("casmParameters")
ps.enable_caching()


#ps.sweep(step_generator,params=[ps.get_parameter("dir")],force_values=True,dependencies=[eci1_param, eci2_param,eci3_param,eci4_param,ps.get_parameter('dir')])

[steps_changed,swap_pairs] = ps.run_process(step_generator,dependencies=[eci1_param, eci2_param,eci3_param,eci4_param,ps.get_parameter('dir')])

traj_buffer=tclient.get_disk_buffer('trajectory_buffer')
tclient.get_parameter("width", "trajectory").value = 0.5
tclient.get_parameter("alpha", "trajectory").value = 0.6
def write_positions(pos):
    import json
    fname = traj_buffer.get_filename_for_writing()

    with open(fname, 'w') as f:
        json.dump(pos, f)
        traj_buffer.done_writing_file(fname)

def correct_arrow(single_coord_list):
    single_arr=np.array(single_coord_list)
    elementary_hop_dist=2.2
    return -single_arr/np.linalg.norm(single_arr)*elementary_hop_dist


import xarray
template_path=ps.get_root_path() + tclient.get_processor("TemplateGenerator").output_files[0]
template_atoms=xarray.open_dataset(template_path)['atoms_var']
def calc_summed_trajectory(value):
    print("registered click")
    [steps_changed,swap_pairs] = ps.run_process(step_generator,dependencies=[eci1_param, eci2_param,eci3_param,eci4_param,ps.get_parameter('dir')])
    arrows=[[[35,35,35]]]
    for i in range(1,len(swap_pairs)):
        disap,ap=swap_pairs[i]
        first=template_atoms[disap].values.tolist()
        second=template_atoms[ap].values.tolist()
        vec=np.array([second[0]-first[0],second[1]-first[1],second[2]-first[2]])
        color=[0.2,0.8,0.0]
        if np.linalg.norm(vec) > 2.5:
            vec=correct_arrow(vec)
            color=[1,1,1]
        arrow=[vec.tolist()]
        arrows.append(arrow)
    print(np.array(arrows).shape)
    write_positions(arrows)
    print("finished")
    return

sum_traj_button=tclient.create_parameter(Trigger,"summed_trajectory","casm")
sum_traj_button.register_callback(calc_summed_trajectory)

time_param=tclient.get_parameter("time")
template_path=ps.get_root_path() + tclient.get_processor("TemplateGenerator").output_files[0]
import xarray
template_atoms=xarray.open_dataset(template_path)['atoms_var']
import time

full_stop=True
def individual_trajectories(v):
    global full_stop
    full_stop=True
    [steps_changed,swap_pairs] = ps.run_process(step_generator,dependencies=[eci1_param, eci2_param,eci3_param,eci4_param,ps.get_parameter('dir')])
    ixes_to_check=[]
    for i in range(len(steps_changed)):
        ixes_to_check.append(len(steps_changed[i]))
    arrows=[]
    for ind,check in enumerate(ixes_to_check):
        if check>0:
            for t in steps_changed[ind]:
                time_param.value=t
                disap,ap=swap_pairs[t-1]
                first=template_atoms[disap].values.tolist()
                second=template_atoms[ap].values.tolist()
                vec=np.array([second[0]-first[0],second[1]-first[1],second[2]-first[2]])
                color=[0.2,0.8,0.0]
                if np.linalg.norm(vec) > 2.5:
                    vec=correct_arrow(vec)
                    color=[1,1,1]
                final=np.array(first)+vec
                arrow=[first,final.tolist(),color]
                arrows.append(arrow)
                write_positions(arrows)
                time.sleep(0.2)
                if not full_stop:
                    full_stop=True
                    break
    return
indiv_traj_button=tclient.create_parameter(Trigger,"individual_trajectory","casm")
indiv_traj_button.register_callback(individual_trajectories)

clear_traj=tclient.create_parameter(Trigger,"clear_trajectories","casm")

full_stop=False

def clear_traj_buffer(value):
    global full_stop
    full_stop=False
    write_positions([])
clear_traj.register_callback(clear_traj_buffer)

def sel_atom(val):
    [steps_changed,swap_pairs] = ps.run_process(step_generator,dependencies=[eci1_param, eci2_param,eci3_param,eci4_param,ps.get_parameter('dir')])
    tclient.get_parameter("currentSelection").value=swap_pairs[steps_changed[val][0]-1][0]
    return
def reset_atom_selection_slider(value):
    [steps_changed,swap_pairs] = ps.run_process(step_generator,dependencies=[eci1_param, eci2_param,eci3_param,eci4_param,ps.get_parameter('dir')])
    all_possible_starts=[]    
    for ix,val in enumerate(steps_changed):
        if len(val)>0:
            all_possible_starts.append(ix)
    print(all_possible_starts)
    selector_slider=tclient.create_parameter(ParameterInt,"moving_atoms","casm",min(all_possible_starts),max(all_possible_starts),all_possible_starts)
    selector_slider.register_callback(sel_atom)
reset_button=tclient.create_parameter(Trigger,"reset_moving_atoms","casm")
reset_button.register_callback(reset_atom_selection_slider)

def calc_traj(value):
    [steps_changed,swap_pairs] = ps.run_process(step_generator,dependencies=[eci1_param, eci2_param,eci3_param,eci4_param,ps.get_parameter('dir')])
    print("grabbing selection")
    idx=tclient.get_parameter("moving_atoms").value
    print("grabbed")
    arrows=[]
    if len(steps_changed[idx])>0:
        for t in steps_changed[idx]:
            time_param.value=t
            disap,ap=swap_pairs[t-1]
            first=template_atoms[disap].values.tolist()
            second=template_atoms[ap].values.tolist()
            vec=np.array([second[0]-first[0],second[1]-first[1],second[2]-first[2]])
            color=[0.2,0.8,0.0]
            if np.linalg.norm(vec) > 2.5:
                vec=correct_arrow(vec)
                color=[1,1,1]
            final=np.array(first)+vec
            arrow=[first,final.tolist(),color]
            arrows.append(arrow)
            write_positions(arrows)
            time.sleep(0.2)
    return
calc_traj_buffer=tclient.create_parameter(Trigger,"calc_single_trajectory","casm")
calc_traj_buffer.register_callback(calc_traj,False)

currSel=tclient.get_parameter("currentSelection")
tclient.print()
tclient.synchronize()


for i,l in enumerate(steps_changed):
    if len(l)>4:
        print(i)
        break

activate_atom_disp.register_callback(calc_traj)

key_in = input("Press Enter to exit")
#time.sleep(100)
