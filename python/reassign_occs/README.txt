reassign_occs notes

Jonas Kaufman
jlk@ucsb.edu

Inputs (supplied by VdV group):
- prim.json relabeled to specify the sites that we care about
- template_POSCAR: a POSCAR containing all possible sites, corresponding to the same supercell used for Monte Carlo. Can be generated using casm super from your PRIM with the same transformation matrix used for Monte Carlo. Site labels are not used, but the site order must be the same as the prim.json.
- final_state.json: casm Monte Carlo output file

Outputs
- altered_POSCAR: POSCAR corresponding to the final_state.json with sites in the same order as the template_POSCAR. Vacancies are included although we normally want to hide these.
- altered_POSCAR_sorted: optional sorted version of the altered_POSCAR

Usage possibilities:
- We may want to keep a single list of site positions from the template_POSCAR and then relabel for each configuration using only the the arrays species_label and species_count from the script. This would reduce the space required compared to storing a separate POSCAR for each data point.
- We could also use this to generate an altered_POSCAR_sorted for each data point and then visualize those

Extended example (directly from CASM output):
The data corresponding to the NaCoO2/P3_a_heating_12 and NaCoO2/P3_a_cooling datasets are on braid in these directories:
/home/jlkaufman/NaxCoO2/P3/monte/a/cooling
/home/jlkaufman/NaxCoO2/P3/monte/a/heating_12

The results.json are under mu_X and the final_state.json are under mu_X/conditions_Y 
I have also added a prim_labels.json and template_POSCAR in each so that you can generate labeled structures.
