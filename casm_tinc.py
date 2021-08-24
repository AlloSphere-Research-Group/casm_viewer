
from tinc import *
import time

tclient = TincClient("host.docker.internal")

# allow time for things to settle. Will be fuxed in future
time.sleep(1)

# Parameters
settings_file_name = tclient.create_parameter(ParameterString, "settings_file_name", default_value="genetic_alg_settings.json")
hall_of_fame_index = tclient.create_parameter(ParameterInt,"hall_of_fame_index", default_value=0)
fit_dir = tclient.create_parameter(ParameterString,"fit_dir", default_value='')
#trigger_check = tclient.create_parameter(Trigger,"trigger_check")

# Processor 1
proc = ProcessorScript("check")
proc.command = "casm-learn"
# Capture the output of the command to file
proc.capture_output()
# Register paramters with processor. Changes trigger computation
proc.register_parameter(hall_of_fame_index)
proc.register_parameter(settings_file_name)
# Define the command line argument template
# names within '%%' that match parameter names will be replaced by their value
proc.set_argument_template("-s %%settings_file_name%% --checkhull --indiv %%hall_of_fame_index%%")

# The output will be managed by a diskbuffer, to update data everywhere
db = DiskBufferText("check_buffer", "check.0", "out/", "/shared")
tclient.register_disk_buffer(db)

# You can set a disk buffer to be the output of a Processor
proc.output_files = [db]

# Because we want to change directory and output file on every run,
# We define a 'prepare' function
def prepare_check(p):
    p.running_dir = fit_dir.value
    db.set_base_filename(f"check.{hall_of_fame_index.value}")
    print(f"Set output to: {p.output_files[0]}")
    #print(p._get_arguments())
    return True
# This function will be called right before calling the command for ProcessorScript
proc.prepare = prepare_check
proc.debug = True

# Alternatively, instead of doing a prepare function, you could update on parameter chagnes
# def set_output(value):
#     db.set_base_filename(f"check.{value}")
# hall_of_fame_index.register_callback(set_output)

# Processor 2
proc2 = ProcessorScript("select")
proc2.command = "casm-learn"
proc2.capture_output()
db2 = DiskBufferText("select_buffer","select_fit_eci.out", "out/", "/shared")
tclient.register_disk_buffer(db2)
proc2.set_argument_template("-s %%settings_file_name%% --select %%hall_of_fame_index%%")
proc2.register_parameter(hall_of_fame_index)
proc2.register_parameter(settings_file_name)

def prepare_select(p):
    p.running_dir = fit_dir.value
    return True


proc2.prepare = prepare_select
proc2.debug = True

# Processor 3
proc3 = ProcessorScript("fit_eci")
proc3.command = "casm"
db3 = DiskBufferText("fit_eci_buffer","select_fit_eci.out", "out/", "/shared")
tclient.register_disk_buffer(db3)
proc3.output_files = [db]
proc3.set_argument_template("query -k comp formation_energy hull_dist clex clex_hulldist -o %%:OUTFILE:0%%")
proc3.register_parameter(hall_of_fame_index)
proc3.register_parameter(settings_file_name)

def prepare_select(p):
    p.running_dir = fit_dir.value
    return True
proc3.prepare = prepare_select
proc3.debug = True

# Now keep the program running until ENTER
print("Press ENTER to exit")
input()

# close client to end threads
tclient.stop()