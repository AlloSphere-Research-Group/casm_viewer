
import sys
import os

from tinc import *
parser = TincArgumentParser(description='Transfmat extraction')

parser.add_argument('__input_dir', type=str, default="./", nargs='?')
parser.add_argument('__output_names', type=str, default="transfmat", nargs='?')
parser.add_argument('__output_dir', type=str, default="", nargs='?')

args = parser.get_args()

import glob
import json
import os

json_files = glob.glob(args['__input_dir'] + "*.json")

found_transfmat = None

for open_file in json_files:
    found_transfmat = False
    with open(open_file) as f:
        if args['__verbose']==   True:
            print(f)
        try:
            j =json.load(f)
            if 'supercell' in j:
                found_transfmat = j['supercell']
        except:
            print("failed to parse json: " + open_file)
    if found_transfmat:
        break
write_name = args['__output_dir'] + args['__output_names'][0]

if os.path.exists(write_name):
    print('Found existing transfmat. Overwriting!')
    #write_name += '_alt'

new_mat_text = ''
if found_transfmat:
    for row in found_transfmat:
        for elem in row:
            new_mat_text += str(elem) + " "
        new_mat_text += '\n'

    print(new_mat_text) 

    with open(write_name, "w") as text_file:
        text_file.write(new_mat_text)
else:
    
    print("supercell not found in json file.")

