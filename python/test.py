
# This script shows usage of TINC's argument parsing
# It allows creating scripts that can be called through command line
# arguments or through json configuration files

# You can use this script by calling it with command line arguments:
# python3 test.py aaa bbb
# where 'aaa' will be interpreted as arg1 and 'bbb' as arg2
# You can also create a json file with:
# {"__tinc_metadata_version": 0.1, "arg1": "aaa", "arg2": "bbb"}
# Then calling this script with:
# python3 test.py config.json

import sys, os
from tinc import *

parser = TincArgumentParser(description='POSCAR Template generator')

# By default only json configs are parsed. You can add arguments using
# the regular argparse functions. For example for positonal arguments use:

parser.add_argument('arg1', type=str, default="arg1_val", nargs='?')
parser.add_argument('arg2', type=str, default="arg2_val", nargs='?')

args = parser.get_args()

print(f'arg1 = {args["arg1"]}' )
print(f'arg2 = {args["arg2"]}' )
