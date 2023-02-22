import sys
import argparse


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('input_filename')
parser.add_argument('output_filename')

args = parser.parse_args()

with open(args.output_filename, "w") as o:
    o.write("OUTPUT")
    
    """
    TODO:
    - Open and parse all GNPS2 Files if it hasn't been done on the first of the month
    - Save Output of Above by Date
    - Clean the above output (This will be the viewable portion)
    - Generate splits based on user input
    """