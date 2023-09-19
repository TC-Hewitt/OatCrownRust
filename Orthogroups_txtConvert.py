#!/usr/bin/env python

import argparse, csv, re

def main():

#parse arguments
    parser = argparse.ArgumentParser(description='converts Orthogroups.txt from OrthoFinder output into long format for use in gggenomes. Prints to STDOUT.')
    parser.add_argument('-i', '--infile', required=True, help='input Orthogroups.txt from OrthoFinder.')
    args = parser.parse_args()

#read infile
    infile = open(args.infile, 'r')
    reader = csv.reader(infile, delimiter = ' ', quoting=csv.QUOTE_NONE)

#transform rows
    print("cluster_id\tfeat_id")
    for row in reader:
        if len(row) >= 3:
            for i in range(1, len(row)):
                print(row[0].strip(":") + "\t" + row[i])
        else:
            print("unassigned\t" + row[1])

if __name__ == '__main__':
    main()
