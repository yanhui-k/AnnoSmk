#!/usr/bin/env python3
import sys
def reformat_to_work_with_MAKER(infile = None, outfile = None):
    infile = open(infile, "r")
    outfile = open(outfile, "w")
    i = 0
    for line in infile:
        if line[0] != "#":
            i = i+1
            outfile.write(line.replace("\n",";ID="))
            outfile.write(str(i))
            outfile.write("\n")
        else:
            outfile.write(line)
    infile.close()
    outfile.close()
reformat_to_work_with_MAKER(sys.argv[1], sys.argv[2])
