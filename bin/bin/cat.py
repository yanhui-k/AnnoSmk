#!/usr/bin/env python3
import sys

def cat(input,output):
    file1=open(input,"r")
    file2=open(output,"w")
    lines=file1.readlines()
    for line in lines:
        if ("LOCUS" in line):
            line1=line.split("       ")[1]
            file2.write(line1)
            file2.write("\n")
    file1.close()
    file2.close()

cat(sys.argv[1], sys.argv[2])