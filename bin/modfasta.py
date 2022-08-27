#!/usr/bin/env python3
import os
import sys

def modfasta(fi = None, fo = None):
    with open(fi, 'r') as f1:
        data = f1.readlines()
        open(fo, "w")
        for line in data:
            if not line.startswith(">"):
                line = line.replace("\n", "")
                with open(fo, "a") as f2:
                    f2.write(line)
                    f2.close()
            else:
                line = line.replace(">", "\n>")
                with open(fo, "a") as f2:
                    f2.write(line)
                    f2.close()
                    
modfasta(sys.argv[1], sys.argv[2])