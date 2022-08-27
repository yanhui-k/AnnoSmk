#!/usr/bin/env python

import argparse
import textwrap
import os

def main():
    prog_name = "HYPHY absrel wrapper v1.0"
    usage = """Output positive selected genes like:
    Positive selected genes are:
    * Zm00008a028894_P01, p-value =  0.00000
    * Node1, p-value =  0.00000
    """

    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("cds", help="cds file")
    parser.add_argument("pep", help="pep file")
#    parser.add_argument("-f", "--flanking", default=10000, type=int, help="flanking distance default (1000)")
    args = parser.parse_args()  

    cds_file = args.cds
    pep_file = args.pep

    path = os.path.abspath(os.curdir) 
    if '/' not in cds_file:
        cds_file = os.path.join(path, cds_file)
        pep_file = os.path.join(path, pep_file)
    pep_align_file = pep_file + ".aln"
    cds_align_file = cds_file + ".aln"
    tree_file = cds_file + ".nwk"
    bs_file = cds_file + ".bs"
    hyphy_log_file = cds_file + ".hyphy_log"

#    flanking_distance = args.flanking
    os.system("muscle -in {} -out {} 2>/dev/null".format(pep_file, pep_align_file))
    os.system("pal2nal.pl -output fasta {} {}  > {} 2>/dev/null".format(pep_align_file, cds_file, cds_align_file))
    os.system("FastTree {} > {} 2>/dev/null".format(pep_align_file, tree_file))
    os.system("hyphy aBSREL --alignment {} --tree {} > {} 2>/dev/null". format(cds_align_file, tree_file, hyphy_log_file))

    flag = 0
    with open(hyphy_log_file) as fh:
        for line in fh:
            if(line.rstrip()== "### Adaptive branch site random effects likelihood test"):
                flag = 1
                print('Positive selected genes are:')
            if(flag == 1 and line.startswith('*')):
                print(line, end = '')

if __name__ == "__main__":
    main()
