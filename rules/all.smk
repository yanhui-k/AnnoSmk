import glob
import re
import os 
import sys
from snakemake.exceptions import WorkflowError

def get_output():
    out = []
    out.append("AnnoSmk/download.log")
    out.extend(expand("{transcriptome}",transcriptome=config["TRAN"]))
    if config["repeat_gff"]:
        out.extend(expand("AnnoSmk/{PREFIX}/evidence/repeat.gff",PREFIX=PREFIX))
    if config["total_est_gff"]:
        out.extend(expand("AnnoSmk/{PREFIX}/evidence/total_est.gff",PREFIX=PREFIX))
    if config["total_pep_gff"]:
        out.extend(expand("AnnoSmk/{PREFIX}/evidence/genblast.gff",PREFIX=PREFIX))
    if int(config["run_maker_Round"]) < 2:
        raise(WorkflowError("the 'run_maker_Round' parameter is < 2"))
    else:
        a = int(1)
        b = int(config["run_maker_Round"])
        for i in range(a,b):
            out.extend(expand("AnnoSmk/{PREFIX}/R{round}/total_master_datastore_index.log",PREFIX=PREFIX,round=i))
            out.extend(expand("AnnoSmk/{PREFIX}/R{round}/ref.fa",PREFIX=PREFIX,round=i))
            out.extend(expand("AnnoSmk/{PREFIX}/R{round}/{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm",PREFIX=PREFIX,round=i))
            out.extend(expand("AnnoSmk/{PREFIX}/R{round}/autoAug/autoAugPred_hints/shells",PREFIX=PREFIX,round=i))
#            out.extend(expand("AnnoSmk/{PREFIX}/R{round}/total.all.maker.transcripts.fasta.busco.embryophyta",PREFIX=PREFIX,round=i))
            out.extend(expand("AnnoSmk/{PREFIX}/R{round}/AED.csv",PREFIX=PREFIX,round=i))
            out.extend(expand("AnnoSmk/{PREFIX}/R{round}/clean_memory.log",PREFIX=PREFIX,round=i))
        out.extend(expand("AnnoSmk/{PREFIX}/R{round}/total_master_datastore_index.log",PREFIX=PREFIX,round=b))
        out.extend(expand("AnnoSmk/{PREFIX}/R{round}/ref.fa",PREFIX=PREFIX,round=b))
#        out.extend(expand("AnnoSmk/{PREFIX}/R{round}/total.all.maker.transcripts.fasta.busco.embryophyta",PREFIX=PREFIX,round=b))
        out.extend(expand("AnnoSmk/{PREFIX}/R{round}/AED.csv",PREFIX=PREFIX,round=b))
        out.extend(expand("AnnoSmk/{PREFIX}/R{round}/clean_memory.log",PREFIX=PREFIX,round=b))
        out.extend(expand("AnnoSmk/{PREFIX}.gff",PREFIX=PREFIX))
#    print(out)
    return out