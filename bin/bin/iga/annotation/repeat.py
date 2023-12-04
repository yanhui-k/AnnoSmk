"""
Repeat annotation utils
"""
import logging
import os

import os.path as op

from iga.apps.base import emain, conda_act, bsub, get_prefix, sh, waitjob


# def repeat_mask():


def clean_fasta_name(fasta=None):
    """
    sed "s/|arrow_np1212//; s/:::fragment//; s/:::debris//" elumb.contig.fa.origin_name > elumb.contig.fa
    :return:
    """
    sedcmd = """sed "s/|arrow_np1212//; s/:::fragment//; s/:::debris//" {0} > {0}.mod""".format(fasta)
    stdout = sh(sedcmd)
    print(stdout)


# 0 ref.fa
# 1 repeat.gff
# 2 threads
repeat_masker_sh = r"""
if [ ! -e {0}-families.fa ];
then
    BuildDatabase -name {0} -engine ncbi {0} 
    RepeatModeler -engine ncbi -pa {1} -database {0}
fi
mkdir -p custom_lib.out
RepeatMasker -lib {0}-families.fa {0} -pa {1} -dir custom_lib.out
"""


# species Viridiplantae
def repeatmasker(genome=None, species='Viridiplantae', denovo='T', threads=30):
    """
    python -m iga.annotation.repeat repeatmasker --species Viridiplantae --denovo F sesep.genome.fa
    python -m iga.annotation.repeat repeatmasker --species '' --denovo T sesep.genome.fa
    Default is custom annotation only
    :param genome:
    :param species: like Viridiplantae
    :param denovo:
    :return:
    """

    cmd = goto_workdir('repeatmask', sample=genome)
    genome = "../{}".format(genome)
    if species != '':
        cmd += "\nmkdir -p species_lib.out && RepeatMasker {0} -species {1} -pa {2} -dir species_lib.out\n".format(
            genome, species, threads
        )
    elif denovo == 'T':
        cmd += repeat_masker_sh.format(genome, threads)
    else:
        logging.error("Either provide a species name or use denovo prediction mode")
        exit(1)
    prefix = get_prefix(genome)
    jobid = bsub(cmd, name="repeat_masker_{}".format(prefix), cpus=threads)
    waitjob(jobid)


post_repeatmasker_sh = r"""
species=Viridiplantae
cd {0}
mkdir -p Full_mask
## unzip
gunzip *lib.out/*.cat.gz
cat *lib.out/*.cat >Full_mask/full_mask.cat
## to mask.out
cd Full_mask
ProcessRepeats -species $species full_mask.cat
## create GFF3
rmOutToGFF3.pl full_mask.out > full_mask.gff3

## isolate complex repeats
grep -v -e "Satellite" -e ")n" -e "-rich" full_mask.gff3 \
     > full_mask.complex.gff3

## reformat to work with MAKER
cat full_mask.complex.gff3 | \
 perl -ane '$id; if(!/^\#/){{@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"}} print $_' \
         > full_mask.complex.reformat.gff3

echo "Repeat GFF file is located in "
echo "$PWD/full_mask.complex.reformat.gff3"

cd ../..
# Step3. Get soft masked genome
"""


# 0 workdir
# 1 genome.fasta
# Output genome.masked.fasta
get_mask_fasta_sh = r"""
bedtools maskfasta -soft -fi {1} -bed {0}/Full_mask/full_mask.complex.reformat.gff3 \
-fo {1}.masked.fa
"""


def post_repeatmasker(dir=None, genome=None):
    """
    Post repeat process
    :param dir:
    :return:
    """
    cmd = post_repeatmasker_sh.format(dir)
    if genome != '':
        cmd += get_mask_fasta_sh.format(dir, genome)
    jobid = bsub(cmd, name='post_repeatmasker')
    waitjob(jobid)


def goto_workdir(program, sample=''):
    """
    :param program:
    :param sample:
    :return:
    """
    folder_name = 'workdir_{}_{}'.format(program, sample)
    cmd = "mkdir -p {0}; cd {0}\n".format(folder_name)
    return cmd


def edta(genome=None, threads=20):
    """
    :param genome:
    :return:
    """
    abs_genome = op.abspath(genome)
    rel_genome = op.basename(genome)
    cmd = conda_act.format('EDTA')
    cmd += goto_workdir('edta', op.basename(rel_genome))
    cmd += 'EDTA.pl --anno 1  --genome {0} --threads {1} \
    >edta_anno.out 2>edta_anno.err'.format(abs_genome, threads)
    bsub(cmd, name='EDTA{}'.format(rel_genome), cpus=threads)


if __name__ == "__main__":
    emain()
