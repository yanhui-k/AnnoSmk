"""
Ortholog calculation related utils
"""
from iga.apps.base import emain, bsub, sh


def ks_dist(prefix1=None, prefix2='', threads=40, use_grid='T', type='prot'):
    """
    :param prefix1: like corne
    :param prefix2: like cusat
    :param threads: 40 default
    :param use_grid: [T/F]
    :param type: [prot|nucl]
    :return: 0
    """
    if prefix2 == '':
        prefix2 = prefix1
    cmd = """
python -m jcvi.compara.catalog ortholog --dbtype {0} {1} {2}
cut -f1,2 {1}.{2}.anchors |grep -v "#" > {1}.{2}.ortho
"""
    if prefix1 == prefix2:
        cmd += """ln -s {1}.cds {1}.{2}.cds; ln -s {1}.pep {1}.{2}.pep"""
    else:
        cmd += """cat {1}.cds {2}.cds > {1}.{2}.cds; cat {1}.pep {2}.pep > {1}.{2}.pep"""
    cmd = cmd.format(type, prefix1, prefix2)
    sh(cmd)
    ortho_file = ".".join([prefix1, prefix2, 'ortho'])
    cds_total = ".".join([prefix1, prefix2, 'cds'])
    pep_total = ".".join([prefix1, prefix2, 'pep'])
    kaks(ortho_file, cds_total, pep_total, threads=threads, use_grid=use_grid)
    return 0


# 0 ortholog
# 1 cds fasta
# 2 pep fasta
# 3 number of threads
kaks_sh = """
#Presequitence 
#1. ParaAT 
#2. KaKsCalculator
#Note:
#gene naming like following will fail
#>genea Orsat
#ParaAT will read as geneaOrsat by default

ORTHO={0}
CDS={1}
PEP={2}
echo {3} > proc
ParaAT.pl -h $ORTHO -n $CDS -a $PEP -p proc -o ParaAT.out -f axt -k
join_kaks.pl ParaAT.out/*.kaks >kaks_result
"""


def kaks(ortho=None, cds=None, pep=None, threads=40, use_grid='T'):
    """
    :param ortho: eg CORNE00006074-t2        CsaV3_1G039430
    :param cds: cds fasta
    :param pep: pep fasta
    :param use_grid: [T/F]
    :return:
    """
    cmd = kaks_sh.format(ortho, cds, pep, threads)
    if use_grid == 'T':
        bsub(cmd, name="kaks", cpus=threads)
    else:
        sh(cmd)


if __name__ == "__main__":
    emain()
