"""
Clean reads
"""

from iga.apps.base import bsub, emain

import logging
import coloredlogs

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)

clean_novogene_sh = """
ADAPTER=/ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/adapters/novogene.fa
TAILCROP=145
HEADCROP=10
LEFT={0}
RIGHT={1}
java -jar  /ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/trimmomatic-0.38.jar PE \
-phred33 $LEFT $RIGHT \
$LEFT.clean.fq.gz $LEFT.clean.unpair.fq.gz \
$RIGHT.clean.fq.gz $RIGHT.clean.unpair.fq.gz \
ILLUMINACLIP:${{ADAPTER}}:0:30:10 \
LEADING:3 TRAILING:3 CROP:$TAILCROP HEADCROP:$HEADCROP SLIDINGWINDOW:1:10 MINLEN:75
"""


def clean_novogene(left=None, right=None):
    """
    Clean fastq for polish purpose
    :param left:
    :param right:
    :return:
    """
    cmd = clean_novogene_sh.format(left, right)
    bsub(cmd, name='Trimmomatic')
    return 0


# 0left
# 1right
# 2 adapter
clean_mgiseq_sh = """
ADAPTER=/ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/adapters/MGISeq.fa
TAILCROP=145
HEADCROP=10
LEFT={0}
RIGHT={1}
java -jar  /ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/trimmomatic-0.38.jar PE \
-phred33 $LEFT $RIGHT \
$LEFT.clean.fq.gz $LEFT.clean.unpair.fq.gz \
$RIGHT.clean.fq.gz $RIGHT.clean.unpair.fq.gz \
ILLUMINACLIP:${{ADAPTER}}:0:30:10 \
LEADING:3 TRAILING:3 CROP:$TAILCROP HEADCROP:$HEADCROP SLIDINGWINDOW:1:10 MINLEN:75
"""


def clean_mgiseq(left=None, right=None, source='bgi'):
    """
    Clean fastq for polish purpose
    :param left:
    :param right:
    :param source: bgi
    :return:
    """
    adapter_dict = {'bgi': '/ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/adapters/MGISeq.bgi.fa'}
    cmd = clean_mgiseq_sh.format(left, right, adapter_dict[source])
    bsub(cmd, name='Trimmomatic')
    return 0


# 0 Adatapter location
# 1 HeadCrop, eg 10
# 2 TailCrop, eg 145
# 3 Left fastq
# 4 right fastq
trimmomatic_head = r"""
ADAPTER={0}
HEADCROP={1}
TAILCROP={2}
LEFT={3}
RIGHT={4}
java -jar  /ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/trimmomatic-0.38.jar PE \
-phred33 $LEFT $RIGHT \
$LEFT.clean.fq.gz $LEFT.clean.unpair.fq.gz \
$RIGHT.clean.fq.gz $RIGHT.clean.unpair.fq.gz \
"""

trimmomatic_head_SE = r"""
ADAPTER={0}
LEFT={1}
java -jar  /ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/trimmomatic-0.38.jar SE \
-phred33 $LEFT  \
$LEFT.clean.fq.gz \
"""

trimmomatic_mode = {'polish': r"""ILLUMINACLIP:${{ADAPTER}}:2:30:10 \
LEADING:3 TRAILING:3 CROP:$TAILCROP HEADCROP:$HEADCROP SLIDINGWINDOW:1:10 MINLEN:75""",
                    'normal': r"""ILLUMINACLIP:${{ADAPTER}}:2:30:10 \
LEADING:3 TRAILING:3 CROP:$TAILCROP HEADCROP:$HEADCROP SLIDINGWINDOW:4:15 MINLEN:75""",
                    'smallrna': r"""ILLUMINACLIP:${{ADAPTER}}:2:30:10 \
SLIDINGWINDOW:4:15 AVGQUAL:25"""
#java -jar /opt/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 24 -phred33 \
                # FemaleMito1.fastq FemaleMito1_noadapters.fastq ILLUMINACLIP:adapters.ultimate.fa:2:30:10 AVGQUAL:25
                    }


def clean_fraser(left=None, right=None, source='hic', headcrop=10, tailcrop=145, mode='normal'):
    """
    Clean fastq for polish purpose. eg
    clean_fraser left.fq right.fq --source hic
    clean_fraser left.fq right.fq --source sgs

    :param left:
    :param right:
    :param source: hic/sgs/lncrna/ssrna/wgbs
    :param headcrop: default 10
    :param tailcrop: default 145
    :param mode: [normal|polish], normal sliding window is 4:15, polish sliding window is 1:10
    :return:
    """
    #TODO: small RNA clean is not available yet (need change SE)
    fraser_adapter = {
        'hic': '/ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/adapters/frasergen/HiC.fa',
        'lncrna': '/ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/adapters/frasergen/lncRNA.fa',
        'sgs': '/ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/adapters/frasergen/sgs.fa',
        'smallrna': '/ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/adapters/frasergen/smallRNA.fa',
        'ssrna': '/ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/adapters/frasergen/ssRNA.fa',
        'wgbs': '/ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/adapters/frasergen/WGBS.fa'
    }

    source = source.lower()
    if source not in fraser_adapter:
        logging.error("{} can't be found. Existing adapters are {}".format(source, str(fraser_adapter.keys())))
        exit(1)
    elif source == 'smallrna':
        mode = 'smallrna'
    if mode == 'smallrna':
        raw_cmd = trimmomatic_head_SE.format(fraser_adapter[source], headcrop, tailcrop, left)
        cmd = raw_cmd.format(fraser_adapter[source], headcrop, tailcrop, left)
    elif mode  in ['normal', 'polish']:
        raw_cmd = trimmomatic_head + trimmomatic_mode[mode]
        cmd = raw_cmd.format(fraser_adapter[source], headcrop, tailcrop, left, right)
    else:
        logging.error("mode could be only normal or polish")
        exit(1)

    bsub(cmd, name='Trimmomatic')
    return 0


if __name__ == "__main__":
    emain()
