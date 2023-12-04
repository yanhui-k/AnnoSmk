#!/usr/bin/env python

import argparse
import os
import shutil
import textwrap
import subprocess
from parse import *
import re
from collections import defaultdict
import iga.apps.cfg
import coloredlogs, logging
import os.path as op
import six
import sys
import time

from iga.apps import cfg
from iga.utils.natsort import natsorted

from rich.logging import RichHandler

"""
The basic library for iga, some of the functions are adapted from jcvi
"""

# Create a logger object.
logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)


## Some examples.
# logger.debug("this is a debugging message")
# logger.info("this is an informational message")
# logger.warning("this is a warning message")
# logger.error("this is an error message")
# logger.critical("this is a critical message")

def mv(oldfile, new_file):
    """
    gnu move, support mv file into directory/
    :param oldfile:
    :param new_file:
    :return:
    """
    return sh('mv {0} {1}'.format(oldfile, new_file))


def mkdir(dirname, overwrite=False):
    """
    Wraps around os.mkdir(), but checks for existence first.
    """
    if op.isdir(dirname):
        if overwrite:
            shutil.rmtree(dirname)
            os.mkdir(dirname)
            logging.debug("Overwrite folder `{0}`.".format(dirname))
        else:
            return False  # Nothing is changed
    else:
        try:
            os.mkdir(dirname)
        except:
            os.makedirs(dirname)
        logging.debug("`{0}` not found. Creating new.".format(dirname))

    return True


def split_fasta(fasta, workdir, chunk=100, bypart='F'):
    """
    :param fasta: fasta
    :param workdir: workdir
    :param chunk: number to be split
    :param bypart: default is F, if it's T, split each '>xx\nATG' into seperate file with fasta header as file name.
    :return:
    """
    if bypart == 'F':
        file_list = sh('split_fastav3.pl {0} {2} && mv {0}._ {1}'.format(
            fasta, workdir, str(chunk))).split()
    else:
        file_list = []
        # https://github.com/uditvashisht/split-fasta/blob/master/splitfasta/__main__.py
        with open(fasta, 'r') as f:
            data = f.read().split('>')
            for i, j in enumerate(data[1:], start=1):
                (header, content) = j.split('\n', 1)
                new_file_name = f'{header}.fasta'
                file_list.append(new_file_name)
                if os.path.exists(new_file_name):
                    continue
                with open(new_file_name, 'w') as f:
                    f.write('>' + j)
    return file_list


def get_prefix(name):
    """
    Input /fjds/g/AB.fds.txt
    return AB
    :param name:
    :return:
    """
    return op.basename(name).split('.')[0]


def sh(cmd, debug=False, parallel='F', cpus=1, warning='T'):
    """
    run command directly with subprocess.run
    :param cmd:
    :param warning: [T/F], wether to output warning with STDOUT
    :return:
    """
    ret = ''
    logger.info(cmd)
    prior_cmd = 'set -eo pipefail\n'
    if parallel == 'T':
        from multiprocessing import Pool
        if type(cmd) != list:
            cmd = cmd.split('\n')
        with Pool(int(cpus)) as p:
            logger.info(p.map(sh, cmd))
    else:
        try:
            ret = subprocess.check_output(prior_cmd + cmd, stderr=subprocess.STDOUT, shell=True).decode()
        except subprocess.CalledProcessError as cpe:
            # cpe.output.decode() +
            logging.error(cpe.output.decode())
            ret = cpe.returncode
        if warning == 'T':
            logger.warning(ret)
    return ret


def rscript(cmd):
    """
    run r code with rscripts
    :param cmd:
    :return:
    """
    logger.info(cmd)
    prior_cmd = 'set -eo pipefail\nRscript '
    rscript_sh = 'plot.' + str(time.time()).replace('.', '') + '.R'
    with open(rscript_sh, 'w') as fh:
        fh.write(cmd)
    cmd_full = prior_cmd + rscript_sh
    try:
        ret = subprocess.check_output(cmd_full,
                                      stderr=subprocess.STDOUT, shell=True).decode()
    except subprocess.CalledProcessError as cpe:
        logger.warning(cpe.output)
        ret = cpe.output
    logger.warning(ret)
    return ret


def qsub(cmd=None, cpus=1, name='output', sub=True, normal='F', node="rock0[12]"):
    """
    submit jobs via qsub
    :param cmd:
    :param cpus:
    :return:
    """
    newqsub = r"""#!/bin/bash
#
#$ -cwd
#$ -N {1}
#$ -j y
#$ -V
#$ -pe smp {0}
#$ -S /bin/bash
#$ -l h="{2}"

#set -exo pipefail
ROOT=$PWD
date
""".format(cpus, name, node)
    if normal == "T":
        newqsub = r"""#!/bin/bash
#
#$ -cwd
#$ -N {1}
#$ -j y
#$ -V
#$ -pe smp {0}
#$ -S /bin/bash

ROOT=$PWD
date
""".format(cpus, name)
    qsub_buff = newqsub + cmd
    qsub_sh = '.'.join(['qsub', name, str(time.time()).replace('.', ''), 'sh'])
    with open(qsub_sh, 'w') as fh:
        fh.write(qsub_buff)
    cmd_full = 'qsub ' + qsub_sh
    # Prepare finished, now submit
    logger.info(cmd_full)
    if sub:
        ret = subprocess.check_output(cmd_full, shell=True).decode()
    # try:
    #     logger.warning(ret)
    #     job_id = parse('Job <{}> is submitted to queue <' + queue + '>.', ret.rstrip())[0]
    #     logger.warning(job_id)
    # except TypeError:
    #     logger.error('submission failed for: {}'.format(cmd_full))
    return ret


def bsub(cmd, queue='Q104C512G_X4', direct_submit='F', cpus=1, name='', submit='T'):
    """
    submit jobs via bsub
    When using variable export in a cmd ,use direct_submit = 'F'
    :param cmd:
    :return:
    """
    bsub_cmd = 'bsub -R "span[hosts=1]" -q {0}  -o output.%J -e error.%J -n {1} '.format(queue, cpus)
    if name != '':
        bsub_cmd += "-J {} ".format(name)
    if direct_submit == 'T':
        prior_cmd = 'set -eo pipefail;'
        cmd_full = bsub_cmd + '"' + prior_cmd + cmd + '"'
    else:
        newbsub = r"""#!/bin/bash
set -exo pipefail
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/ActivePerl-5.24/bin:$PATH;
ROOT=$PWD
date
"""
        bsub_buff = newbsub + cmd
        rand_digit = str(time.time()).replace('.', '')[7:11] + str(time.time()).replace('.', '')[-2:]
        bsub_sh = ".".join(['bsub', name, rand_digit, 'sh'])
        with open(bsub_sh, 'w') as fh:
            fh.write(bsub_buff)
        cmd_full = bsub_cmd + '< ' + bsub_sh
    # Prepare finished, now submit
    logger.info(cmd_full)
    if submit != 'T':
        return 0
    # ret = subprocess.check_output(bsub_cmd + '"' + prior_cmd + cmd + '"', shell=True).decode()
    # Incase queue has trainling characters like -m 'node02'
    queue = re.sub(r'\s.*', '', queue)
    ret = subprocess.check_output(cmd_full, shell=True).decode()
    try:
        logger.info(ret)
        job_id = parse('Job <{}> is submitted to queue <' + queue + '>.', ret.rstrip())[0]
        # logger.info(job_id)
    except TypeError:
        logger.error('submission failed for: {}'.format(cmd_full))
    # queue jumper
    # sh('btop {}'.format(job_id))
    return job_id


def is_job_finished(joblist=None):
    """
    job list in submited jobs in list format
    will check whether job is finished,
    if finished return 1
    else return 0
    :param joblist:
    :return:
    """
    if type(joblist) == str:
        if ',' in joblist:
            joblist = joblist.split(',')
        else:
            joblist = [joblist]
    for j in joblist:
        status = sh("bjobs {}".format(j))
        # logger.warning(status)
        if re.search(r'{}  yitings DONE'.format(j), status) or \
                re.search(r'{}  yitings EXIT'.format(j), status) or \
                re.search(r'Job .* is not found', status):
            if 'EXIT' in status:
                logger.error("Job {} finished with error!".format(j))
            continue
        else:
            return 0
    return 1


def waitjob(joblist=None):
    """
    waitjob "101,102" or 101
    a while loop for check_job_status,
    won't leave unless all jobs are finished
    :param joblist:
    :return:
    """
    while not is_job_finished(joblist):
        time.sleep(10)
    return 1


conda_act = r"""
source ~/lh/anaconda3/etc/profile.d/conda.sh
conda activate {}
"""

workdir_sh = r"""
mkdir -p {0}
cd {0}
"""


###JCVI code start

def splitall(path):
    allparts = []
    while True:
        path, p1 = op.split(path)
        if not p1:
            break
        allparts.append(p1)
    allparts = allparts[::-1]
    return allparts


class ActionDispatcher(object):
    """
    jcvi
    This class will be invoked
    a) when the base package is run via __main__, listing all MODULESs
    a) when a directory is run via __main__, listing all SCRIPTs
    b) when a script is run directly, listing all ACTIONs
    This is controlled through the meta variable, which is automatically
    determined in get_meta().
    """

    def __init__(self, actions):

        self.actions = actions
        if not actions:
            actions = [(None, None)]
        self.valid_actions, self.action_helps = zip(*actions)

    def get_meta(self):
        args = sys.argv[0].split('/')[-3:]
        args[-1] = args[-1].replace(".py", "")
        if args[-2] == "iga":
            meta = "MODULE"
        elif args[-1] == "__main__":
            meta = "SCRIPT"
        else:
            meta = "ACTION"
        return meta, args

    def print_help(self):
        meta, args = self.get_meta()
        if meta == "MODULE":
            del args[0]
            args[-1] = meta
        elif meta == "SCRIPT":
            args[-1] = meta
        else:
            args[-1] += " " + meta

        help = "Usage:\n    python -m {0}\n\n\n".format(".".join(args))
        help += "Available {0}s:\n".format(meta)
        max_action_len = max(len(action) for action, ah in self.actions)
        for action, action_help in sorted(self.actions):
            action = action.rjust(max_action_len + 4)
            help += (
                    " | ".join((action, action_help[0].upper() + action_help[1:])) + "\n"
            )
        help += "\n"

        sys.stderr.write(help)
        sys.exit(1)

    def dispatch(self, globals):
        from difflib import get_close_matches

        meta = "ACTION"  # function is only invoked for listing ACTIONs
        if len(sys.argv) == 1:
            self.print_help()

        action = sys.argv[1]

        if not action in self.valid_actions:
            print("[error] {0} not a valid {1}\n".format(action, meta), file=sys.stderr)
            alt = get_close_matches(action, self.valid_actions)
            print(
                "Did you mean one of these?\n\t{0}\n".format(", ".join(alt)),
                file=sys.stderr,
            )
            self.print_help()

        globals[action](sys.argv[2:])


def dmain(mainfile, type="action"):
    cwd = op.dirname(mainfile)
    pyscripts = (
        [x for x in glob(op.join(cwd, "*", "__main__.py"))]
        if type == "module"
        else glob(op.join(cwd, "*.py"))
    )
    actions = []
    for ps in sorted(pyscripts):
        action = (
            op.basename(op.dirname(ps))
            if type == "module"
            else op.basename(ps).replace(".py", "")
        )
        if action[0] == "_":  # hidden namespace
            continue
        pd = get_module_docstring(ps)
        action_help = (
            [
                x.rstrip(":.,\n")
                for x in pd.splitlines(True)
                if len(x.strip()) > 10 and x[0] != "%"
            ][0]
            if pd
            else "no docstring found"
        )
        actions.append((action, action_help))

    a = ActionDispatcher(actions)
    a.print_help()


def get_module_docstring(filepath):
    "Get module-level docstring of Python module at filepath, e.g. 'path/to/file.py'."
    co = compile(open(filepath).read(), filepath, "exec")
    if co.co_consts and isinstance(co.co_consts[0], six.string_types):
        docstring = co.co_consts[0]
    else:
        docstring = None
    return docstring


def glob(pathname, pattern=None):
    """
    Wraps around glob.glob(), but return a sorted list.
    """
    import glob as gl

    if pattern:
        pathname = op.join(pathname, pattern)
    return natsorted(gl.glob(pathname))


### jcvi code end

class DictDb:
    """
    The data structure for storing ctf files
    Could be understand as a sub class of Config class
    Change log:
    0106: Allowing lines like <<include ideogram.conf>> to be read into dict, the values are None
    """

    def __init__(self):
        self.dictdb = defaultdict(dict)
        """
        whether this dictionary is nested (has sections) like 
        [abc]
        a=1
        b=2
        """
        self.has_section = False

    def update_val(self, key, val, section=''):
        """
        Add tag, value to dict
        """
        if section != '':
            self.dictdb[section][key] = val
            self.has_section = True
        else:
            self.dictdb[key] = val

    def get_dict_text(self, seperator='='):
        """
        print this dict to text
        Change in 0118: Allow <<index>> lines in config files, key was assigned lines while values are None
        """
        return_text = ''
        if self.has_section:
            # logger.debug('sections are {}'.format(self.dictdb.keys()))
            for section in self.dictdb:
                return_text += ('[{}]'.format(section) + "\n")
                # logger.debug('section is {}'.format(section))
                # logger.debug('section content is {}'.format(self.dictdb[section]))
                # logger.debug(self.dictdb[section].keys())
                for k in self.dictdb[section]:
                    if self.dictdb[section][k] is not None:
                        return_text += ('{}{}{}'.format(k, seperator, self.dictdb[section][k]) + "\n")
                    else:
                        return_text += k + "\n"
        else:
            for k in self.dictdb:
                if self.dictdb[k] is not None:
                    return_text += ('{}{}{}'.format(k, seperator, self.dictdb[k]) + "\n")
                else:
                    return_text += k + "\n"
        return return_text


class VersatileTable:
    """
    Store common ctl table like aa=bb using DictDb
    eg:
        vt = VersatileTable('[section]
        aa=bb
        ee=gg')
        vt.update('[section]aa=ee')
        vt.get_text()
        vt.write_to_file()
    """

    def __init__(self, content=''):
        self.seperator = '='
        self.dictdb = DictDb()
        self.content
        if content != '':
            self.content = content
            self.load()

    def load(self):
        """
        Preprocess config file, then call update to update all keys and values
        :return:
        """
        # Seperator for tag and value, like tag=value is default
        this_list = self.content.splitlines()
        section = ''
        for line in this_list:
            if line.rstrip() == '':
                # Skip blank lines
                continue
            elif line.startswith('['):
                # Finding section definition
                section = parse('[{}]', line)[0]
            elif line.startswith('#') or line.startswith(';'):
                # Finding section annotation
                section_annotation = re.sub(r'^#+', '', line)
            else:
                # Findng value assignments
                trimmed_line = re.sub('#.*', '', line).strip()
                if section != '':
                    trimmed_line = '[{}]{}'.format(section, trimmed_line)
                self.update(trimmed_line, use_semicolon_sep=False)

    def update(self, args, use_semicolon_sep=True):
        """
        Different from change_val in dictdb, this allows input like :
        "[general]genomesize=12M;[general]threads=11"
        """
        # logger.debug('use_semicolon_sep: {}'.format(use_semicolon_sep))
        if use_semicolon_sep:
            mylist = re.split(r'[\n;]', args)
        else:
            mylist = args.split('\n')
        for this_arg in mylist:
            if this_arg.strip() == '':
                continue
            section = ''
            content = this_arg
            if '[' in this_arg:
                (section, content) = parse('[{}]{}', this_arg)
            try:
                # logger.warning(content)
                (key, value) = content.split(self.seperator)
                key = key.strip()
                value = value.strip()
            except ValueError as e:
                # In case no seperator was found in the line
                # logger.debug("Treating {} and {} as non-seperator lines".format(this_arg))
                key = content.strip()
                value = None
            # if section if null charactor, will not used by dictdb.change_val
            self.dictdb.update_val(key=key, val=value, section=section)

    def get_text(self):
        """
        Return the table as text
        :return:
        """
        return self.dictdb.get_dict_text(self.seperator)

    def write_to_file(self, output_file):
        """Write the table to files"""
        with open(output_file, 'w') as fh:
            fh.write(self.dictdb.get_dict_text(self.seperator))


def grep(args, content=''):
    """
    GNU grep like function:
    eg: grep("-v ^#", content)
    :param args:
    :param content:
    :return:
    """
    keep = True
    key = args
    result = ''
    if '-v' in args:
        key = args.split()[1]
        keep = False
    for i in content.splitlines():
        if re.search(key, i) and keep or \
                not re.search(key, i) and not keep:
            result += i + "\n"
    return i


def debug(level=logging.DEBUG):
    """
    Turn on the debugging
    """
    #    borrowed from jcvi
    from rich.console import Console
    logging.basicConfig(
        level=level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(console=Console(stderr=True))],
    )


debug()


class Config(VersatileTable):
    """
    Eg:
        cfg = Config('maker')
        cfg.update('est_gff={};protein_gff={};rm_gff={}'.format(estgff, pepgff, rmgff))
        cfg.write_to_file("maker_round1.opts")
        cfg = Config('circos')
        cfg.update('karyotype={};plots.plot={};rm_gff={}'.format(estgff, pepgff, rmgff))
        cfg.write_to_file("maker_round1.opts")
    Usage:
        1. Give arguments, will return the specific CTLs object
        2. The object could be directly printed with get_text() or write to file with write_to_file("filename")
    Further explanation:
        A sub class inherated from Versatile Table that extend usage in following aspecects:
        1. Allow reading config template from cfg dictionary
        2. Allow handling html type tables (Currently only circos used this type, so not seperating it to another class)
    """

    def __init__(self, cfg_type):
        """
        :param cfg_type: the name of the config file, maybe falcon or maker; can also be the config file path
        :param format: 
        """
        # logger.debug(cfg_type)
        # logger.debug(cfg.cfg[cfg_type])
        self.multiple_section_list = ['highlight', 'plot', 'link', 'zoom', 'tick',
                                      'pairwise', 'rule', 'axis', 'background']
        self.format = ''
        try:
            self.content = cfg.cfg[cfg_type]
        except KeyError as e:
            logger.error("Unknown type of cfg file: {}, try running in external load".format(cfg_type))
            if op.exists(cfg_type):
                with open(cfg_type) as fh:
                    self.content = fh.read()
            else:
                logger.error("External config load error, exiting")
                return 1
        if 'circos' in cfg_type:
            # circos format is a bit weird, so handle it seperately
            self.format = 'circos'

        if self.format == '':
            super().__init__(self.content)
        elif self.format == 'circos':
            super().__init__()
            self.section_seperator = '-+-'
            self.load_circos()

    def load_circos(self):
        # Seperator for tag and value, like tag=value is default
        # logger.debug(self.content)
        self.insert_block(self.content)

    def update(self, args, use_semicolon_sep=True):
        """
        Different from change_val in dictdb, this allows input like :
        "[general]genomesize=12M;[general]threads=11"
        Or for circos:
        "plots.plot.color=1"
        """
        if self.format == '':
            super().update(args, use_semicolon_sep)
        elif self.format == 'circos':
            # like plots.plot.
            mylist = args.split(';')
            section = ''
            for this_arg in mylist:
                if re.search(r'{}'.format(self.section_seperator), this_arg):
                    section = this_arg.split(self.section_seperator)
                    key_value = section.pop()
                    (key, value) = key_value.split(self.seperator)
                    # if len(section) >= 1 and section[-1] in self.multiple_section_list:
                    #     #Automatically rename
                    #     self.multiple_section_list[section[-1]] += 1
                    #     section[-1] += '-{}'.format(self.multiple_section_list[section[-1]])
                    section = self.section_seperator.join(section)
                else:
                    content = this_arg
                    try:
                        (key, value) = content.split(self.seperator)
                        key = key.strip()
                        value = value.strip()
                    except ValueError as e:
                        # In case no seperator was found in the line
                        logger.warning("Treating {} and {} as non-seperator lines".format(this_arg, content))
                        key = content.strip()
                        value = None
                # logger.debug(this_arg)
                # logger.debug(content)
                # logger.debug('key{}val{}section{}'.format(key, value, section))
                if section != '':
                    key = "{}{}{}".format(section, self.section_seperator, key)
                self.dictdb.update_val(key=key, val=value)

    def insert_block(self, block):
        """
        foramt "<abc>a=b</abc>" to abc-+-a=b" where -+- is seperator
        :param block:
        :return:
        """
        # logger.debug(block)
        this_list = block.splitlines()
        section = []
        count_tag = defaultdict(int)
        for line in this_list:
            # logger.debug(line)
            if line.rstrip() == '':
                # Skip blank lines
                continue
            elif line.startswith('</'):
                # End of xx block
                # logger.debug(line)
                # logger.debug(section)
                section.pop()
            elif line.startswith('<') \
                    and not line.startswith('<<'):
                tag = parse('<{}>', line)[0]
                if tag in self.multiple_section_list:
                    # Those allow multiple section with same name, so I'll rename them
                    # so plots.plot will be automatically stored as plots.plot-1 for the first time.
                    count_tag[tag] += 1
                section.append("{}-{}".format(tag, count_tag[tag]))
            elif line.startswith('#') or line.startswith(';'):
                # Finding section annotation
                section_annotation = re.sub(r'^#+', '', line)
            else:
                if len(section) > 0:
                    line = self.section_seperator.join(section) + self.section_seperator + line
                # logger.debug(line)
                self.update(line)

    def get_text(self):
        """
        get_text:
        1. if the format is '', use VersatileTable.get_text()
        2. if the format if 'circos',
        :return:
        """
        if self.format == '':
            result_text = super().get_text()
        elif self.format == 'circos':
            result_text = ''
            last_section = []
            # logger.debug(self.dictdb.dictdb.keys())
            for param, val in self.dictdb.dictdb.items():
                if self.section_seperator in param:
                    this_section = param.split(self.section_seperator)
                    param = this_section.pop()
                else:
                    this_section = []
                    # First condition test on whether to print </tag> lines
                result_text += self.get_circos_section(last_section, this_section)
                last_section = this_section
                # logger.debug(self.dictdb.dictdb)
                if val is None:
                    # <<include ideogram.conf>> lines
                    result_text += "{}\n".format(param)
                else:
                    result_text += "{}{}{}\n".format(param, self.seperator, val)
            if len(last_section) > 0:
                result_text += self.get_circos_section(last_section, [])
        return result_text

    @staticmethod
    def get_circos_section(last_section, this_section):
        # logger.debug('last_section')
        # logger.debug(last_section)
        # logger.debug('this_section')
        # logger.debug(this_section)
        result_text = ''
        for last_index, last_s in enumerate(last_section):
            # logger.debug('printing last section')
            if len(this_section) > last_index and \
                    last_section[last_index] == this_section[last_index]:
                # "plots plots"
                pass
            else:
                # plots.plot1.rule1 plots.plot2
                # print(</rule1>\n</plot1>\n)
                for i in range(len(last_section) - 1, last_index - 1, -1):
                    last_section[i] = re.sub(r'[-\s].*', '', last_section[i])
                    result_text += "</{}>\n".format(last_section.pop(i))
                break
        # Second condition test on whether to print <tag> lines
        for this_index, this_s in enumerate(this_section):
            # logger.debug('printing this section')
            if len(last_section) > this_index and \
                    last_section[this_index] == this_section[this_index]:
                # "plots plots"
                pass
            else:
                # plots.plot1.rule1 plots.plot2.rule2
                # print(<plot2>\n<rule2>\n)
                for i in range(this_index, len(this_section)):
                    tmp = re.sub(r'[-].*', '', this_section[i])
                    result_text += "<{}>\n".format(tmp)
                break
        return result_text

    def write_to_file(self, output_file):
        """Write to files"""
        with open(output_file, 'w') as fh:
            fh.write(self.get_text())


def abspath_list(file_list):
    for i, v in enumerate(file_list):
        file_list[i] = op.abspath(v)
    return file_list


def fmain(func_name, args):
    """
    execute functions directly via command line interface
    :param func_name: the name of the function
    :param args: the args, usually sys.argv[2:]
    :return:
    """
    # parser = argparse.ArgumentParser(
    #     prog=prog_name,
    #     formatter_class=argparse.RawDescriptionHelpFormatter,
    #     description=textwrap.dedent(usage),
    #     epilog="")
    # parser.add_argument("GENOME", help="Genome to be evalutated in fasta format")
    # parser.add_argument("-t", "--threads", default=64, type=int, help="flanking distance default (1000)")
    # args = parser.parse_args()

    # p = argparse.ArgumentParser(prog=func_name, usage=func_doc)
    position_arg = []
    keyword_arg = {}
    number_args = 1
    import sys
    import inspect
    # 下面两行命令用于在函数内部得到函数的名称和文档
    # func_name = sys._getframe().f_code.co_name
    # func_doc = sys._getframe().f_code.co_consts[0]
    # 下面命令用于把字符串的函数名称转换成对象
    # func_name = 'isoseq_'
    # logger.debug(args)
    object_pointer = getattr(sys.modules['__main__'], func_name)
    p = argparse.ArgumentParser(prog=func_name, usage=object_pointer.__doc__)
    # 下面的两个命令用于从函数对象中调取形参的名字和默认值（空值用Nonetype表示），用来转换成parse_args
    if len(inspect.getfullargspec(object_pointer).args) == 0:
        # Incase a function dont' have an arg
        object_pointer()
        return 0
    for kw, kw_defaults in zip(inspect.getfullargspec(object_pointer).args,
                               inspect.getfullargspec(object_pointer).defaults):
        if kw_defaults == None:
            position_arg.append(kw)
        else:
            keyword_arg[kw] = kw_defaults
    if len(position_arg) == 1:
        # If only one input arg is needed for the function, allow multiple files as input
        number_args = '+'
    for k in position_arg:
        p.add_argument(k, help=k, nargs=number_args)
    for k, v in keyword_arg.items():
        p.add_argument("--" + k, default=v, help="default: '%(default)s'")

    real_arg = p.parse_args(args)

    # Results for storing arguments after running parse_args
    position_result = []
    keyword_result = {}

    for k in keyword_arg:
        # logger.debug("{}\n{}\n{}".format(k, real_arg, getattr(real_arg, k)))
        keyword_result[k] = getattr(real_arg, k)
        # print(keyword_result)

    for k in position_arg:
        k_arg = getattr(real_arg, k)
        if number_args == '+' and len(k_arg) > 1:
            # passing spaced args as a list
            position_result.append(k_arg)
        else:
            position_result.append(k_arg[0])

    # used to debug
    # logger.debug(position_result)
    # logger.debug(keyword_result)
    object_pointer(*position_result, **keyword_result)

    # if(number_args == 1):
    #     position_result = getattr(p, position_arg[0])
    #     object_pointer(**position_result, **keyword_result)
    # else:
    #     for k in position_arg:
    #         position_result.append(getattr(p, k))
    # object_pointer(p.__dict__)


def emain():
    """
    A common wrapper for all main function of xx.py under iga
    :return:
    """
    # iga_prior = 'iga v1.0'
    # actions are available functions, like [maker, isoseq]
    actions = []
    # this list are available functions with function pointer like [[maker, pointerxxxx], [isoseq, pointerxxx]]
    actions_with_help = []
    from inspect import getmembers, isfunction
    functions_list = [o for o in getmembers(sys.modules['__main__']) if isfunction(o[1])]
    # logger.warning(functions_list)
    for f in functions_list:
        if f[1].__module__ == "__main__" and f[0] != 'main':
            actions.append(f[0])
            actions_with_help.append([f[0], str(f[1].__doc__).replace('\n', ' ').lstrip()[:50]])
    # actions = ['isoseq', 'fastq2gff', 'isoseq_pb', 'maker_round1']
    if len(sys.argv) > 1 and sys.argv[1] in actions:
        action = sys.argv[1]
        if len(sys.argv) > 2:
            args = sys.argv[2:]
        else:
            args = []
        fmain(action, args)
    else:
        # Print help, from jcvi
        help_print = "Usage:\n    python -m {0}\n\n\n".format(sys.argv[0].replace('.py', ''))
        help_print += "Available {0}s:\n".format('action')
        max_action_len = max(len(action) for action, ah in actions_with_help)
        for action, action_help in sorted(actions_with_help):
            action = action.rjust(max_action_len + 4)
            help_print += (
                    " | ".join((action, action_help.title())) + "\n"
            )
        help_print += "\n"
        print(help_print)


if __name__ == "__main__":
    emain()
