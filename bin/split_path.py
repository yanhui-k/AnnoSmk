#!/usr/bin/env python3
import os
import sys

fullpath = sys.argv[1]

def split_path(fullpath = None):
    '''
    用于分离文件路径和文件名的脚本
    :param fullpath: 需要分离的完整路径
    :return:文件路径

    '''
    (filepath, filename) = os.path.split(fullpath)
    print(filepath)

split_path(fullpath)