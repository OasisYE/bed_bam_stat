# -*- coding: utf-8 -*-
import os
import time
import pandas as pd
import logging
logger = logging.getLogger(__name__)

DEFAULT_CONFIG = {
    'sep': '\t'
}


class BedBamStat():
    def __init__(self, config=None):
        self.sep = config.get('sep')

    def create_folder(self, folder):
        '''
        创建目录
        :param folder 需要创建的目录
        '''
        if os.path.exists(folder) is False:
            logger.debug('创建目录 %s' % folder)
            os.makedirs(folder)
        else: # pragma: no cover
            logger.debug('目录 %s 已存在' % folder)

    def write_file(self, path, data):
        '''
        写入文件
        '''
        with open(path, 'w') as f:
            f.write(data)

    def bed2nparray(self, bed_file):
        '''
        构造bed seg 数组区间
        :param bed_file bed文件
        :return:
            [['chrY' '2649862' '2651063']
            ...
            ['chrY' '28817931' '28817932']
        '''
        assert os.path.exists(bed_file)

        df = pd.read_csv(bed_file, sep=self.sep, header=None) # 暂时不支持有头注释的bed文件 TODO
        target_area = df.values
        return target_area

    def qcbedbam(self, bam_file:str, bed_file:str):
        '''
        统计bed文件整个区间的 平均深度，覆盖度，4X覆盖度，10X覆盖度，捕获率
        :param bam_file: bam文件，注意构建index
        :param bed_file: bed文件
        :return: print 屏幕输出
        '''
        bedArray = self.bed2nparray(bed_file)

