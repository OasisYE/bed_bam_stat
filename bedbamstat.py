# -*- coding: utf-8 -*-
import os
import time
import pandas as pd
import numpy as np
import logging
import pysam
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
        assert os.path.exists(bed_file)
        assert os.path.exists(bam_file)

        bedArray = self.bed2nparray(bed_file)
        samfile = pysam.AlignmentFile(bam_file)

        starts = list(map(int, bedArray[:, 1]))
        ends = list(map(int, bedArray[:, 2]))
        target_size = np.sum(np.array(ends) - np.array(starts))
        all_alignments = samfile.mapped + samfile.unmapped
        coverage_arr = []
        capture_alignments = 0

        for one_region in bedArray:
            ACGT = samfile.count_coverage(str(one_region[0]), int(one_region[1]), int(one_region[2]))
            coverage_arr.append(sum(np.array(ACGT)))
            capture_alignments += samfile.count(str(one_region[0]), int(one_region[1]), int(one_region[2]))

        # 捕获率
        capture_rate = float(capture_alignments) / all_alignments * 100  ##
        target_base = np.sum(list(map(np.sum, coverage_arr)))
        # 平均深度
        depth_in_target = float(target_base) / target_size
        # 覆盖度
        capture_size = sum(map(lambda x: len(x[x > 0]), coverage_arr))
        target_coverage = float(capture_size) / target_size * 100
        # 4X覆盖度
        X4_capture_size = sum(map(lambda x: len(x[x >= 4]), coverage_arr))
        target_4X = float(X4_capture_size) / target_size * 100
        # 10X覆盖度
        X10_capture_size = sum(map(lambda x: len(x[x >= 10]), coverage_arr))
        target_10X = float(X10_capture_size) / target_size * 100

        return 'TargetMeanDepth TargetCov TargetCov_4X TargetCov_10X CaptureRate\n{} {} {} {} {}'.format(depth_in_target,
                                                                                                     target_coverage,
                                                                                                     target_4X,
                                                                                                     target_10X,
                                                                                                     capture_rate)

    def onebaitqc(self, samfile, region):
        '''
        统计一个bait区间的平均深度，覆盖度，4X覆盖度，10X覆盖度
        :return:
        '''
        contig, starts, ends = str(region[0]), int(region[1]), int(region[2])
        target_size = np.sum(np.array(ends) - np.array(starts))
        baitid = '\-'.join(str(region[0]), str(region[1]), str(region[2]))

        coverage_arr = []

        ACGT = samfile.count_coverage(contig, starts, ends)
        coverage_arr = sum(np.array(ACGT))
        target_base = np.sum(list(map(np.sum, coverage_arr)))

        # 平均深度
        depth_in_target = float(target_base) / target_size
        # 覆盖度
        capture_size = sum(map(lambda x: len(x[x > 0]), coverage_arr))
        target_coverage = float(capture_size) / target_size * 100
        # 4X覆盖度
        X4_capture_size = sum(map(lambda x: len(x[x >= 4]), coverage_arr))
        target_4X = float(X4_capture_size) / target_size * 100
        # 10X覆盖度
        X10_capture_size = sum(map(lambda x: len(x[x >= 10]), coverage_arr))
        target_10X = float(X10_capture_size) / target_size * 100

        return baitid, depth_in_target, target_coverage, target_4X, target_10X

    def baitsqcbedbam(self, bam_file:str, bed_file:str, result_file:str):
        '''
        统计bed文件每个bait区间的平均深度，覆盖度，4X覆盖度，10X覆盖度
        :param result_file: 文件输出
        '''
        assert os.path.exists(bed_file)
        assert os.path.exists(bam_file)

        bedArray = self.bed2nparray(bed_file)
        samfile = pysam.AlignmentFile(bam_file)

        starts = list(map(int, bedArray[:, 1]))
        ends = list(map(int, bedArray[:, 2]))
        target_size = np.sum(np.array(ends) - np.array(starts))
        all_alignments = samfile.mapped + samfile.unmapped
        coverage_arr = []
        capture_alignments = 0

        out = []

        for one_bait in bedArray:
            out.append('\t'.join(self.onebaitqc(samfile, one_bait)))


        self.write_file(result_file, '\n'.join(out))










