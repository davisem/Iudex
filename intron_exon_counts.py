#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""
__modname__ = 'hasher.py'

import argparse
from pybedtools import BedTool


def parse_cmdline_params(cmdline_params):

    info = "Count gene trap insertions in samples"

    parser = argparse.ArgumentParser(description=info)

    parser.add_argument('-i', '--intron_bed', type=str, required=True,
                        help='Provide the path to a gene intron annotation'))

	parser.add_argument('-e', '--exon_bed', type=str, required=True,
						help='Provide the path to a gene exon annotation')

	parser.add_argument('-b', '--bed_file', type=str, required=True,
						help='Provide the path to the bed file for intersection')

	parser.add_argument('-o', '--output_file', type=str, require=True,
						help="Provide an output filename")

    return parser.parse_args(cmdline_params)


class InsertionTableBuilder(object):

	def __init__(self, intron_path, exon_path):
		self._intron_bedtool = BedTool(intron_path)
		self._exon_bedtool = BedTool(exon_path)

	def intersectIntronsAntisense(self, input_bedtool):
		input_bedtool.intersect(self._intron_bedtool, wa=True, S=True)

	def intersectIntronsAntisense(self, input_bedtools):
		input_bedtool.intersect(self._intron_bedtools, wa=True, S=True)

	def intersectExons(self, input_bedtool):
		input_bedtool.intersect(self._exon_bedtool, wa=True)

	def build(self, input_bed_path):








