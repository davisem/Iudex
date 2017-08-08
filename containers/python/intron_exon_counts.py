#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""

import sys
import argparse
import pandas as pd
from pybedtools import BedTool
import math



def parse_cmdline_params(cmdline_params):

    info = "Count gene trap insertions in samples"

    parser = argparse.ArgumentParser(description=info)

    parser.add_argument('-i', '--intron_bed', type=str, required=True,
                        help='Provide the path to a gene intron annotation')

    parser.add_argument('-e', '--exon_bed', type=str, required=True,
                        help='Provide the path to a gene exon annotation')

    parser.add_argument('-b', '--bed_file', type=str, required=True,
                        help='Provide the path to the bed file for intersection')

    parser.add_argument('-o', '--output_file', type=str, required=True,
                        help="Provide an output filename")

    return parser.parse_args(cmdline_params)


class InsertionRecord(object):
    
    FIELD_NAMES = ['chrom', 'start', 'insertion_orientation', 'gene']
    GENE_INDEX = 9
    GENE_STRAND_IDX = 11
    PLUS_STRAND = '+'
    
    ORIENTATION_SENSE_KEY = 'SENSE'
    ORIENTATION_ANTISENSE_KEY = 'ANTISENSE'

    __slots__ = ['chrom', 'start', 'orientation', 'gene']

    def __init__(self, bedtool_record):
        """
        Init method for class
        :param BedTool: A line in a BedTool object
        """
        self.chrom = bedtool_record.chrom
        self.start = self.getInsertionStart(bedtool_record)
        self.orientation = self.getOrientation(bedtool_record)
        self.gene = bedtool_record.fields[self.GENE_INDEX]

    def toIterable(self):
        """
        Retruns an iterable representation of the record
        :return iter(str): A string generator of the record
        """
        return (self.chrom, self.start, self.orientation, self.gene)

    def getOrientation(self, bedtool_line):
        """
        Gets the orientation of the gene_trap insertion with resepect to the gene
        :return str: The string key of the orientation.
        """
        return self.ORIENTATION_SENSE_KEY if bedtool_line.strand == bedtool_line.fields[self.GENE_STRAND_IDX] else self.ORIENTATION_ANTISENSE_KEY

    def getInsertionStart(self, bedtool_line):
        """
        Gets the start coordinate of the gene_trap insertion
        :return int: The starting coordinate relative to the gene trap.
        """
        return bedtool_line.start if bedtool_line.strand == self.PLUS_STRAND else bedtool_line.stop

    def __repr__(self):
        """
        Repr of this class
        :return str:
        """

        return "<{}>:{},{},{},{}".format(self.__class__.__name__, self.chrom, self.start, self.orientation, self.gene)

    def __eq__(self, other):
        """
        Defines equalit:y behavior of this class
        :return boolean: Is this equal to another record.
        """

        return all([self.chrom == other.chrom, self.start == other.start, self.orientation == other.orientation, self.gene==other.gene])

    def __hash__(self):
        """
        Implements hashing behavior of this class
        :return int: Hash value of this class.
        """

        return hash(self.__repr__())


class InsertionTableBuilder(object):
    
    SENSE_KEY = 'sense'
    ANTISENSE_KEY = 'antisense'

    KEYS = [SENSE_KEY, ANTISENSE_KEY, 'exon']

    def __init__(self, intron_path, exon_path):
        """
        Init method for class
        :param str intron_path: Path to an bed annotation of introns
        :param str exon_path: Path to a bed annotation of exons
        """
        self._intron_bedtool = BedTool(intron_path)
        self._exon_bedtool = BedTool(exon_path)

    def intersectIntronsSense(self, input_bedtool):
        """
        Intersect bedtool object with intron table
        :param BedTool: A BedTool object
        """
        return input_bedtool.intersect(self._intron_bedtool, wa=True, wb=True, s=True)

    def intersectIntronsAntisense(self, input_bedtool):
        """
        Intersect bedtool with antisense table
        :param input_bedtool: A BedTool object
        """
        return input_bedtool.intersect(self._intron_bedtool, wa=True, wb=True, S=True)

    def intersectExons(self, input_bedtool):
        """
        Intersect with exon table
        :param input_bedtool: A BedTool object
        :
        """
        return input_bedtool.intersect(self._exon_bedtool, wa=True, wb=True)

    def toDataFrame(self, input_bedtool):
        """
        Convert to a dataframe
        :param input_bedtool: A BedTool object
        :return DataFrame: A dataframe containing insertion information.
        """
        unique_insertions = {InsertionRecord(line) for line in input_bedtool}
        
        return pd.DataFrame.from_records([x.toIterable() for x in unique_insertions], columns=InsertionRecord.FIELD_NAMES)

    def calcGsp(self, row):
        """
        Calculates the GSP score for each gene in the dataset
        :param row: A dataframe row
        :return row: A 
        """

        sense_smooth = row[self.SENSE_KEY] + 1
        antisense_smooth = row[self.ANTISENSE_KEY] + 1
        
        row['GSP'] = math.log((sense_smooth / antisense_smooth), 2) * \
                math.log((sense_smooth * antisense_smooth), 10)
        
        return row

    def buildTable(self, input_file, output_file):
        """
        Builds a table containing the insertion counts.
        """
        
        input_bedtool = BedTool(input_file)
        input_bed = input_bedtool.bam_to_bed()

        # get sense strand insertions in gene introns
        sense = self.intersectIntronsSense(input_bed)
        
        # get antisense insertions in gene introns
        antisense = self.intersectIntronsAntisense(input_bed)
        
        # get both sense and antisense insertions in gene exons
        exon = self.intersectExons(input_bed)
        
        frames = map(self.toDataFrame, [sense, antisense, exon])

        count_dfs = []
        for i, curr_frame in enumerate(frames):
            counts = curr_frame.groupby(['chrom', 'gene', 'insertion_orientation'])['start'].apply(lambda x:x.nunique()).reset_index()
            counts.set_index('gene', inplace=True)
            counts.drop(['chrom', 'insertion_orientation'], inplace=True, axis=1)
            counts.columns = [self.KEYS[i]]
            count_dfs.append(counts)
        
        joined = pd.concat(count_dfs, axis=1)
        filled = joined.fillna(value=0)
        # Adds the GSP score to the table
        gsp_filled = filled.apply(self.calcGsp, axis=1)
        
        gsp_filled.to_csv(output_file, sep='\t')


if __name__ == '__main__':

    opts= parse_cmdline_params(sys.argv[1:])

    table_builder = InsertionTableBuilder(opts.intron_bed, opts.exon_bed)
    table_builder.buildTable(opts.bed_file, opts.output_file)

