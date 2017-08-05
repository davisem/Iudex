#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""

import pandas as pd
import unittest
import tempfile

from pybedtools import BedTool
from Iudex.containers.python.intron_exon_counts import InsertionTableBuilder, InsertionRecord

INPUT_BED=\
"""chr1	101	102	BLAH	0	+
chr1	201	202	AK	0	-
"""

INTRON_BED=\
"""chr1	100	200	FOO	0	+
chr1	200	300	BAR	0	+
"""

EXON_BED=\
"""chr1	200	300	FOO	0	+
"""


DUPLICATE_BED=\
"""chr1	100	300	BLAH	0	+	chr1	0	500	GENE1	0	+
chr1	200	300	BLAH	0	+	chr1	0	500	GENE1	0	+
chr1	200	400	BLAH	0	-	chr1	0	500	GENE1	0	+
chr1	300	400	BLAH	0	-	chr1	0	500	GENE1	0	+
"""


class InsertionTableBuilderTests(unittest.TestCase):

	def setUp(self):
		introns = tempfile.NamedTemporaryFile(delete=False)
		introns.write(INTRON_BED)
		introns.close()
		self.introns = introns.name

		exons = tempfile.NamedTemporaryFile(delete=False)
		exons.write(EXON_BED)
		exons.close()
		self.exons = exons.name

		input_bed = tempfile.NamedTemporaryFile(delete=False)
		input_bed.write(INPUT_BED)
		input_bed.close()
		self.input_bed = input_bed.name

	def test_Init(self):
		"""Test that the instance holds the correct attrs"""
		build = InsertionTableBuilder(self.introns, self.exons)
		self.assertTrue(hasattr(build, '_intron_bedtool'))
		self.assertTrue(hasattr(build, '_exon_bedtool'))

	def testIntersectIntronsSense(self):
		"""Test that we only intersect with the sense intron"""
		build = InsertionTableBuilder(self.introns, self.exons)
		input_bedtool = BedTool(self.input_bed)
		intersected = build.intersectIntronsSense(input_bedtool)
		self.assertEquals(intersected[0].name, 'BLAH')

	def testIntersectIntronsAntisense(self):
		"""Test that we only intersect with the correct antisesense intron"""
		build = InsertionTableBuilder(self.introns, self.exons)
		input_bedtool = BedTool(self.input_bed)
		intersected = build.intersectIntronsAntisense(input_bedtool)
		self.assertEquals(intersected[0].name, 'AK')

	def testIntersectExons(self):
		"""Test that we intersect with the correct exon"""
		build = InsertionTableBuilder(self.introns, self.exons)
		input_bedtool = BedTool(self.input_bed)
		intersected = build.intersectExons(input_bedtool)
		self.assertEquals(intersected[0].name, 'AK')

	def test_toDataFrame(self):
		"""Test that we correctly return a valid dataframe"""
		build = InsertionTableBuilder(self.introns, self.exons)
		input_bedtool = BedTool(self.input_bed)
		intersected = build.intersectExons(input_bedtool)
		self.assertTrue(isinstance(build.toDataFrame(intersected), pd.DataFrame))

	def test_buildTable(self):
		"""Test that we correctly make unique indices from the gene names"""
		expect_index = ['FOO', 'BAR']
		build = InsertionTableBuilder(self.introns, self.exons)
		mydf = build.buildTable(self.input_bed)
		self.assertListEqual(sorted(list(mydf.index.values)), sorted(expect_index))


class InsertionRecordTests(unittest.TestCase):

	def setUp(self):
		dups = tempfile.NamedTemporaryFile(delete=False)
		dups.write(DUPLICATE_BED)
		dups.close()
		self.dups = dups.name

	def test_getInsertionStart(self):
		"""Test that we get the correct start according to the insertion orientation"""
		expect_starts = [100, 200, 400, 400]
		bedtool = BedTool(self.dups)
		records = [InsertionRecord(x) for x in bedtool]
		self.assertListEqual(expect_starts, [x.start for x in records])

	def test_getOrientation(self):
		"""Test that orientations are correct given the gene and insertion orientations"""
		expect_orientations = ['SENSE', 'SENSE', 'ANTISENSE', 'ANTISENSE']
		bedtool = BedTool(self.dups)
		records = [InsertionRecord(x) for x in bedtool]
		self.assertListEqual(expect_orientations, [x.orientation for x in records])

	def test_toIterable(self):
		"""Test that this correctly returns an iterable of the first record"""
		expect = ('chr1', '100', 'SENSE', 'GENE1')
		bedtool = BedTool(self.dups)
		record = InsertionRecord(bedtool[0])
		record_iter = record.toIterable()
		self.assertTrue(all([item[0]==str(item[1]) for item in zip(expect, record_iter)]))

	def test_setBehavior(self):
		"""Test that we can uniqueify the records based on their start, chrom , gene and orientation"""
		bedtool = BedTool(self.dups)
		records = [InsertionRecord(x) for x in bedtool]
		unique = list(set(records))
		self.assertTrue(all([unique[i] != unique[i+1] for i in range(len(unique)-1)]))

	def test_equality(self):
		"""Test that the equlity behavior is correct"""
		bedtool = BedTool(self.dups)
		records = [InsertionRecord(x) for x in bedtool]
		self.assertTrue(records[0] == records[0])
		self.assertFalse(records[0] == records[1])
