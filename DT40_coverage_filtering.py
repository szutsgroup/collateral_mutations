# This Python script will check each mutation for signs of sequencing artefacts:
# 1. Extract mapping quality for all reads overlapping the mutated site
# 2. Mark those mutations that have sudden jumps in coverage on both sides in a +/- 200 bps window
# Sudden jump limit is max(5, median_coverage_for_window / 6)

import pysam
import matplotlib.pyplot as plt
import numpy
import sys

args = list(sys.argv)
infile = args[1]
bamfolder = args[2]
outfile = args[3]

out = open(outfile, "w")

with open(infile, "r") as f:
    for line in f:
        pos = []
        cov = []
        words=line.strip().split("\t")
        if words[0] == "X.sample_name":
            continue
        samn=words[10].split("_")[0]
        filename=bamfolder + "/" + samn + ".bam"
        bamfile = pysam.AlignmentFile(filename, "rb")
        for pileupcolumn in bamfile.pileup(str(words[1]), int(words[2])-200, int(words[2]) + 200, stepper = "samtools", min_mapping_quality = 0, flag_filter = 0, ignore_overlaps = False, ignore_orphans = False):
            cov.append(pileupcolumn.n)
            pos.append(pileupcolumn.pos)
        if int(words[2]) in pos:
            middle = pos.index(int(words[2]))
        else:
            out.write("\t".join(words + ["NA\n"]))
            continue
	# extract coverage values and coverage differences
        cov = cov[middle-100:middle+100]
        pos = pos[middle-100:middle+100]
        dcov = [i - j for i, j in zip(cov[1:] + [0], cov)]
        dcov = dcov[:-1]
        limit = max(numpy.median(cov)/6, 5)

	# get median coverage in the window around the mutation
        mqs = []
        for read in bamfile.fetch(str(words[1]), int(words[2]), int(words[2])+1):
            mqs.append(read.mapping_quality)

        if len([x for x in dcov[0:99] if x >= limit]) > 0 and len([x for x in dcov[101:] if x <= -limit]) > 0:
            out.write("\t".join(words + ["TRUE", str(numpy.average(mqs)), "\n"]))
            print("\t".join(words + ["TRUE", str(numpy.average(mqs)), "\n"]))
        else:
            out.write("\t".join(words + ["FALSE", str(numpy.average(mqs)), "\n"]))
            print("\t".join(words + ["FALSE", str(numpy.average(mqs)), "\n"]))

out.close()
