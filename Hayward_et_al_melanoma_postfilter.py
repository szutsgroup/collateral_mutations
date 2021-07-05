
# This Python script will extract some relevant parameters of close mutational pairs in the melanoma dataset,
# that will help to decide whether the close pair can be considered a collateral pair.
# 1. Mean MAPQ for both mutations in the tumor and matching normal samples
# 2. Number of supporting reads or both mutations in the tumor and matching normal samples
# 3. Number of reads supporting both mutations in the tumor sample
# 4. Fisher test's P value and odds ratio for the difference between the allele frequencies of the two mutations in the tumor sample

import pysam
import numpy as np
import scipy.stats
import sys

args = list(sys.argv)
infile = args[1]
bam_folder = args[2]
outfile = open(args[3], "w")


outfile.write("\t".join(map(str, ["ID", "Mut1", "Mut2", "Dist", "Chr",
               "Pos1", "Pos2", "Sample", "Trip1", "Trip2",  
               "Cov1_t", "Supp1_t", "Cov2_t", "Supp2_t",
               "Cov1_n", "Supp1_n", "Cov2_n", "Supp2_n",
               "Fisher.p", "Fisher.oddsr", "Spanning_read_n", 
               "MAPQ1_t", "MAPQ2_t", "MAPQ1_n", "MAPQ2_n"])) + "\n")

with open(infile, "r") as f:
    tab = []
    first_line = f.readline()
    for line in f:
        ID,Mut1,Mut2,Dist,Chr,Pos1,Pos2,Sample,Trip1,Trip2 = line.strip().split('\t')
        tumor_bam = bam_folder + "/" + Sample
        normal_bam = bam_folder + "/" + "_".join(Sample.split('_')[:2] + ["normal_remap_RMdup.bam"]) 
        print(tumor_bam + normal_bam)
        tumbamfile = pysam.AlignmentFile(tumor_bam, "rb")
        normbamfile = pysam.AlignmentFile(normal_bam, "rb")
        supp_t1 = 0
        mapq1_t = []
        reads1 = []

        if Mut1 == "NA":
            supp_t1 = "NA"
        else:
            for pileupcolumn in tumbamfile.pileup(Chr, int(Pos1), int(Pos1)+1, stepper = "samtools", min_base_quality = 0, ignore_orphans = False, ignore_overlaps = False):
                if pileupcolumn.pos == int(Pos1) - 1:
                    cov1 = pileupcolumn.nsegments
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            if pileupread.alignment.query_sequence[pileupread.query_position] == Mut1.split(">")[1]:
                                supp_t1 += 1
                                reads1.append(pileupread.alignment.query_name)
                                mapq1_t.append(pileupread.alignment.mapping_quality)
        supp_n1 = 0
        mapq1_n = []
        reads1n = []
        if Mut1 == "NA":
            supp_n1 = "NA"
        else:
            for pileupcolumn in normbamfile.pileup(Chr, int(Pos1), int(Pos1)+1, stepper = "samtools", min_base_quality = 0, ignore_orphans = False, ignore_overlaps = False):
                if pileupcolumn.pos == int(Pos1) - 1:
                    cov1n = pileupcolumn.nsegments
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            if pileupread.alignment.query_sequence[pileupread.query_position] == Mut1.split(">")[1]:
                                supp_n1 += 1
                                reads1n.append(pileupread.alignment.query_name)
                                mapq1_n.append(pileupread.alignment.mapping_quality)
        supp_t2 = 0
        mapq2_t = []
        reads2 = []
        if Mut2 == "NA":
            supp_t2 = "NA"
        else:
            for pileupcolumn in tumbamfile.pileup(Chr, int(Pos2), int(Pos2)+1, stepper = "samtools", min_base_quality = 0, ignore_orphans = False, ignore_overlaps = False):
                if pileupcolumn.pos == int(Pos2) - 1:
                    cov2 = pileupcolumn.nsegments
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            if pileupread.alignment.query_sequence[pileupread.query_position] == Mut2.split(">")[1]:
                                supp_t2 += 1
                                reads2.append(pileupread.alignment.query_name)
                                mapq2_t.append(pileupread.alignment.mapping_quality)
        supp_n2 = 0
        mapq2_n = []
        reads2n = []
        if Mut2 == "NA":
            supp_n2 = "NA"
        else:
            for pileupcolumn in normbamfile.pileup(Chr, int(Pos2), int(Pos2)+1, stepper = "samtools", min_base_quality = 0, ignore_orphans = False, ignore_overlaps = False):
                if pileupcolumn.pos == int(Pos2) - 1:
                    cov2n = pileupcolumn.nsegments
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            if pileupread.alignment.query_sequence[pileupread.query_position] == Mut2.split(">")[1]:
                                supp_n2 += 1
                                reads2n.append(pileupread.alignment.query_name)
                                mapq2_n.append(pileupread.alignment.mapping_quality)
       
        if supp_t1 != "NA" and supp_t2 != "NA":                            
            oddsratio, pvalue = scipy.stats.fisher_exact([[cov1, supp_t1], [cov2, supp_t2]])
        else:
            oddsratio, pvalue = "NA", "NA"
        reads1 = set(reads1)
        reads2 = set(reads2)
        mapq1_t = np.mean(mapq1_t)
        mapq2_t = np.mean(mapq2_t)
        mapq1_n = np.mean(mapq1_n)
        mapq2_n = np.mean(mapq2_n)
        try:
            outfile.write("\t".join(map(str, [ID,Mut1,Mut2,Dist,Chr,Pos1,Pos2,Sample,Trip1,Trip2, 
               cov1, supp_t1, cov2, supp_t2, cov1n, supp_n1, cov2n, supp_n2, 
               round(pvalue, 3), round(oddsratio, 2), 
               len(reads1.intersection(reads2)), mapq1_t, mapq2_t, mapq1_n, mapq2_n])) + "\n")
            print("\t".join(map(str, [ID,Mut1,Mut2,Dist,Chr,Pos1,Pos2,Sample,Trip1,Trip2,
               cov1, supp_t1, cov2, supp_t2, cov1n, supp_n1, cov2n, supp_n2,
               round(pvalue, 3), round(oddsratio, 2),
               len(reads1.intersection(reads2)), mapq1_t, mapq2_t, mapq1_n, mapq2_n])) + "\n")
        except TypeError:
            outfile.write("\t".join(map(str, [ID,Mut1,Mut2,Dist,Chr,Pos1,Pos2,Sample,Trip1,Trip2, 
               cov1, supp_t1, cov2, supp_t2, cov1n, supp_n1, cov2n, supp_n2, 
               pvalue, oddsratio, len(reads1.intersection(reads2)), mapq1_t, mapq2_t, mapq1_n, mapq2_n])) + "\n")
            print("\t".join(map(str, [ID,Mut1,Mut2,Dist,Chr,Pos1,Pos2,Sample,Trip1,Trip2,
               cov1, supp_t1, cov2, supp_t2, cov1n, supp_n1, cov2n, supp_n2,
               round(pvalue, 3), round(oddsratio, 2),
               len(reads1.intersection(reads2)), mapq1_t, mapq2_t, mapq1_n, mapq2_n])) + "\n")
outfile.close()
tumbamfile.close()
