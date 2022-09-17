import pandas as pd
import shutil
import os
import re
import glob
import getpass
import sys 
import argparse

configfile: './clusterConfig/slurmConfig.json'

#### NOTE: that you must use python 3.8 or higher

# Set Genome and Bin Size (edit this section if you're not working with dm6):
GenomeAssembly = 'dm6'
restrictionDigestPATH = str('/proj/mckaylab/users/astutzman/hic-analyses/dm6_arima_digest.bed')
chrsizePATH = str('/proj/mckaylab/users/astutzman/hic-analyses/juicer/AS1-21_juicer-output/dm6.chrom.sizes')
########################

# Set Module Versions:

bwaVer = str('bwa/0.7.17')

#########################

# Load Required Input Files:
BWAIndexPATH = str('/proj/seq/data/' + GenomeAssembly + '_UCSC/Sequence/BWAIndex/genome.fa')
#BWAIndexPATH = str('bwa-index/genome.fa')
GenomePATH = str('/proj/seq/data/' + GenomeAssembly + '_UCSC/Sequence/WholeGenomeFasta/genome.fa')
sampleSheetPath = str('master-samplesheet-hic.csv')
sampleDF = pd.read_csv(sampleSheetPath, comment = '#')

techList = list(set(sampleDF.techName))
fastqList= list(set(sampleDF.fastqName))
readNumList = list(sampleDF.readNum)


########################
localrules: all

rule all:
	input:
		expand("Fastq/{fastq}.fastq.gz", fastq=fastqList),
		expand("Fastq/{tech}_R1.fastq.gz", tech=techList),
		expand("Fastq/{tech}_R2.fastq.gz", tech=techList),
		directory(expand('fanc-{tech}', tech=techList)),
		expand("fanc-{tech}/hic/binned/{tech}_10kb.hic", tech=techList),
		expand("Compartments/{tech}_AB_compartmental_domains.bed", tech=techList),
		expand("Compartments/{tech}_AB_matrix.ab", tech=techList)
		

rule copyFiles:
	input:
		lambda x: list(sampleDF.htsfFile)
	output:
		expand('Fastq/{fastq}.fastq.gz', fastq = list(sampleDF.fastqName))
	message: "Copying files to Fastq directory with corrected file name"
	run:
		for htsf in list(sampleDF.htsfFile):
			outFileFilt = sampleDF [ sampleDF.htsfFile == htsf ] 
                        outFileBase = list(outFileFilt.fastqName)[0]
                        outFile = 'Fastq/{fastq}.fastq.gz'.format(fastq = outFileBase)
                        #print('THERE')
                        #print(outFile)
                        shutil.copyfile(htsf, outFile)
                        print('copied file')

rule fancAuto: 
	input:
		index = BWAIndexPATH,
		restriction = restrictionDigestPATH,
		fastq_R1 = "Fastq/{tech}_R1.fastq.gz",
		fastq_R2 = "Fastq/{tech}_R2.fastq.gz"
	output:
		outDir = directory('fanc-{tech}')
	params:
		moduleVer = bwaVer
	message: "Running Fan-C"
	shell:
		"""
		module load {params.moduleVer}
		fanc auto {input.fastq_R1} {input.fastq_R2} {output.outDir} \
		-i {input.index} \
		-n {wildcards.tech} \
		-b 5Mb, 1Mb, 500kb, 250kb, 100kb, 50kb, 25kb, 10kb, 5kb, 2kb, 1kb \
		-g {input.restriction} \
		"""

rule fancCompartments:
	input:
		hicMat = "fanc-{tech}/hic/binned/{tech}_10kb.hic",
		genome = GenomePATH
	output:
		domains = "Compartments/{tech}_AB_compartmental_domains.bed",
		abMatrixFile = "Compartments/{tech}_AB_matrix.ab",
		enrichmentPlot = "Compartments/{tech}_AB_enrichmentProfile.png"
	message: "Calling compartments from 10kb matrices"
	shell:
		"""
		fanc compartments -d {output.domains} -e {output.enrichmentPlot} -g {input.genome} --eigenvector-index 2 
		"""
