# CNVlocator
Accurate CNV (Amplification and Deletion) calling by comparing test and control samples (BAM files)

Required tools to run the CNV locator: python2.7, R, genomeCoverageBed, samtools

Required python dependencies: pybedtools, pandas, sklearn, numpy, pysal, scipy

Required R libraries: ggplot2

The CNVlocator contains four step:

1) Prepocessbam - Process bam file to generate bedgraph file contains read accumulation for every base and count file contains read count for every chromosome

2) Intervalcalls - The read accumulation for user defined window/intervals and sliding windows. For better result, use the median read length as the window and 1/5th of the median read length as the sliding window

3) cnvlocator - This will identifies CNVs by comparing the test sample vs control sample and other user defined parameters

4) cnvplot - This will help to visualize the CNVs across the genome or user defined location
