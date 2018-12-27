# CNVlocator
Accurate CNV (Amplification and Deletion) calling by comparing test and control samples (BAM files)

Required tools to run the CNV locator: python2.7, R, genomeCoverageBed, samtools

Required python dependencies: pybedtools, pandas, sklearn, numpy, pysal, scipy

Required R libraries: ggplot2

This CNVlocator contains four step:

1) Prepocessbam - process bam file to generate bedgraph file contains read accumulation for every base and count file contains read count for every chromosome along with mean read length and total reads

2) Intervalcalls - The read accumulation for user defined intervals and sliding windows. For better result the 
