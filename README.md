# CNVlocator

Author: Mohamood Adhil and Ramveer Choudhary

Accurate CNV (Amplification and Deletion) calling by comparing test and control samples (BAM files)

Required tools to run the CNV locator: python2.7, R, genomeCoverageBed, samtools

Required python dependencies: pybedtools, pandas, sklearn, numpy, pysal, scipy

Required R libraries: ggplot2

Input files required to run this tool: 

1) sorted bam test and control files & 2) (Optional) Interval file containing three columns chromosome, start and end

Output files:

1) Text files containing most probable CNV locations & 2) CNV plot 

The CNVlocator contains four step:

1) Prepocessbam - Process bam file to generate bedgraph file contains read accumulation for every base and count file contains read count for every chromosome. Input should be sorted bam file.

2) Intervalcalls - The read accumulation for user defined window/intervals and sliding windows. For better result, use the median read length as the window size and 1/5th of the median read length as the sliding window size

3) cnvlocator - This will identifies CNVs by comparing the test sample vs control sample and other user defined parameters

4) cnvplot - This will help to visualize the CNVs across the genome or user defined location

![test_cnvplot](https://user-images.githubusercontent.com/18418058/51907337-c2c99b80-23c6-11e9-936e-78b3451a7c81.jpeg)

