#!/bin/sh

#Author:Mohamood Adhil, Date:10/10/2017
function usage()
{
    echo "This script does the preprocessing sorted bam file for cnvcalls"
    echo ""
    echo "-h --help"
    echo "Main path were the bam files are present (Mandatory) --mainpath=$mainpath"
    echo "Path for interval file with tab seperated three column containing chromosome, start and end (Optional) --intervalfile=$intervalfile"
    echo "Comma seperated bam file names (Mandatory) --files=$files"
    echo "Output Path (Mandatory) --outpath=$outpath"
}

while [ "$1" != "" ]; do
    PARAM=`echo $1 | awk -F= '{print $1}'`
    VALUE=`echo $1 | awk -F= '{print $2}'`
    case $PARAM in
        -h | --help)
            usage
            exit
            ;;
        --mainpath)
            mainpath=$VALUE
            ;;
	--intervalfile)
	    intervalfile=$VALUE
	    ;;
        --files)
            files=$VALUE
            ;;
        --outpath)
            outpath=$VALUE
            ;;
        *)
            echo "ERROR: unknown parameter \"$PARAM\""
            usage
            exit 1
            ;;
    esac
    shift
done
mkdir -p $outpath

for i in $(echo $files | sed "s/,/ /g")
do
echo $i
j=${i%.bam}
genomeCoverageBed -d -ibam $mainpath/$j.bam | awk -F "\t" '{print $1"\t"$2"\t"$2"\t"$3}' > $outpath/$j.bdg
if [ "$intervalfile" != "" ] ; then
    bedtools intersect -a $outpath/$j.bdg -b $intervalfile > $outpath/$j\_temp.bdg
    mv $outpath/$j\_temp.bdg $outpath/$j.bdg
fi
samtools view -F 256 $mainpath/$j.bam | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c | awk -F " " '{print $1"\t"$2"\t"$1*$2}' | awk -F "\t" '{ sum1+= $1; sum3+= $3 } END { print sum3/sum1"\t"sum1 }' > $outpath/$j.cnt
samtools idxstats $mainpath/$j.bam >> $outpath/$j.cnt
done
