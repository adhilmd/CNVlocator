#!/bin/sh

#
function usage()
{
    echo "This script does the preprocessing of bam for cnvcalls"
    echo ""
    echo "-h --help"
    echo "Main path were the bam files are present --mainpath=$mainpath"
    echo "Comma seperated bam file names --files=$files"
    echo "Output Path --outpath=$outpath"
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
echo "genomeCoverageBed -d -ibam $mainpath/$j.bam -g /lustre/home/amohamme/data/ChIA-PET/whole-chr-romans.txt | awk -F "\t" '{print $1"\t"$2"\t"$2"\t"$3}' > $outpath/temp.bdg"
genomeCoverageBed -d -ibam $mainpath/$j.bam -g /lustre/home/amohamme/data/ChIA-PET/whole-chr-romans.txt | awk -F "\t" '{print $1"\t"$2"\t"$2"\t"$3}' > $outpath/temp.bdg
echo "/storage/bioinfo/bin/samtools view -F 256 $mainpath/$j.bam | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c | awk -F " " '{print $1"\t"$2"\t"$1*$2}' | awk -F "\t" '{ sum1+= $1; sum3+= $3 } END { print sum3/sum1"\t"sum1 }' > $outpath/$j.cnt"
/storage/bioinfo/bin/samtools view -F 256 $mainpath/$j.bam | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c | awk -F " " '{print $1"\t"$2"\t"$1*$2}' | awk -F "\t" '{ sum1+= $1; sum3+= $3 } END { print sum3/sum1"\t"sum1 }' > $outpath/$j.cnt
echo "/storage/bioinfo/bin/samtools idxstats $mainpath/$j.bam >> $outpath/$j.cnt"
/storage/bioinfo/bin/samtools idxstats $mainpath/$j.bam >> $outpath/$j.cnt
echo "/lustre/home/amohamme/My_Tools/CNVpipeline/chrsortfile.sh $outpath/temp.bdg /lustre/home/amohamme/My_Tools/CNVpipeline/chrinfo_sort.txt $outpath/$j.bdg"
/lustre/home/amohamme/My_Tools/CNVpipeline/chrsortfile.sh $outpath/temp.bdg /lustre/home/amohamme/My_Tools/CNVpipeline/chrinfo_sort.txt $outpath/$j.bdg
echo "rm $outpath/temp.bdg"
rm $outpath/temp.bdg
done
