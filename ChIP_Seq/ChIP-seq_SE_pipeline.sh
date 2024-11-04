#!/bin/sh

#SBATCH --job-name=fqToWig
#SBATCH --output=%x-%j_"$sample_name".out      # Output to a file with a name constructed from the job name (%x)
                                               # and the job ID number (%j) (Default: slurm-%j.out)
#SBATCH --time=48:00:00
#SBATCH --mem=20G
#SBATCH  -N 1
#SBATCH  -n 23
#SBATCH --mail-type=END,FAIL

dataset1=$1
bwaIndex=$2
chromSizes=$3
prefixLess=$4
outDir=$5

module add java/1.8.0_121 trimmomatic/0.36 fastqc/0.11.5 bwa/0.7.15 samtools/1.3.1 bedtools/2.26.0 scipy-stack/2018b python/3.5.4 perl/5.22.4 picard/2.18.9
module del mugqic/python/2.7.14

datasetLess1=`echo $dataset1 | tr '/' ' ' |awk '{print $NF}'|sed 's#.fastq.gz##'`
prefix=$outDir"/"$prefixLess
prefix_1=$outDir"/"$prefixLess"_1"

mkdir -p $outDir
cd $outDir
#trim et qc
echo "###########################
starting fastqc
###########################"
fastqc $dataset1 --outDir $outDir/
rename $datasetLess1 ${prefixLess}_1 $outDir/$datasetLess1*
echo "###########################
starting trim
###########################"
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar SE \
-threads 20 \
$dataset1 \
${prefix_1}_trimmed.fastq \
ILLUMINACLIP:/nfs3_ib/ip29/ip29/jacques_group/Marc-Antoine/pipeline_ChIPseq/adapters_all_0.32.fa:2:30:15 \
LEADING:30 \
TRAILING:30 \
MINLEN:23

echo "###########################
starting fastqc
###########################"
fastqc ${prefix_1}_trimmed.fastq --outDir $outDir/

#alignement
echo "###########################
starting bwa
###########################"
bwa mem -t 20 -split $bwaIndex ${prefix_1}_trimmed.fastq | samtools view -hb - > $prefix.bam
samtools view -hb  -F 256 -F 2048 $prefix.bam | samtools sort -T ${prefix}_temp_bam -o ${prefix}_F.bam

echo "###########################
starting aligment stat
###########################"
samtools flagstat $prefix.bam > ${prefix}_flagstat.txt
samtools flagstat ${prefix}_F.bam > ${prefix}_F_flagstat.txt


bam to bedgraph
echo "###########################
starting genomecov
###########################"
bedtools genomecov -bga -trackline -trackopts name=$prefixLess.bg -ibam ${prefix}_F.bam -g $chromSizes > $prefix.bg
head -n1 $prefix.bg > $prefix.bg.header
tail -n+2 $prefix.bg|sort -k1,1 -k2,2n |cat $prefix.bg.header - > $prefix.temp
mv $prefix.temp $prefix.bg
rm $prefix.bg.header

#bin bg 10pb
echo "###########################
starting scl
###########################"
perl /nfs3_ib/ip29/ip29/jacques_group/Marc-Antoine/Bachand_F/normalisationCE_2.0/bedgraph_to_wig.pl --bedgraph $prefix.bg --wig $prefix.wig --step 10


#scale = scl (normalize) files
python3 /nfs3_ib/ip29/ip29/jacques_group/Marc-Antoine/Bachand_F/normalisationCE_2.0/sclOnTotal.py $prefix.wig $prefix.scl.wig
python3 /nfs3_ib/ip29/ip29/jacques_group/Marc-Antoine/Bachand_F/normalisationCE_2.0/wigToBg.py $prefix.wig $prefix.bg
python3 /nfs3_ib/ip29/ip29/jacques_group/Marc-Antoine/Bachand_F/normalisationCE_2.0/wigToBg.py $prefix.scl.wig $prefix.scl.bg

/nfs3_ib/ip29/ip29/jacques_group/Marc-Antoine/Bachand_F/normalisationCE_2.0/wigToBigWig $prefix.wig $chromSizes $prefix.bw
/nfs3_ib/ip29/ip29/jacques_group/Marc-Antoine/Bachand_F/normalisationCE_2.0/wigToBigWig $prefix.scl.wig $chromSizes $prefix.scl.bw

echo "###########################
starting compressing
###########################"
gzip ${prefix_1}_trimmed.fastq
