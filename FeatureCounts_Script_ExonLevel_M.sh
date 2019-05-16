#!/bin/bash
#
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 32768
#SBATCH -t 1-0 # time (D-HH:MM)
#SBATCH -o featureCounts.%N.%j.out
#SBATCH -e featureCounts.%N.%j.err
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=khelfri@live.unc.edu

module load subread
sleep 8 

Bam1=PREGB6-A1.bam
Bam2=PREGB6-A2.bam
Bam3=PREGB6-A3.bam
Bam4=PREGB6-A4.bam
Bam5=PREGB6-A5.bam
Bam6=PREGB6-A6.bam
Bam7=PREGB6-A7.bam
Bam8=PREGB6-A8.bam
Bam9=PREGB6-B1.bam
Bam10=PREGB6-B2.bam
Bam11=PREGB6-B4.bam
Bam12=PREGB6-B5.bam
Bam13=PREGB6-B6.bam
Bam14=PREGB6-B7.bam
Bam15=PREGB6-B8.bam
Bam16=PREGB6-B9.bam

gtf=/proj/seq/data/MM10_UCSC/Annotation/Genes/genes.gtf

# Print date of script start
date
echo "Running featureCounts"

featureCounts -p -B -f -s 0 -t exon -a $gtf -o MatMLiver_MDvsEtOH_EXON.counts -T 8 -g gene_id $Bam1 $Bam2 $Bam3 $Bam4 $Bam5 $Bam6 $Bam7 $Bam8 $Bam9 $Bam10 $Bam11 $Bam12 $Bam13 $Bam14 $Bam15 $Bam16
# Make sure that the files are in the order you will want the final columns in the table

# Show date and that that its finished
date
echo "featureCounts finished"