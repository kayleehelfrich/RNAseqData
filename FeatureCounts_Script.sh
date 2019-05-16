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

Bam1=PREGB6-E1.bam
Bam2=PREGB6-E2.bam
Bam3=PREGB6-E3.bam
Bam4=PREGB6-E4.bam
Bam5=PREGB6-E5.bam
Bam6=PREGB6-E6.bam
Bam7=PREGB6-E7.bam
Bam8=PREGB6-E8.bam
Bam9=PREGB6-F1.bam
Bam10=PREGB6-F2.bam
Bam11=PREGB6-F4.bam
Bam12=PREGB6-F5.bam
Bam13=PREGB6-F6.bam
Bam14=PREGB6-F7.bam
Bam15=PREGB6-F8.bam
Bam16=PREGB6-F9.bam

gtf=/proj/seq/data/MM10_UCSC/Annotation/Genes/genes.gtf

# Print date of script start
date
echo "Running featureCounts"

featureCounts -p -B -a $gtf -o MatMLiver_MdvsEtOH.counts -T 8 -g gene_id $Bam1 $Bam2 $Bam3 $Bam4 $Bam5 $Bam6 $Bam7 $Bam8 $Bam9 $Bam10 $Bam11 $Bam12 $Bam13 $Bam14 $Bam15 $Bam16
# Make sure that the files are in the order you will want the final columns in the table

# Show date and that that its finished
date
echo "featureCounts finished"