#!/bin/bash
#
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 32768
#SBATCH -t 2-0 # time (D-HH:MM)
#SBATCH -o bam_sort.%N.%j.out
#SBATCH -e bam_sort.%N.%j.err
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=khelfri@live.unc.edu

# Print date of script start
date
echo "Running samtools"

module load samtools
sleep 10

for file in *.sam
do
	echo $file
	describer=$(echo ${file} | sed 's/.sam//')
	echo $describer

	# Run samtools
	samtools view -b $file > ${describer}.uns.bam 
	sleep 2
	samtools sort ${describer}.uns.bam -o ${describer}.bam
	sleep 2
	samtools index ${describer}.bam
	sleep 2
done

# Show date and that that its finished
date
echo "samtools finished"