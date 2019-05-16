#!/bin/bash

# Check for inputs
if [[ $# -eq 0 ]] ; then
    echo 'Please enter the reference path'
    exit 2
fi

# Print date of script start
date
echo "Running bowtie2"

# Assign Variables
refpath=$1

module load bowtie2
sleep 10

# Run Bowtie2

for file in $(ls -1 *R1*.fastq.gz); do
	file_input1=$file
	file_input2=$(echo $file_input1 | sed -r 's/_R1_/_R2_/g')
	file_output=$(echo $file_input1 | sed -r 's/_.*/.sam/g')
	slurm_output=$(echo $file_input1 | sed -r 's/_.*/_slurm/g')

	# Test
	echo "bowtie2 --non-deterministic --met-stderr . --sensitive -x $refpath -1 $file_input1 -2 $file_input2 -S $file_output"

	# Run Bowtie 2
	sbatch -p general -N 1 -n 8 --mem 32768 -t 4-0 -o $slurm_output.out -e $slurm_output.err --mail-type=END,FAIL --mail-user=khelfri@live.unc.edu --wrap="bowtie2 --non-deterministic --met-stderr $errorpath --sensitive -x $refpath -1 $file_input1 -2 $file_input2 -S $file_output"
	

	sleep 10
done


# Show date and that that its finished
date
echo "ending bowtie2"
