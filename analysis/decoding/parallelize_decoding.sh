#!/usr/bin/bash

# this program calls the decoding_objects.py script for all subjects listed in the freesurfing directory
# flags for decoding_objects.py are:
# --sub -s: 		subject on which the analysis should be applied to (default sub-01)
# --kernels -k: 	how many processes should be paralellized  (default 1)
# --algorithm -a: 	what decoder should be used (SVM or LDA) default (LDA)
# --cutoff -c: 		data outside the cutoff boundary (in std) are excluded (default Inf)
# --smooth: 		which smoothed data should be used (default 0)
# --feature -f: 	which feature reduction should be used

# define the directory where the raw data is stored
FSDATA=$HOME/Documents/Magic_fMRI/DATA/MRI/derivatives/freesurfer
# loop over every subject in RAWDATA and call qsub with docker_fs.sh to reconstruct the subject
for subject in ${FSDATA}/sub-*
do
	x=$(basename "$subject")
	cores=10
	echo $x
	qsub -l nodes=1:ppn=$cores \
		-l mem=6 \
		-l walltime=4:00:00:00 \
		-F "--sub $x --runs post --over tricks --perms 1000 --kernels $cores" \
		decode_magic_effects.py -q shared
done
