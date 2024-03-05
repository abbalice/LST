#####################################
#!/bin/bash

#PBS -N unit
#PBS -o logfiles
#PBS -j oe
#PBS -l walltime=72:00:00
#PBS -J 0-99
#PBS -l mem=8G
##PBS -l select=1:ncpus=1:mem=8G
##PBS -q mercalli_gpu

source /soft/centos7/anaconda/2021.05/etc/profile.d/conda.sh
conda activate my_environment

cd $PBS_O_WORKDIR

for (( iint = 0 ; iint <= NCELLY ; iint++ )); do
  python3 test_step2.py --workdir $(pwd) --event_loc EVENT_LOC --cellX ${PBS_ARRAY_INDEX} --cellY $iint   
done





