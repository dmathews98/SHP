#!/bin/bash -l
# Batch script to run an MPI parallel job with the upgraded software
# stack under SGE with Intel MPI.

# 1. Force bash as the executing shell.
#$ -S /bin/bash

# 2. Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=0:10:00

# Budget setting
#$ -P Gold
#$ -A UKCP_ED_P

# 3. Request 1 gigabyte of RAM per process.
#$ -l mem=1G

# 4. Request 15 gigabyte of TMPDIR space per node (default is 10 GB)
#$ -l tmpfs=15G

# 5. Set the name of the job.
#$ -N vasptest1

# 6. Select the MPI parallel environment and 16 processes.
#$ -pe mpi 24

# 7. Set the working directory to somewhere in your scratch space.  This is
# a necessary step with the upgraded software stack as compute nodes cannot
# write to $HOME.
## Replace "<your_UCL_id>" with your UCL user ID :
##$ -wd /home/mmm0587/Scratch/output
# Set working directory to current directory
#$ -cwd

# 8. Run our MPI job.  GERun is a wrapper that launches MPI jobs on our clusters.
VASPBIN=$HOME/VASPinitial/vasp_std

gerun $VASPBIN > out-vasp
