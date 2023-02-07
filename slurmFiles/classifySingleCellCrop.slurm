#!/bin/bash
#SBATCH -J classifySingleCellCrop		# Job name
#SBATCH -o %x-%j.out                    	# Name of stdout output file (%x is the job name and %j is the jobId)
#SBATCH -p gpu-a100                             # Queue name
#SBATCH -N 1                                    # Total number of nodes requested
#SBATCH -n 1                                    # Total number of threas tasks requested (128 per node)
#SBATCH -t 03:00:00                             # Run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=tyler_jost@utexas.edu       # Address email notifications
#SBATCH --mail-type=all                         # Email at begin and end of job
#SBATCH -A DMS21043

# Activate conda environment
source /work/08002/tjost/ls6/miniconda3/etc/profile.d/conda.sh
conda activate computerVision
cd /scratch/08002/tjost/cellMorph/notebooks/
python ./predictBBPheno.py $SLURM_JOB_ID
cd /work/08002/tjost/ls6/slurmFiles/cellMorph
scp *.out $POD:/stor/work/Brock/Tyler/cellMorph/results/TJ2201SingleCellCrop/
