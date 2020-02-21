#!/bin/bash

#SBATCH -J CMS                  # Name of the job
#SBATCH -o CMS.out              # Name of file that will have program output
#SBATCH -e CMS.err              # Name of the file that will have job errors, if any
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 1                    # Number of cores
#SBATCH -p normal               # Partition
#SBATCH --mail-user=andrew.fields@tamucc.edu
#SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         # email me when the job ends
#SBATCH --time=96:00:00

module load cms

echo $FILE

mpirun -np 1 cms $FILE
