#!/bin/bash

#SBATCH -J CMS                  # Name of the job
#SBATCH -o CMS.out              # Name of file that will have program output
#SBATCH -e CMS.err              # Name of the file that will have job errors, if any
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 20                   # Number of cores
#SBATCH -p normal               # Partition
#SBATCH --mail-user=andrew.fields@tamucc.edu
#SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         #email me when the job ends
#SBATCH --time=96:00:00

module load cms

mpirun -np 1 cms group01 2>1 | tee group01.log &
mpirun -np 1 cms group02 2>1 | tee group02.log &
mpirun -np 1 cms group03 2>1 | tee group03.log &
mpirun -np 1 cms group04 2>1 | tee group04.log &
mpirun -np 1 cms group05 2>1 | tee group05.log &
mpirun -np 1 cms group06 2>1 | tee group06.log &
mpirun -np 1 cms group07 2>1 | tee group07.log &
mpirun -np 1 cms group08 2>1 | tee group08.log &
mpirun -np 1 cms group09 2>1 | tee group09.log &
mpirun -np 1 cms group10 2>1 | tee group10.log &
mpirun -np 1 cms group11 2>1 | tee group11.log &
mpirun -np 1 cms group12 2>1 | tee group12.log &
mpirun -np 1 cms group13 2>1 | tee group13.log &
mpirun -np 1 cms group14 2>1 | tee group14.log &
mpirun -np 1 cms group15 2>1 | tee group15.log &
mpirun -np 1 cms group16 2>1 | tee group16.log &
mpirun -np 1 cms group17 2>1 | tee group17.log &
mpirun -np 1 cms group18 2>1 | tee group18.log &
mpirun -np 1 cms group19 2>1 | tee group19.log

mkdir logfiles
mv *.log logfiles/
