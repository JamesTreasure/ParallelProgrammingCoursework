#!/bin/sh
# Account & partition (must have these)
#SBATCH --account=cm30225
#SBATCH --partition=teaching
# Name of job (optional)
#SBATCH --job-name=Test_Hello_World
# one node
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:15:00
# any normal shell stuff
pwd
# send mail to this address
#SBATCH --mail-user=jat52@bath.ac.uk
# Load dependant modules
module load openmpi/gcc
# Run the program
mpirun -np 62 ./distributedMemory
mpirun -np 63 ./distributedMemory
mpirun -np 64 ./distributedMemory