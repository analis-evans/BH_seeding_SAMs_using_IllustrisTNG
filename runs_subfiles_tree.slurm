#!/usr/bin/bash
#SBATCH -o run%J.out
#SBATCH -e run%J.err
#SBATCH -p 
#SBATCH -n 11 #example
#SBATCH -N 
#SBATCH -t 
#SBATCH --mem=

mpiexec -n 11 --oversubscribe --mca orte_base_help_aggregate 0 python TNGsubfilestree_descendant_ID_snapshot_table_files_parallel_MPI.py $1

