#!/bin/bash/

for i in $(seq 0 1 7); do sbatch runs_tree.slurm $i; done



