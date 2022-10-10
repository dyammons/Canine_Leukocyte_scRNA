#!/usr/bin/env bash

### This script will generate a reference to be used for cellranger counts ###

### NOTE: ###
# --genome is output directory (so change if you want a different name)
# --fasta is the path to whole genome fasta file
# --genes is a link to the .gtf (use primary assembly from ensembel - dont use all assemblys)
# --memgb specifices RAM; optional, but can be useful for optimization of each genome

export PATH=/projects/$USER/software/cellranger-6.1.2:$PATH ### Ensure path is correct ###

#ensure all paths/file names are correct!
cellranger mkref --genome=canine_ref_genome \
                 --fasta=/projects/$USER/references/wholeK9genome.fasta \
                 --genes=/projects/$USER/references/CanFam_filtered.gtf
