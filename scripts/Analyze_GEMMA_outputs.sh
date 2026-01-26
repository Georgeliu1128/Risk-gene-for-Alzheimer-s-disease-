#!/bin/bash
#SBATCH --mem=96G
spack load r

Rscript Analyze_GEMMA_outputs.r
