#!/bin/bash
#SBATCH --job-name=gemma_DLPFC
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=./SLURM_OUT/%x_%A_%a.out
#SBATCH --error=./SLURM_OUT/%x_%A_%a.err

module load gemma

input_dir="/home/zliu9/zliu9/ROSMAP2025/proteimics_input_for_GEMMA"

out_dir="/home/zliu9/zliu9/ROSMAP2025/GEMMA_Data_out"
quant_trait="Spinal_cord_protein"

cd ${input_dir}

## Generate Kinship matrix
gemma -g ${quant_trait}_dt.tsv.gz -p ${quant_trait}.cogn_global.tsv -gk 2 -notsnp -o ${quant_trait}.cov_mat

## Run LMM
for pheno in "gait_speed"  "motor_dexterity"  "motor_gait"  "motor_handstreng"  "motor10"  "bradysc"  "gaitsc" "parksc"  "rigidsc"  "tremsc"  "LewyBody"  "gpath"  "amyloid"  "tangles"  "cogn_global"  "ADD"  "caa" ; do

	echo Analyze $pheno vs. $quant_trait

	gemma -g ${quant_trait}_dt.tsv.gz -a ${quant_trait}_anno.tsv -k ./output/${quant_trait}.cov_mat.sXX.txt -c ${quant_trait}_cov_bim.tsv -p ${quant_trait}.${pheno}.tsv -lmm 4 -notsnp -o ${quant_trait}.${pheno}.lmm

	cut -f1-3,8-15 ./output/${quant_trait}.${pheno}.lmm.assoc.txt > ${out_dir}/${quant_trait}.${pheno}.lmm.assoc.txt
done

exit


