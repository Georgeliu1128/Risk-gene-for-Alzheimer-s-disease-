#!/bin/bash
#SBATCH --job-name=preranked
#SBATCH --output=logs/gsea_preranked_%A_%a.log
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --array=1

files_dir="/home/zliu9/zliu9/ROSMAP2025/gsea_analysis_01082026/prernk_files"
out_dir="/home/zliu9/zliu9/ROSMAP2025/gsea_analysis_01082026/gsea_output/fourth_run"

# Sample settings
high_sig=(
  "DLPFC.gpath" "DLPFC.amyloid" "DLPFC_protein.gpath"
  "DLPFC_protein.amyloid" "DLPFC_protein.tangles"
  "DLPFC_protein.caa" "Vh_protein.LewyBody" "Vh.tremsc" "Quad.tremsc" "DLPFC.cogn_global"
  "DLPFC_protein.cogn_global" "DLPFC.ADD"
  "DLPFC_protein.tremsc" "Sma.gpath" "Sma.tangles" "Sma.cogn_global" "Sma.motor_dexterity"
  "Sma.ADD" "Vh.gpath" "Vh.tangles" "Vh.cogn_global" "Quad.gait_speed"
  "Quad.ADD" "DLPFC.tangles" "DLPFC.parksc" "DLPFC.motor_dexterity"
  "DLPFC_protein.motor_gait" "DLPFC_protein.motor10"
  "Quad_protein.gait_speed" "Quad_protein.motor_handstreng"
  "Vh_protein.cogn_global"
)


###############################################
# Analyze HIGH significant gene lists
###############################################
for hs in "${high_sig[@]}"; do
    
    echo "[`date`] Processing high_sig: $hs"
    

    out_subdir_c2="${out_dir}/C2/${hs}"
    out_subdir_hallmark="${out_dir}/hallmark/${hs}"
    mkdir -p "$out_subdir_hallmark"
    mkdir -p "$out_subdir_c2"

    ranked_file="${files_dir}/${hs}.rnk"
    echo "Ranked file: $ranked_file"

    if [[ ! -f "$ranked_file" ]]; then
        echo "ERROR: Ranked file not found: $ranked_file"
        continue
    fi

    #analysis with C2 databse
    /nfs/yangfss2/data/commons/sshe227/Software/GSEA_Linux_4.4.0/gsea-cli.sh GseaPreranked \
        -rnk  "$ranked_file" \
        -gmx  /nfs/yangfss2/data/commons/sshe227/bcBiomarker/gsea/gmt/c2.all.v2025.1.Hs.symbols.gmt \
        -collapse false \
        -nperm 1000 \
        -rnd_seed 12345 \
        -out  "$out_subdir_c2"
    
    #analysis with hallmark database
    /nfs/yangfss2/data/commons/sshe227/Software/GSEA_Linux_4.4.0/gsea-cli.sh GseaPreranked \
        -rnk  "$ranked_file" \
        -gmx  /nfs/yangfss2/data/commons/sshe227/bcBiomarker/gsea/gmt/h.all.v2025.1.Hs.symbols.gmt \
        -collapse false \
        -nperm 1000 \
        -rnd_seed 12345 \
        -out  "$out_subdir_hallmark"


    echo "Completed: $hs "
done


echo "All analyses completed"