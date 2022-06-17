#!/bin/bash -l

#SBATCH --partition=brc,shared
#SBATCH --job-name=vep_chr
#SBATCH --time=06:00:00
#SBATCH --mem=48G
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/vep/vep_%A_%a.out
#SBATCH --array=[1-22]%6

#############################################################
# set variables

base_dir=/scratch/users/k2142172
out_dir=${base_dir}/tests/vep
bgzip=${base_dir}/packages/anaconda3/envs/peddy/bin/bgzip
bcftools=${base_dir}/packages/bcftools-1.14/bcftools

echo $SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID

# conda vep environment
conda activate vep

# remove perl path due to version conflict
unset PERL5LIB

ASSEMBLY=GRCh${2:-38}
OUT_FORMAT=${3:-vcf}
VEP_DIR="/scratch/users/k2142172/resources/${ASSEMBLY}/vep"
PLUGIN_DIR="/scratch/users/k2142172/resources/vep_plugins"

#############################################################

#Shortcut flag to switch on all of the following:
# --sift b, --polyphen b, --ccds, --uniprot, --hgvs, --symbol,
#--numbers, --domains, --regulatory, --canonical, --protein,
#--biotype, --uniprot, --tsl, --appris, --gene_phenotype
#--af, --af_1kg, --af_esp, --af_gnomad, --max_af, --pubmed,
#--variant_class

##############################################################

run_vep(){
VARIANT_FILE=$1
vep -i $VARIANT_FILE \
--everything \
--per_gene \
--pick_order rank \
--cache \
--dir_cache $VEP_DIR \
--cache_version 105 \
--output_file ${VARIANT_FILE}.batch.vep.${OUT_FORMAT} \
--assembly ${ASSEMBLY} \
--no_stats \
--dir_plugin $PLUGIN_DIR \
--plugin SpliceRegion \
--plugin MaxEntScan,${PLUGIN_DIR}/maxentscan/fordownload \
--plugin CADD,${VEP_DIR}/whole_genome_SNVs.tsv.gz \
--plugin Conservation,/scratch/users/k2142172/resources/GRCh38/vep/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
--fasta /scratch/users/k2142172/resources/${ASSEMBLY}/Homo_sapiens.${ASSEMBLY}.dna.primary_assembly.fa \
--force_overwrite \
--${OUT_FORMAT}
}

run_vep ${out_dir}/ofg_chr${i}.vcf.gz
