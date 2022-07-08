# srun -p brc,shared --mem=20G --pty /bin/bash
base_dir=/scratch/users/k2142172
exome_dir=/scratch/datasets/ukbiobank/exomes
out_dir=${base_dir}/outputs/ofg/top_hits_ukbb

plink2=${base_dir}/packages/plink2
bcftools=${base_dir}/packages/bcftools-1.14/bcftools

for i in `seq 1 22`;
do
  exome_bed=${exome_dir}/ukb23155_c${i}_b0_v1.bed
  exome_bim=${exome_dir}/UKBexomeOQFE_chr${i}.bim
  exome_fam=${exome_dir}/ukb23155_c20_b0_v1_s200632.fam
  var_id=${out_dir}/vars_as_ukbb_ids_${i}.txt

  $plink2 --bed $exome_bed --bim $exome_bim --fam $exome_fam --extract $var_id \
  --export vcf --out ${out_dir}/vars_from_ukbb_${i}
done
