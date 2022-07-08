# srun -p brc,shared --mem=20G --pty /bin/bash
base_dir=/scratch/users/k2142172
exome_dir=/scratch/datasets/ukbiobank/exomes
out_dir=${base_dir}/outputs/ofg/top_hits_ukbb

plink2=${base_dir}/packages/plink2
bcftools=${base_dir}/packages/bcftools-1.14/bcftools

exome_bed=${exome_dir}/ukb23155_c11_b0_v1.bed
exome_bim=${exome_dir}/UKBexomeOQFE_chr11.bim
exome_fam=${exome_dir}/ukb23155_c20_b0_v1_s200632.fam

#var_id=${base_dir}/tests/vars_as_ukbb_ids_11.txt

for i in `seq 1 22`;
do
$plink2 --bed $exome_bed --bim $exome_bim --fam $exome_fam --extract ${out_dir}/vars_as_ukbb_ids_${i}.txt \
--export vcf --out ${out_dir}/vars_from_ukbb_${i}
done

