base_dir=/scratch/users/k2142172
out_dir=${base_dir}/tests/vep
bgzip=${base_dir}/packages/anaconda3/envs/peddy/bin/bgzip
bcftools=${base_dir}/packages/bcftools-1.14/bcftools

#/mnt/lustre/groups/prescott_lab/OFG_Sanger_Samples/Sanger_SFTP
#/mnt/lustre/groups/prescott_lab/OFG_Sanger_Samples/Sanger_SFTP/ofg_samples.extended_set.vcf.gz

for i in `seq 11 22`;
do
$bcftools view ofg_samples.extended_set.vcf.gz --regions chr${i} -O z -o ofg_chr${i}.vcf.gz
$bcftools index -c ofg_chr${i}.vcf.gz
done

$bcftools view -h  ofg_chr22.vcf.gz.batch.vep.vcf | grep CSQ > ofg_vep_anno_csq.txt
$bcftools view -h ofg_chr22.vcf.gz.batch.vep.vcf | tail -n 1 > ofg_vep_anno_headers.txt

for i in `seq 1 22`;
do
  sbatch vcf_to_tsv.sh ofg_chr${i}.vcf.gz.batch.vep.vcf ofg_vep_anno_headers.txt ofg_vep_anno_csq.txt
done

conda activate r4
