# set variable paths
base_dir=/scratch/users/k2142172
out_dir=${base_dir}/outputs/ofg
bgzip=${base_dir}/packages/anaconda3/envs/peddy/bin/bgzip
bcftools=${base_dir}/packages/bcftools-1.14/bcftools

#/mnt/lustre/groups/prescott_lab/OFG_Sanger_Samples/Sanger_SFTP
#/mnt/lustre/groups/prescott_lab/OFG_Sanger_Samples/Sanger_SFTP/ofg_samples.extended_set.vcf.gz

# split vcf by chromosome number
# split multiallelic entries into their own separate rows
# only keep polymorphic variants
# index produced vcf
for i in `seq 1 22`;
do
  $bcftools view ${out_dir}/ofg_samples.extended_set.vcf.gz --regions chr${i} -O u | \
  $bcftools norm --multiallelics - -O u | $bcftools view -c1.0:minor -o ${out_dir}/ofg_chr${i}.vcf.gz
  $bcftools index -c ${out_dir}/ofg_chr${i}.vcf.gz
done

# annotate all chr vcfs with prefix provided as arg on command line
sbatch ${base_dir}/scripts/misc_projects/vep_chr.sh /scratch/users/k2142172/outputs/ofg/ofg_chr

# extract the CSQ and header sections from a vcf
$bcftools view -h  ${out_dir}/ofg_chr22.vcf.gz.batch.vep.vcf | grep CSQ > ${out_dir}/ofg_vep_anno_csq.txt
$bcftools view -h ${out_dir}/ofg_chr22.vcf.gz.batch.vep.vcf | tail -n 1 > ${out_dir}/ofg_vep_anno_headers.txt

# create tsv file convert from vcf for each chr with headers and csq as additional args
for i in `seq 1 22`;
do
  sbatch ${base_dir}/scripts/misc_scripts/vcf_to_tsv.sh ${out_dir}/ofg_chr${i}.vcf.gz.batch.vep.vcf ${out_dir}/ofg_vep_anno_headers.txt ${out_dir}/ofg_vep_anno_csq.txt
done

# merge chr data into RData file
sbatch ${base_dir}/scripts/ofg_wes/run_ofg_merge.sh

# filter chr data for variants of interest
sbatch ${base_dir}/scripts/ofg_wes/run_ofg_filter.sh
