#!/bin/bash

# Usage explanation function
usage() {
    echo "Usage: $0 -i <input_file> -s <sample_file> -l <liftover_chain> -o <output_dir>"
    echo ""
    echo "Arguments:"
    echo "  -i <input_file>       Input file with candidate SNPs"
    echo "  -s <sample_file>      File containing list of sample IDs"
    echo "  -l <liftover_chain>   Liftover chain file for coordinate conversion (optional)"
    echo "  -o <output_dir>       Output directory"
    echo "  -b <batch_jobs>       Flag will submit as SLURM batch jobs (optional). Specify the time for each job"
    exit 1
}

# Argument parsing
batch_time=false
while getopts ":i:s:v:l:o:b" opt; do
    case ${opt} in
        i)
            input_file=$OPTARG
            ;;
        s)
            sample_file=$OPTARG
            ;;
        l)
            liftover_chain=$OPTARG
            ;;
        o)
            output_dir=$OPTARG
            ;;
        b)
            batch_time=true
            ;;
        \?)
            echo "Invalid option: $OPTARG" 1>&2
            usage
            ;;
        :)
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done

##for testing
#liftover_chain=/data/abattle4/aomdahl1/reference_data/liftOver_chains/hg19ToHg38.over.chain.gz
#input_file=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.250kb.0.2r2.prune.in
#output_dir=test
#sample_file=/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/selection_stats/Fst/unrelated_samples.txt
#bash /data/abattle4/aomdahl1/snp_utils/Fst_from_1KG.sh -s EAS_samples.txt -i /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.250kb.0.2r2.prune.in  -l /data/abattle4/aomdahl1/reference_data/liftOver_chains/hg19ToHg38.over.chain.gz -o EAS_test


# Check if mandatory arguments are provided
if [ -z "${input_file}" ] || [ -z "${sample_file}" ]  || [ -z "${output_dir}" ]; then
    usage
fi


#Check- if 
if [ ${sample_file} = "UNRELATED" ]; then
    echo "Using the default list of unrelated samples"
   sample_file=/data/abattle4/lab_data/1000Genomes/unrelated_sample_ids.initial.txt
fi

#User report
nsamps=`wc -l ${sample_file}| cut -f 1 -d " " `
echo "Going to only extract the ${nsamps} samples given"

# Create output directory if it doesn't exist
mkdir -p "${output_dir}"


#Build the bed file for the query:
ml bedtools
awk -F ":" '(NF>1){print "chr"$1"\t"$2-1"\t"$2"\t"$0}' ${input_file} > ${output_dir}/bed_file_query.bed
bedtools sort -i ${output_dir}/bed_file_query.bed > ${output_dir}/t && mv ${output_dir}/t ${output_dir}/bed_file_query.bed
snp_query_file="${output_dir}/bed_file_query.bed"


# If liftover is specified, perform liftover
if [ ! -z "${liftover_chain}" ]; then
    #Build the bed file
    echo "Automatically lifting over to specificed chain"
    ~/.bin/liftOver "${output_dir}/bed_file_query.bed" "${liftover_chain}" "${output_dir}/lifted.over.variants.bed" "${output_dir}/unlifted_variants.bed"
    snp_query_file="${output_dir}/lifted.over.variants.bed"
fi

# Split SNPs by chromosome, and make sure teh files are sorted
for i in {1..22}; do
    grep -e "chr${i}\s" ${snp_query_file} > ${output_dir}/chr${i}_snps.tmp
    bedtools sort -i ${output_dir}/chr${i}_snps.tmp > ${output_dir}/chr${i}_snps.bed
    rm ${output_dir}/chr${i}_snps.tmp
done

# Loop through each chromosome file and extract SNPs for the specified samples
#Accept an option to push this to slurm?
ml bcftools
for chr_file in ${output_dir}/chr*.bed; do
    chr=$(basename "${chr_file}" .snps)
    output_vcf="${output_dir}/${chr}.vcf.gz"
    # Extract SNPs for specified samples using bcftools
    chr_assign=`echo ${chr_file} | cut -f 2 -d "/" | cut -f 1 -d "_"`
    vcf_file=/data/abattle4/lab_data/1000Genomes/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_${chr_assign}.recalibrated_variants.vcf.gz
    echo "Extracting SNPs for ${chr_file}"
    echo """ bcftools view -S "${sample_file}" -R "${chr_file}" -o "${output_vcf}" -O z "${vcf_file}" """
    if [ "$batch_time" = false ]; then
        bcftools view -S "${sample_file}" -R "${chr_file}" -o "${output_vcf}" -O z "${vcf_file}"
    fi
    echo "Finished SNPs for ${chr_file}"
done

echo "Processing completed. Output files are in ${output_dir}"
