import argparse
import os
import glob

def parse_vcf(vcf_file):
    total_variants = 0
    total_snps = 0
    total_non_snps = 0
    non_snp_categories = {}

    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            total_variants += 1

            data = line.strip().split('\t')
            ref, alt = data[3], data[4]

            if len(ref) == 1 and len(alt) == 1:
                total_snps += 1
            else:
                total_non_snps += 1
                variant_type = "Indel" if "INDEL" in line else "Complex"
                non_snp_categories[variant_type] = non_snp_categories.get(variant_type, 0) + 1

    return total_variants, total_snps, total_non_snps, non_snp_categories

def write_analysis(output_dir, vcf_filename, analysis_data):
    os.makedirs(output_dir, exist_ok=True)
    output_filename = os.path.join(output_dir, os.path.basename(vcf_filename).replace('.vcf', '_analysis.txt'))

    with open(output_filename, 'w') as output_file:
        for line in analysis_data:
            output_file.write(line + "\n")

def main():
    parser = argparse.ArgumentParser(description='VCF File Analysis')
    parser.add_argument('-i', '--input', required=True, help='Input directory containing VCF files')
    parser.add_argument('-o', '--output', required=True, help='Output directory for analysis files')
    args = parser.parse_args()

    for vcf_file in glob.glob(os.path.join(args.input, '*.vcf')):
        total_variants, total_snps, total_non_snps, non_snp_categories = parse_vcf(vcf_file)
        analysis_data = [
            f"Total Variants: {total_variants}",
            f"Total SNPs: {total_snps}",
            f"Total Non-SNP Variants: {total_non_snps}"
        ] + [f"{category}: {count}" for category, count in non_snp_categories.items()]

        write_analysis(args.output, vcf_file, analysis_data)

if __name__ == "__main__":
    main()
