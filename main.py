import configparser
import os
import csv
import argparse
import pysam
from src.vcf_parser import parse_vcf
from src.plotting import plot_pre_sam_histogram, plot_post_sam_histogram, plot_hq_histogram

def main():
    parser = argparse.ArgumentParser(description="Process SBC data from VCF and SAM files.")
    parser.add_argument('--config', type=str, required=True, help='Path to the config.ini file')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory path')
    args = parser.parse_args()

    config_path = args.config
    output_dir = args.output_dir

    config = configparser.ConfigParser()
    if not config.read(config_path):
        raise FileNotFoundError(f"Config file not found at {config_path}")
    
    vcf_path = config['paths']['vcf']
    bam1_path = config['paths']['bam1']
    bam2_path = config['paths']['bam2']
    reference_fasta = config['paths'].get('reference_fasta', 'hg38.fa')
    sample1_id = config['samples']['sample1']
    sample2_id = config['samples']['sample2']
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print("Parsing VCF...")
    variants = parse_vcf(vcf_path, sample1_id, sample2_id)
    print(f"Total SBCs: {len(variants)}")
    
    quality_counts_pre = {'high_quality': 0, 'moderate_quality': 0, 'unclassified': 0}
    for variant in variants:
        quality_counts_pre[variant.initial_quality] += 1
    print("\n=== Before SAM Validation ===")
    print(f"Total SBCs: {len(variants)}")
    print(f"High-quality SBCs: {quality_counts_pre['high_quality']}")
    print(f"Moderate-quality SBCs: {quality_counts_pre['moderate_quality']}")
    print(f"Unclassified SBCs: {quality_counts_pre['unclassified']}")
    
    print("\nPlotting pre-SAM validation histogram...")
    plot_pre_sam_histogram(variants, output_dir)
    
    print("\nValidating with SAM files...")
    for variant in variants:
        try:
            bam1_file = pysam.AlignmentFile(bam1_path, "rb")
            bam2_file = pysam.AlignmentFile(bam2_path, "rb")
            reads1 = bam1_file.fetch(variant.chrom, max(1, variant.position - 21), variant.position + 20)
            reads2 = bam2_file.fetch(variant.chrom, max(1, variant.position - 21), variant.position + 20)
            variant.validate_with_sam(sample1_id, reads1, sample1_id, sample2_id, output_dir, bam1_path, reference_fasta)
            variant.validate_with_sam(sample2_id, reads2, sample1_id, sample2_id, output_dir, bam2_path, reference_fasta)
            bam1_file.close()
            bam2_file.close()
        except Exception as e:
            print(f"Error validating {variant.chrom}:{variant.position} - {str(e)}")
    
    quality_counts_post = {'high_quality': 0, 'moderate_quality': 0, 'unclassified': 0}
    for variant in variants:
        quality_counts_post[variant.quality_category] += 1
    print("\n=== After SAM Validation ===")
    print(f"Total SBCs: {len(variants)}")
    print(f"High-quality SBCs: {quality_counts_post['high_quality']}")
    print(f"Moderate-quality SBCs: {quality_counts_post['moderate_quality']}")
    print(f"Unclassified SBCs: {quality_counts_post['unclassified']}")
    
    print("\nPlotting post-SAM validation histogram...")
    try:
        plot_post_sam_histogram(variants, output_dir)
    except Exception as e:
        print(f"Error plotting post-SAM histogram: {str(e)}")
    
    print("\nPlotting High-quality histogram...")
    try:
        plot_hq_histogram(variants, output_dir)
    except Exception as e:
        print(f"Error plotting HQ histogram: {str(e)}")
    
    print("\nGenerating TSV file...")
    tsv_path = os.path.join(output_dir, 'sbc_data.tsv')
    with open(tsv_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['chrom', 'position', 'ref>alt', 'GT1', 'GT2', 'initial_quality', 'quality',
                        'sample1_bases', 'sample2_bases'])
        for variant in variants:
            writer.writerow([variant.chrom, variant.position, f"{variant.ref}>{variant.alt}",
                            str(variant.sample1_gt), str(variant.sample2_gt),
                            variant.initial_quality, variant.quality_category,
                            str(variant.sample1_base_counts), str(variant.sample2_base_counts)])
    print(f"TSV file '{tsv_path}' generated in output directory.")

if __name__ == "__main__":
    main()