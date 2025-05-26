import configparser
import os
import csv
import argparse
import pysam
from src.vcf_parser import parse_vcf
from src.sam_validator import validate_sam
from src.plotting import plot_pre_sam_histogram, plot_post_sam_histogram, plot_hq_histogram

def main():
    # argparse로 인자 처리
    parser = argparse.ArgumentParser(description="Process SBC data from VCF and SAM files.")
    parser.add_argument('--config', type=str, required=True, help='Path to the config.ini file')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory path')
    args = parser.parse_args()

    config_path = args.config
    output_dir = args.output_dir

    # 설정 파일 읽기
    config = configparser.ConfigParser()
    if not config.read(config_path):
        raise FileNotFoundError(f"Config file not found at {config_path}")
    
    vcf_path = config['paths']['vcf']
    bam1_path = config['paths']['bam1']
    bam2_path = config['paths']['bam2']
    reference_fasta = config['paths'].get('reference_fasta', 'hg38.fa')  # 추가: 참조 서열 파일 경로
    sample1_id = config['samples']['sample1']
    sample2_id = config['samples']['sample2']
    
    # 출력 디렉토리 확인 및 생성
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # VCF 파싱 및 변이 추출
    print("Parsing VCF...")
    variants = parse_vcf(vcf_path, sample1_id, sample2_id)
    print(f"Total SBCs: {len(variants)}")
    
    # SAM 검증 전 분류
    quality_counts_pre = {'high_quality': 0, 'low_quality': 0, 'incorrect': 0}
    for variant in variants:
        quality_counts_pre[variant.quality_category] += 1
    print("\n=== Before SAM Validation ===")
    print(f"Total SBCs: {len(variants)}")
    print(f"High-quality SBCs: {quality_counts_pre['high_quality']}")
    print(f"Low-quality SBCs: {quality_counts_pre['low_quality']}")
    print(f"Incorrect SBCs: {quality_counts_pre['incorrect']}")
    
    # SAM 검증 전 히스토그램
    print("\nPlotting pre-SAM validation histogram...")
    plot_pre_sam_histogram(variants, output_dir)
    
    # SAM 검증
    print("\nValidating with SAM files...")
    variants = validate_sam(variants, bam1_path, bam2_path, sample1_id, sample2_id)
    
    # SAM 검증 후 분류
    quality_counts_post = {'high_quality': 0, 'low_quality': 0, 'incorrect': 0}
    for variant in variants:
        variant.quality_category = variant.classify_quality(use_sam=True)
        quality_counts_post[variant.quality_category] += 1
    print("\n=== After SAM Validation ===")
    print(f"Total SBCs: {len(variants)}")
    print(f"High-quality SBCs: {quality_counts_post['high_quality']}")
    print(f"Low-quality SBCs: {quality_counts_post['low_quality']}")
    print(f"Incorrect SBCs: {quality_counts_post['incorrect']}")
    
    # SAM 검증 후 히스토그램
    print("\nPlotting post-SAM validation histogram...")
    plot_post_sam_histogram(variants, output_dir)
    
    # High-quality만의 히스토그램
    print("\nPlotting High-quality histogram...")
    plot_hq_histogram(variants, output_dir)
    
    # CSV 파일 생성
    print("\nGenerating CSV file...")
    csv_path = os.path.join(output_dir, 'sbc_data.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['chrom', 'position', 'ref', 'alt', 'gt1', 'gt2', 'dp1', 'dp2', 'gq1', 'gq2', 'ad1', 'ad2', 'ab1', 'ab2', 'sample1_bases', 'sample2_bases', 'sample1_alt_ratio', 'sample2_alt_ratio', 'classification'])
        for variant in variants:
            writer.writerow([variant.chrom, variant.position, variant.ref, variant.alt, str(variant.sample1_gt), str(variant.sample2_gt), variant.sample1_dp, variant.sample2_dp, variant.sample1_gq, variant.sample2_gq, ':'.join(map(str, variant.sample1_ad)), ':'.join(map(str, variant.sample2_ad)), str(variant.sample1_ab), str(variant.sample2_ab), str(variant.sample1_base_counts), str(variant.sample2_base_counts), str(variant.sample1_alt_ratio), str(variant.sample2_alt_ratio), variant.quality_category])
    print(f"CSV file '{csv_path}' generated in output directory.")

    # SAM 검증 후 High-quality 스냅샷 생성
    for variant in variants:
        if variant.quality_category == "high_quality":
            bam1_file = pysam.AlignmentFile(bam1_path, "rb")
            bam2_file = pysam.AlignmentFile(bam2_path, "rb")
            reads1 = bam1_file.fetch(variant.chrom, max(1, variant.position - 21), variant.position + 20)
            reads2 = bam2_file.fetch(variant.chrom, max(1, variant.position - 21), variant.position + 20)
            variant.validate_with_sam(sample1_id, reads1, sample1_id, sample2_id, output_dir, reference_fasta)
            variant.validate_with_sam(sample2_id, reads2, sample1_id, sample2_id, output_dir, reference_fasta)
            bam1_file.close()
            bam2_file.close()

if __name__ == "__main__":
    main()