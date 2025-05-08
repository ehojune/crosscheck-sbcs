import configparser
from src.vcf_parser import parse_vcf
from src.sam_validator import validate_sam
from src.plotting import plot_chromosome_counts

# 설정 파일 읽기
config = configparser.ConfigParser()
config.read('config/config.ini')
vcf_path = config['paths']['vcf']
bam1_path = config['paths']['bam1']
bam2_path = config['paths']['bam2']
sample1_id = config['samples']['sample1']
sample2_id = config['samples']['sample2']
working_dir = config['paths']['working_dir']

# 출력 파일 경로
output_all = f"{working_dir}/sbc_all.txt"
output_high = f"{working_dir}/sbc_high_quality.txt"
output_low = f"{working_dir}/sbc_low_quality.txt"
plot_all = f"{working_dir}/sbc_all_chrom.png"
plot_high = f"{working_dir}/sbc_high_quality_chrom.png"

# VCF 파싱 및 SBC 추출
variants = parse_vcf(vcf_path, sample1_id, sample2_id)

# 검증 전 품질 분류
high_quality_pre = [v for v in variants if v.quality_category == "high_quality"]
low_quality_pre = [v for v in variants if v.quality_category == "low_quality"]
incorrect_pre = [v for v in variants if v.quality_category == "incorrect"]

print("=== Before SAM Validation ===")
print(f"Total SBCs: {len(variants)}")
print(f"High-quality SBCs: {len(high_quality_pre)}")
print(f"Low-quality SBCs: {len(low_quality_pre)}")
print(f"Incorrect SBCs: {len(incorrect_pre)}")

# SAM 파일로 검증
variants = validate_sam(variants, bam1_path, bam2_path, sample1_id, sample2_id)

# 검증 후 품질 재분류
for variant in variants:
    variant.quality_category = variant.classify_quality()

high_quality_post = [v for v in variants if v.quality_category == "high_quality"]
low_quality_post = [v for v in variants if v.quality_category == "low_quality"]
incorrect_post = [v for v in variants if v.quality_category == "incorrect"]

print("\n=== After SAM Validation ===")
print(f"Total SBCs: {len(variants)}")
print(f"High-quality SBCs: {len(high_quality_post)}")
print(f"Low-quality SBCs: {len(low_quality_post)}")
print(f"Incorrect SBCs: {len(incorrect_post)}")

# 파일로 저장
with open(output_all, 'w') as f:
    for v in variants:
        f.write(str(v) + '\n')
with open(output_high, 'w') as f:
    for v in high_quality_post:
        f.write(str(v) + '\n')
with open(output_low, 'w') as f:
    for v in low_quality_post:
        f.write(str(v) + '\n')

# 플롯 생성
plot_chromosome_counts(variants, plot_all, "All SBCs per Chromosome")
plot_chromosome_counts(high_quality_post, plot_high, "High-quality SBCs per Chromosome")