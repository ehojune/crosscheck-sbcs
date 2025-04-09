import configparser
from src.vcf_parser import parse_vcf
from src.plotting import plot_chromosome_counts

# 설정 파일 읽기
config = configparser.ConfigParser()
config.read('config/config.ini')
vcf_path = config['paths']['vcf']
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

# 품질별 분류
high_quality = [v for v in variants if v.quality_category == "high_quality"]
low_quality = [v for v in variants if v.quality_category == "low_quality"]
incorrect = [v for v in variants if v.quality_category == "incorrect"]

# 결과 출력
print(f"Total SBCs: {len(variants)}")
print(f"High-quality SBCs: {len(high_quality)}")
print(f"Low-quality SBCs: {len(low_quality)}")
print(f"Incorrect SBCs: {len(incorrect)}")

# 파일로 저장
with open(output_all, 'w') as f:
    for v in variants:
        f.write(str(v) + '\n')
with open(output_high, 'w') as f:
    for v in high_quality:
        f.write(str(v) + '\n')
with open(output_low, 'w') as f:
    for v in low_quality:
        f.write(str(v) + '\n')

# 플롯 생성
plot_chromosome_counts(variants, plot_all, "All SBCs per Chromosome")
plot_chromosome_counts(high_quality, plot_high, "High-quality SBCs per Chromosome")