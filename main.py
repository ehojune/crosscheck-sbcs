import configparser
from src.vcf_parser import parse_vcf
from src.sam_validator import validate_sam
from src.variant import Variant






def main():
    # 설정 파일 읽기
    config = configparser.ConfigParser()
    config.read('config/config.ini')
    bam_path = config['paths']['bam']

    # VCF 파싱 및 Variant 객체 생성
    vcf_path = "path/to/vcf"
    variants = parse_vcf(vcf_path)  # Variant 객체 리스트 반환

    # BAM 파일로 검증
    validated_results = validate_sam(variants, bam_path)

    # 결과 출력
    for variant in validated_results:
        print(f"Variant {variant.position}: {'Valid' if variant.validate_variant() else 'Invalid'}")


if __name__ == '__main__':
    main()