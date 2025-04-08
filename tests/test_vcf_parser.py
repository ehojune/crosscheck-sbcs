from src.vcf_parser import parse_vcf

def test_parse_vcf():
    # 실제 VCF 파일 대신 더미 데이터로 테스트
    variants = parse_vcf("dummy_path")
    assert len(variants) == 2  # 예시에서 2개 반환
    assert variants[0].position == 100
    assert variants[1].alt == "C"