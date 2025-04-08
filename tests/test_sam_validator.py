import pysam
from src.sam_validator import validate_sam
from src.variant import Variant

def test_validate_sam():
    # 테스트용 BAM 파일이 필요하므로 더미 데이터로 간단히
    variants = [Variant(chrom="chr1", position=100, ref="A", alt="T")]
    # 실제 BAM 파일 경로 대신 가짜 경로로 테스트 (실제로는 파일 필요)
    validated = validate_sam(variants, "/data2/Twin/scripts/sample.temp.bam")
    assert len(validated) == 1
    # 더미 로직에서는 검증 로직이 없으므로 추가 테스트 필요