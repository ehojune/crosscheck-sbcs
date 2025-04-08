def parse_vcf(vcf_path):
    """
    VCF 파일을 파싱해서 Variant 객체 리스트를 반환.
    실제로는 pysam이나 vcfpy 같은 라이브러리를 사용할 가능성이 높음.
    """
    from .variant import Variant  # 같은 패키지 내 variant 모듈에서 가져옴
    
    # 예시로 더미 데이터 반환
    variants = [
        Variant(chrom="chr1", position=100, ref="A", alt="T"),
        Variant(chrom="chr1", position=200, ref="G", alt="C")
    ]
    return variants