from .variant import Variant  # 같은 패키지 내 variant 모듈에서 가져옴

def check_SBC(line):
	if line != variant.twin2_GT:
		return True
	return False







def parse_vcf(vcf_path):
    """
    VCF 파일을 파싱해서 Variant 객체 리스트를 반환.
    """

    vcf = open(vcf_path, 'r')


	for line in vcf_path











    # 예시로 더미 데이터 반환
    variants = [
        Variant(chrom="chr1", position=100, ref="A", alt="T"),
        Variant(chrom="chr1", position=200, ref="G", alt="C")
    ]
    return variants






