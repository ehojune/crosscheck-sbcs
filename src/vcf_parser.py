import pysam
from .variant import Variant

def parse_vcf(vcf_file, sample1_id, sample2_id):
	"""VCF 파일에서 쌍둥이 간 SBC를 추출."""
	vcf_in = pysam.VariantFile(vcf_file)
	variants = []

	for record in vcf_in:
		sample1 = record.samples[sample1_id]
		sample2 = record.samples[sample2_id]
		gt1, gt2 = sample1['GT'], sample2['GT']

		# 쌍둥이 간 GT가 다르면 SBC로 간주
		if gt1 != gt2 and '.' not in gt1 and '.' not in gt2:
			variant = Variant(
				chrom=record.chrom,
				position=record.pos,
				ref=record.ref,
				alt=record.alts[0],  # 단일 alt 가정
				sample1_gt=gt1,
				sample2_gt=gt2,
				sample1_dp=sample1.get('DP', 0),
				sample2_dp=sample2.get('DP', 0),
				sample1_dp=sample1.get('AD', 0),
				sample2_dp=sample2.get('AD', 0),
				sample1_gq=sample1.get('GQ', 0),
				sample2_gq=sample2.get('GQ', 0),
				sample1_ab=sample1.get('AB', [0])[0] if 'AB' in sample1 else None,
				sample2_ab=sample2.get('AB', [0])[0] if 'AB' in sample2 else None
			)
			variants.append(variant)

	return variants