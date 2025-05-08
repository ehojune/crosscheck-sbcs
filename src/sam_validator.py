import pysam
from .variant import Variant

def validate_sam(variants, bam1_path, bam2_path, sample1_id, sample2_id):
    """샘플별 BAM 파일에서 변이 검증."""
    try:
        bam1_file = pysam.AlignmentFile(bam1_path, "rb")
        bam2_file = pysam.AlignmentFile(bam2_path, "rb")
    except Exception as e:
        print(f"Error opening BAM files: {e}")
        return []

    for variant in variants:
        # 샘플1 검증
        expected_base1 = variant.ref if variant.sample1_gt == (0, 0) else variant.alt
        reads1 = bam1_file.fetch(variant.chrom, variant.position - 1, variant.position)
        variant.validate_with_sam(sample1_id, reads1, expected_base1)

        # 샘플2 검증
        expected_base2 = variant.ref if variant.sample2_gt == (0, 0) else variant.alt
        reads2 = bam2_file.fetch(variant.chrom, variant.position - 1, variant.position)
        variant.validate_with_sam(sample2_id, reads2, expected_base2)

    bam1_file.close()
    bam2_file.close()
    return variants