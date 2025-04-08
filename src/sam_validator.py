import pysam
from .variant import Variant

def validate_sam(variants, bam_path):
    """
    BAM 파일에서 Variant 위치에 해당하는 alignment를 확인해 검증.
    variants: Variant 객체 리스트
    bam_path: BAM 파일 경로 (인덱스된 .bai 파일이 있어야 함)
    """
    # BAM 파일 열기 (인덱스가 필요하므로 .bai 파일이 있어야 함)
    bam_file = pysam.AlignmentFile(bam_path, "rb")

    validated_variants = []
    for variant in variants:
        # 특정 위치의 읽기 가져오기
        reads = bam_file.fetch(variant.chrom, variant.position - 1, variant.position)
        # 참고: fetch는 0-based 좌표를 사용하므로 position - 1

        # 읽기에서 변이 확인 (간단한 예시로 ref/alt 비교)
        is_valid = False
        for read in reads:
            if read.reference_start <= variant.position - 1 <= read.reference_end:
                # 해당 위치의 base 확인 (read.query_sequence 사용)
                pos_in_read = variant.position - 1 - read.reference_start
                if pos_in_read < len(read.query_sequence):
                    base = read.query_sequence[pos_in_read]
                    if base == variant.alt:
                        is_valid = True
                        break  # 일치하면 더 볼 필요 없음

        variant.is_valid = is_valid
        validated_variants.append(variant)

    bam_file.close()
    return validated_variants