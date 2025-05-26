from collections import Counter

class Variant:
    def __init__(self, chrom, position, ref, alt, sample1_gt, sample2_gt, 
                 sample1_dp, sample2_dp, sample1_gq, sample2_gq, sample1_ad, sample2_ad):
        self.chrom = chrom
        self.position = position
        self.ref = ref
        self.alt = alt if isinstance(alt, str) else alt[0]  # Multi-allelic ALT의 첫 번째만 사용
        self.alt2 = alt[1] if isinstance(alt, tuple) and len(alt) > 1 else None  # 두 번째 ALT
        self.sample1_gt = sample1_gt
        self.sample2_gt = sample2_gt
        self.sample1_dp = sample1_dp
        self.sample2_dp = sample2_dp
        self.sample1_gq = sample1_gq
        self.sample2_gq = sample2_gq
        self.sample1_ad = sample1_ad
        self.sample2_ad = sample2_ad
        self.sample1_ab = self.calculate_ab(sample1_gt, sample1_ad, sample1_dp)
        self.sample2_ab = self.calculate_ab(sample2_gt, sample2_ad, sample2_dp)
        self.sample1_validated = False
        self.sample2_validated = False
        self.sample1_read_count = 0
        self.sample2_read_count = 0
        self.sample1_base_counts = {}  # SAM base 분포 (예: {'G': 20, 'A': 1})
        self.sample2_base_counts = {}
        self.sample1_alt_ratio = 0.0
        self.sample2_alt_ratio = 0.0
        self.quality_category = self.classify_quality(use_sam=False)
        print(f"Initialized {self.chrom}:{self.position} - AB1={self.sample1_ab}, AB2={self.sample2_ab}, GT1={self.sample1_gt}, GT2={self.sample2_gt}")

    def calculate_ab(self, gt, ad, dp):
        """AD와 DP를 기반으로 AB 계산."""
        if dp == 0 or gt == (0, 0) or not ad or len(ad) < 2:
            return None
        try:
            if gt in [(0, 1), (1, 0), (1, 1)]:
                return ad[1] / dp
            if gt in [(0, 2), (1, 2)] and len(ad) > 2:
                return ad[2] / dp
        except (TypeError, IndexError):
            return None
        return None

    def classify_quality(self, use_sam=True):
        """SBC를 High-quality, Moderate-quality, Low-quality로 분류."""
        # AB 조건
        ab_condition_high = (
            (self.sample1_gt != (0, 0) and self.sample1_ab is not None and self.sample1_ab >= 0.2) or  # het
            (self.sample2_gt != (0, 0) and self.sample2_ab is not None and self.sample2_ab >= 0.2) or  # het
            (self.sample1_gt == (0, 0) and self.sample1_ab is not None and self.sample1_ab in [0, 1]) or  # hom
            (self.sample2_gt == (0, 0) and self.sample2_ab is not None and self.sample2_ab in [0, 1])  # hom
        )
        ab_condition_moderate = (
            (self.sample1_gt != (0, 0) and self.sample1_ab is not None and self.sample1_ab >= 0.1) or  # het
            (self.sample2_gt != (0, 0) and self.sample2_ab is not None and self.sample2_ab >= 0.1) or  # het
            (self.sample1_gt == (0, 0) and self.sample1_ab is not None and self.sample1_ab in [0, 1]) or  # hom
            (self.sample2_gt == (0, 0) and self.sample2_ab is not None and self.sample2_ab in [0, 1])  # hom
        )

        # SAM 기반 조건
        sam_condition_high = (
            (self.sample1_gt != (0, 0) and self.sample1_alt_ratio >= 0.2 and self.sample1_major_count >= 5 and
            self.sample2_alt_ratio >= 0.2 and self.sample2_major_count >= 5 and self.sample2_minor_count == 0) or  # het
            (self.sample1_gt == (0, 0) and self.sample1_alt_ratio in [0, 1] and self.sample1_major_count >= 5 and
            self.sample2_alt_ratio not in [0, 1] and self.sample2_major_count >= 5 and self.sample2_minor_count >= 1) or  # hom vs het
            (self.sample1_gt == (0, 0) and self.sample1_alt_ratio in [0, 1] and self.sample1_major_count >= 5 and
            self.sample2_alt_ratio in [0, 1] and self.sample2_major_count >= 5 and self.sample2_minor_count == 0)  # hom vs hom
        )
        sam_condition_moderate = (
            (self.sample1_gt != (0, 0) and self.sample1_alt_ratio >= 0.1 and self.sample1_major_count >= 1 and
            self.sample2_alt_ratio >= 0.1 and self.sample2_major_count >= 1) or  # het
            (self.sample1_gt == (0, 0) and self.sample1_alt_ratio in [0, 1] and self.sample1_major_count >= 1 and
            self.sample2_alt_ratio not in [0, 1] and self.sample2_major_count >= 1) or  # hom vs het
            (self.sample1_gt == (0, 0) and self.sample1_alt_ratio in [0, 1] and self.sample1_major_count >= 1 and
            self.sample2_alt_ratio in [0, 1] and self.sample2_major_count >= 1)  # hom vs hom
        )

        # High-quality
        if (self.sample1_dp and self.sample1_dp >= 20 and self.sample2_dp and self.sample2_dp >= 20 and
            self.sample1_gq and self.sample1_gq >= 30 and self.sample2_gq and self.sample2_gq >= 30 and
            ab_condition_high):
            if not use_sam or (self.sample1_validated and self.sample2_validated and sam_condition_high):
                print(f"{self.chrom}:{self.position} classified as high_quality")
                return "high_quality"

        # Moderate-quality
        if (self.sample1_dp and self.sample1_dp >= 10 and self.sample2_dp and self.sample2_dp >= 10 and
            self.sample1_gq and self.sample1_gq >= 20 and self.sample2_gq and self.sample2_gq >= 20 and
            ab_condition_moderate):
            if not use_sam or (self.sample1_validated and self.sample2_validated and sam_condition_moderate):
                print(f"{self.chrom}:{self.position} classified as moderate_quality")
                return "moderate_quality"

        print(f"{self.chrom}:{self.position} classified as low_quality")
        return "low_quality"

    def validate_with_sam(self, sample_id, reads, sample1_id, sample2_id):
        """샘플별 SAM 읽기로 변이 검증 (base 분포 분석)."""
        base_counts = Counter()
        read_count = 0
        for read in reads:
            if read.mapping_quality >= 20:
                if read.reference_start <= self.position - 1 <= read.reference_end:
                    pos_in_read = self.position - 1 - read.reference_start
                    if 0 <= pos_in_read < len(read.query_sequence):
                        base_qual = read.query_qualities[pos_in_read] if read.query_qualities else 0
                        if base_qual >= 20:
                            base = read.query_sequence[pos_in_read]
                            base_counts[base] += 1
                        else:
                            print(f"Skipped low quality base at {self.chrom}:{self.position} for {sample_id}: qual={base_qual}")
                else:
                    print(f"Read out of range at {self.chrom}:{self.position} for {sample_id}: {read.reference_start}-{read.reference_end}")
            else:
                print(f"Skipped low MAPQ read at {self.chrom}:{self.position} for {sample_id}: MAPQ={read.mapping_quality}")
            read_count += 1  # 총 reads 수

        # Major/minor read 계산
        sorted_bases = base_counts.most_common()
        total_bases = sum(base_counts.values())
        major_count = sorted_bases[0][1] if sorted_bases else 0
        minor_count = sorted_bases[1][1] if len(sorted_bases) > 1 else 0
        major_ratio = major_count / total_bases if total_bases > 0 else 0.0
        minor_ratio = minor_count / total_bases if total_bases > 0 else 0.0

        # ALT 비율 계산
        alt_count = base_counts.get(self.alt, 0)
        alt2_count = base_counts.get(self.alt2, 0) if self.alt2 else 0
        alt_ratio = (alt_count + alt2_count) / total_bases if total_bases > 0 else 0.0

        if sample_id == sample1_id:
            self.sample1_read_count = read_count
            self.sample1_major_count = major_count
            self.sample1_minor_count = minor_count
            self.sample1_major_ratio = major_ratio
            self.sample1_minor_ratio = minor_ratio
            self.sample1_validated = alt_count > 0 or alt2_count > 0 or base_counts.get(self.ref, 0) > 0
            self.sample1_base_counts = dict(base_counts)
            self.sample1_alt_ratio = alt_ratio
        elif sample_id == sample2_id:
            self.sample2_read_count = read_count
            self.sample2_major_count = major_count
            self.sample2_minor_count = minor_count
            self.sample2_major_ratio = major_ratio
            self.sample2_minor_ratio = minor_ratio
            self.sample2_validated = alt_count > 0 or alt2_count > 0 or base_counts.get(self.ref, 0) > 0
            self.sample2_base_counts = dict(base_counts)
            self.sample2_alt_ratio = alt_ratio

        print(f"Validated {self.chrom}:{self.position} for {sample_id} - Bases: {base_counts}, "
            f"Major: {major_count} ({major_ratio:.2f}), Minor: {minor_count} ({minor_ratio:.2f}), "
            f"ALT_ratio: {alt_ratio}, Total_reads: {read_count}")

    def __str__(self):
        return (f"{self.chrom}:{self.position} {self.ref}>{self.alt} "
                f"GT1={self.sample1_gt} GT2={self.sample2_gt} "
                f"Quality={self.quality_category} "
                f"Sample1_Bases={self.sample1_base_counts} Sample2_Bases={self.sample2_base_counts}")