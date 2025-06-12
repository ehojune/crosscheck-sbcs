from collections import Counter
import os
from igv_reports.report import create_report

class Variant:
    def __init__(self, chrom, position, ref, alt, sample1_gt, sample2_gt, 
                 sample1_dp, sample2_dp, sample1_gq, sample2_gq, sample1_ad, sample2_ad):
        self.chrom = chrom
        self.position = position
        self.ref = ref
        self.alt = alt if isinstance(alt, str) else alt[0]
        self.alt2 = alt[1] if isinstance(alt, tuple) and len(alt) > 1 else None
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
        self.sample1_base_counts = {}
        self.sample2_base_counts = {}
        self.sample1_alt_ratio = 0.0
        self.sample2_alt_ratio = 0.0
        self.sample1_major_count = 0
        self.sample1_minor_count = 0
        self.sample1_major_ratio = 0.0
        self.sample1_minor_ratio = 0.0
        self.sample2_major_count = 0
        self.sample2_minor_count = 0
        self.sample2_major_ratio = 0.0
        self.sample2_minor_ratio = 0.0
        self.initial_quality = self.classify_quality(use_sam=False)
        print(f"Initialized {self.chrom}:{self.position} - AB1={self.sample1_ab}, AB2={self.sample2_ab}, GT1={self.sample1_gt}, GT2={self.sample2_gt}")

    def calculate_ab(self, gt, ad, dp):
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
        # VCF 기반 조건
        ab_condition_high = (
            (self.sample1_gt != (0, 0) and self.sample1_ab is not None and self.sample1_ab >= 0.2) or
            (self.sample2_gt != (0, 0) and self.sample2_ab is not None and self.sample2_ab >= 0.2) or
            (self.sample1_gt == (0, 0) and self.sample1_ab is not None and self.sample1_ab in [0, 1]) or
            (self.sample2_gt == (0, 0) and self.sample2_ab is not None and self.sample2_ab in [0, 1])
        )
        ab_condition_moderate = (
            (self.sample1_gt != (0, 0) and self.sample1_ab is not None and self.sample1_ab >= 0.1) or
            (self.sample2_gt != (0, 0) and self.sample2_ab is not None and self.sample2_ab >= 0.1) or
            (self.sample1_gt == (0, 0) and self.sample1_ab is not None and self.sample1_ab in [0, 1]) or
            (self.sample2_gt == (0, 0) and self.sample2_ab is not None and self.sample2_ab in [0, 1])
        )

        if use_sam and (self.sample1_base_counts and self.sample2_base_counts):
            # SAM 기반 조건
            def is_valid_homozygous(sample_gt, base_counts, sample_id):
                if sample_gt == (0, 0):
                    return base_counts.get(self.ref, 0) == sum(base_counts.values()) and base_counts.get(self.alt, 0) == 0 and base_counts.get(self.alt2, 0) == 0
                elif sample_gt == (1, 1):
                    return base_counts.get(self.alt, 0) == sum(base_counts.values()) and base_counts.get(self.ref, 0) == 0
                elif sample_gt == (2, 2) and self.alt2:
                    return base_counts.get(self.alt2, 0) == sum(base_counts.values()) and base_counts.get(self.ref, 0) == 0 and base_counts.get(self.alt, 0) == 0
                return False

            def is_valid_heterozygous(sample_gt, base_counts, sample_id):
                if sample_gt in [(0, 1), (1, 0)]:
                    return (base_counts.get(self.ref, 0) > 0 and base_counts.get(self.alt, 0) > 0 and
                            sum(base_counts.values()) == base_counts.get(self.ref, 0) + base_counts.get(self.alt, 0))
                return False

            sam_condition_high = (
                (is_valid_homozygous(self.sample1_gt, self.sample1_base_counts, "sample1") and
                 is_valid_heterozygous(self.sample2_gt, self.sample2_base_counts, "sample2") and
                 self.sample1_major_count >= 10 and self.sample2_major_count >= 10) or
                (is_valid_heterozygous(self.sample1_gt, self.sample1_base_counts, "sample1") and
                 is_valid_homozygous(self.sample2_gt, self.sample2_base_counts, "sample2") and
                 self.sample1_major_count >= 10 and self.sample2_major_count >= 10) or
                (is_valid_homozygous(self.sample1_gt, self.sample1_base_counts, "sample1") and
                 is_valid_homozygous(self.sample2_gt, self.sample2_base_counts, "sample2") and
                 self.sample1_major_count >= 10 and self.sample2_major_count >= 10)
            )

            if (self.sample1_dp and self.sample1_dp >= 50 and self.sample2_dp and self.sample2_dp >= 50 and
                self.sample1_gq and self.sample1_gq >= 50 and self.sample2_gq and self.sample2_gq >= 50 and
                ab_condition_high and sam_condition_high):
                print(f"{self.chrom}:{self.position} classified as high_quality")
                return "high_quality"

            sam_condition_moderate = (
                (is_valid_homozygous(self.sample1_gt, self.sample1_base_counts, "sample1") or
                 is_valid_heterozygous(self.sample1_gt, self.sample1_base_counts, "sample1")) and
                (is_valid_homozygous(self.sample2_gt, self.sample2_base_counts, "sample2") or
                 is_valid_heterozygous(self.sample2_gt, self.sample2_base_counts, "sample2")) and
                self.sample1_major_count >= 5 and self.sample2_major_count >= 5
            )

            if (self.sample1_dp and self.sample1_dp >= 20 and self.sample2_dp and self.sample2_dp >= 20 and
                self.sample1_gq and self.sample1_gq >= 30 and self.sample2_gq and self.sample2_gq >= 30 and
                ab_condition_moderate and sam_condition_moderate):
                print(f"{self.chrom}:{self.position} classified as moderate_quality")
                return "moderate_quality"
        else:
            if (self.sample1_dp and self.sample1_dp >= 20 and self.sample2_dp and self.sample2_dp >= 20 and
                self.sample1_gq and self.sample1_gq >= 30 and self.sample2_gq and self.sample2_gq >= 30 and
                ab_condition_high):
                print(f"{self.chrom}:{self.position} classified as high_quality (pre-SAM)")
                return "high_quality"

            if (self.sample1_dp and self.sample1_dp >= 10 and self.sample2_dp and self.sample2_dp >= 10 and
                self.sample1_gq and self.sample1_gq >= 20 and self.sample2_gq and self.sample2_gq >= 20 and
                ab_condition_moderate):
                print(f"{self.chrom}:{self.position} classified as moderate_quality (pre-SAM)")
                return "moderate_quality"

        print(f"{self.chrom}:{self.position} classified as unclassified")
        return "unclassified"

    def validate_with_sam(self, sample_id, reads, sample1_id, sample2_id, output_dir=None, bam_path=None, reference_fasta=None):
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
            read_count += 1

        sorted_bases = base_counts.most_common()
        total_bases = sum(base_counts.values())
        major_count = sorted_bases[0][1] if sorted_bases else 0
        minor_count = sorted_bases[1][1] if len(sorted_bases) > 1 else 0
        major_ratio = major_count / total_bases if total_bases > 0 else 0.0
        minor_ratio = minor_count / total_bases if total_bases > 0 else 0.0

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

        if output_dir and self.quality_category == "high_quality" and reference_fasta:
            self.generate_igv_snapshot(sample_id, bam_path, output_dir, sample1_id, sample2_id, reference_fasta)

        print(f"Validated {self.chrom}:{self.position} for {sample_id} - Bases: {base_counts}, "
              f"Major: {major_count} ({major_ratio:.2f}), Minor: {minor_count} ({minor_ratio:.2f}), "
              f"ALT_ratio: {alt_ratio}, Total_reads: {read_count}")

    def generate_igv_snapshot(self, sample_id, bam_path, output_dir, sample1_id, sample2_id, reference_fasta):
        region = f"{self.chrom}:{max(1, self.position - 20)}-{self.position + 20}"
        output_html = os.path.join(output_dir, f"{self.chrom}_{self.position}_{self.ref}>{self.alt}.html")
        tracks = [bam_path] if sample_id == sample1_id else [bam_path.replace(sample1_id, sample2_id)]
        
        create_report(
            sites=[{"chrom": self.chrom, "start": max(1, self.position - 20), "end": self.position + 20}],
            fasta=reference_fasta,
            tracks=tracks,
            flanking=20,
            output=output_html
        )
        print(f"IGV snapshot generated for {sample_id} at {output_html}")

    def __str__(self):
        return (f"{self.chrom}\t{self.position}\t{self.ref}>{self.alt}\t"
                f"GT1={self.sample1_gt}\tGT2={self.sample2_gt}\t"
                f"Initial_Quality={self.initial_quality}\tQuality={self.quality_category}\t"
                f"Sample1_Bases={self.sample1_base_counts}\tSample2_Bases={self.sample2_base_counts}")