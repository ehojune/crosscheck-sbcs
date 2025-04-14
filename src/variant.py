class Variant:
    def __init__(self, chrom, position, ref, alt, sample1_gt, sample2_gt, 
                 sample1_dp, sample2_dp, sample1_gq, sample2_gq, sample1_ad, sample2_ad):
        self.chrom = chrom
        self.position = position
        self.ref = ref
        self.alt = alt
        self.sample1_gt = sample1_gt
        self.sample2_gt = sample2_gt
        self.sample1_dp = sample1_dp
        self.sample2_dp = sample2_dp
        self.sample1_gq = sample1_gq
        self.sample2_gq = sample2_gq
        self.sample1_ad = sample1_ad  # AD 필드 추가
        self.sample2_ad = sample2_ad
        self.sample1_ab = self.calculate_ab(sample1_gt, sample1_ad, sample1_dp)
        self.sample2_ab = self.calculate_ab(sample2_gt, sample2_ad, sample2_dp)
        self.quality_category = self.classify_quality()

    def calculate_ab(self, gt, ad, dp):
        """AD와 DP를 기반으로 AB 계산."""
        if dp == 0 or gt == (0, 0):  # hom-ref는 AB 계산 불가
            return None
        if gt in [(0, 1), (1, 0)] and len(ad) >= 2:  # het
            return ad[1] / dp
        if gt == (1, 1) and len(ad) >= 2:  # hom-alt
            return ad[1] / dp
        return None

    def classify_quality(self):
        """논문 기준으로 SBC를 High-quality, Low-quality, Incorrect로 분류."""
        # High-quality: DP ≥ 20, GQ ≥ 30, AB ≥ 0.2
        if (self.sample1_dp >= 20 and self.sample2_dp >= 20 and
            self.sample1_gq >= 30 and self.sample2_gq >= 30 and
            (self.sample1_ab is not None and self.sample1_ab >= 0.2 or
             self.sample2_ab is not None and self.sample2_ab >= 0.2)):
            return "high_quality"
        
        # Low-quality: DP ≥ 10, GQ ≥ 20, AB ≥ 0.1
        elif (self.sample1_dp >= 10 and self.sample2_dp >= 10 and
              self.sample1_gq >= 20 and self.sample2_gq >= 20 and
              (self.sample1_ab is not None and self.sample1_ab >= 0.1 or
               self.sample2_ab is not None and self.sample2_ab >= 0.1)):
            return "low_quality"
        
        # Clear incorrect
        return "incorrect"

    def __str__(self):
        return (f"{self.chrom}:{self.position} {self.ref}>{self.alt} "
                f"GT1={self.sample1_gt} GT2={self.sample2_gt} "
                f"Quality={self.quality_category}")