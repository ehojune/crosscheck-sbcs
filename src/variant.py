class Variant:
    def __init__(self, chrom: str, pos: int, ref: str, alt: str, twin1_GT: list, twin2_GT: list, 
                    twin1_AD: list, twin2_AD: list, twin1_MAF: float, twin2_MAF: float ):
        """Variant 객체 초기화."""
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.is_valid = None  # 검증 결과 저장용
        self.twin1_GT = twin1_GT
        self.twin2_GT = twin2_GT
        self.twin1_AD = twin1_AD
        self.twin2_AD = twin2_AD
        self.twin1_MAF = twin1_MAF
        self.twin2_MAF = twin2_MAF

    def validate_variant(self):
        """이 변이가 SAM 파일과 일치하는지 검증하는 메서드."""
        # 실제로는 SAM 데이터와 비교 로직이 들어감
        # 지금은 더미로 True 반환
        self.is_valid = True
        return self.is_valid

    @classmethod
    def from_dict(cls, variant_dict):
        """딕셔너리에서 Variant 객체를 생성하는 클래스 메서드 예시."""
        return cls(
            chrom=variant_dict["chrom"],
            position=variant_dict["pos"],
            ref=variant_dict["ref"],
            alt=variant_dict["alt"]
        )

    def __str__(self):
        return f"{self.chrom}:{self.position} {self.ref}>{self.alt}"