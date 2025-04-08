from src.variant import Variant

def test_variant_initialization():
    v = Variant(chrom="chr1", position=100, ref="A", alt="T")
    assert v.chrom == "chr1"
    assert v.position == 100
    assert v.ref == "A"
    assert v.alt == "T"
    assert v.is_valid is None

def test_validate_variant():
    v = Variant(chrom="chr1", position=100, ref="A", alt="T")
    result = v.validate_variant()
    assert result is True  # 더미 로직이 True 반환
    assert v.is_valid is True