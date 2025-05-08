import matplotlib.pyplot as plt
from collections import defaultdict
import os

def plot_pre_sam_histogram(variants, output_dir):
    """SAM 검증 전 염색체별 SBC 분류 개수를 히스토그램으로 시각화."""
    chr_counts = defaultdict(lambda: {'high_quality': 0, 'low_quality': 0, 'incorrect': 0, 'total': 0})
    for variant in variants:
        chr = variant.chrom
        quality = variant.quality_category
        chr_counts[chr][quality] += 1
        chr_counts[chr]['total'] += 1
    
    chromosomes = sorted(chr_counts.keys())
    hq_counts = [chr_counts[chr]['high_quality'] for chr in chromosomes]
    lq_counts = [chr_counts[chr]['low_quality'] for chr in chromosomes]
    incorrect_counts = [chr_counts[chr]['incorrect'] for chr in chromosomes]
    
    fig, ax = plt.subplots(figsize=(8, 6))
    bar_width = 0.35
    index = range(len(chromosomes))
    
    p1 = ax.bar(index, hq_counts, bar_width, label='High-quality', color='#1f77b4')
    p2 = ax.bar(index, lq_counts, bar_width, bottom=hq_counts, label='Low-quality', color='#ff7f0e')
    p3 = ax.bar(index, incorrect_counts, bar_width, bottom=[x + y for x, y in zip(hq_counts, lq_counts)], label='Incorrect', color='#2ca02c')
    
    ax.set_xlabel('Chromosomes')
    ax.set_ylabel('Number of SBCs')
    ax.set_title('SBC Classification by Chromosome (Before SAM Validation)')
    ax.set_xticks(index)
    ax.set_xticklabels(chromosomes)
    ax.legend()
    
    output_path = os.path.join(output_dir, 'pre_sam_classification_histogram.png')
    plt.savefig(output_path)
    plt.close()

def plot_post_sam_histogram(variants, output_dir):
    """SAM 검증 후 염색체별 SBC 분류 개수를 히스토그램으로 시각화."""
    chr_counts = defaultdict(lambda: {'high_quality': 0, 'low_quality': 0, 'incorrect': 0, 'total': 0})
    for variant in variants:
        chr = variant.chrom
        quality = variant.quality_category
        chr_counts[chr][quality] += 1
        chr_counts[chr]['total'] += 1
    
    chromosomes = sorted(chr_counts.keys())
    hq_counts = [chr_counts[chr]['high_quality'] for chr in chromosomes]
    lq_counts = [chr_counts[chr]['low_quality'] for chr in chromosomes]
    incorrect_counts = [chr_counts[chr]['incorrect'] for chr in chromosomes]
    
    fig, ax = plt.subplots(figsize=(8, 6))
    bar_width = 0.35
    index = range(len(chromosomes))
    
    p1 = ax.bar(index, hq_counts, bar_width, label='High-quality', color='#1f77b4')
    p2 = ax.bar(index, lq_counts, bar_width, bottom=hq_counts, label='Low-quality', color='#ff7f0e')
    p3 = ax.bar(index, incorrect_counts, bar_width, bottom=[x + y for x, y in zip(hq_counts, lq_counts)], label='Incorrect', color='#2ca02c')
    
    ax.set_xlabel('Chromosomes')
    ax.set_ylabel('Number of SBCs')
    ax.set_title('SBC Classification by Chromosome (After SAM Validation)')
    ax.set_xticks(index)
    ax.set_xticklabels(chromosomes)
    ax.legend()
    
    output_path = os.path.join(output_dir, 'post_sam_classification_histogram.png')
    plt.savefig(output_path)
    plt.close()