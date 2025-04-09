import matplotlib.pyplot as plt
from collections import defaultdict

def plot_chromosome_counts(variants, output_image, title):
    """Chromosome별 변이 수를 플롯."""
    chrom_counts = defaultdict(int)
    for variant in variants:
        chrom_counts[variant.chrom] += 1

    chroms = list(chrom_counts.keys())
    counts = list(chrom_counts.values())

    plt.figure(figsize=(10, 6))
    plt.bar(chroms, counts, color='blue')
    plt.xlabel('Chromosome')
    plt.ylabel('SBC Count')
    plt.title(title)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_image)
    plt.close()