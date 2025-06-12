import matplotlib.pyplot as plt
import numpy as np
import os

def plot_pre_sam_histogram(variants, output_dir):
    qualities = [v.initial_quality for v in variants]
    unique_quals, counts = np.unique(qualities, return_counts=True)
    plt.bar(unique_quals, counts)
    plt.title('Pre-SAM Classification Histogram')
    plt.xlabel('Quality')
    plt.ylabel('Count')
    plt.savefig(os.path.join(output_dir, 'pre_sam_classification_histogram.png'))
    plt.close()

def plot_post_sam_histogram(variants, output_dir):
    qualities = [v.quality_category for v in variants]
    unique_quals, counts = np.unique(qualities, return_counts=True)
    plt.bar(unique_quals, counts)
    plt.title('Post-SAM Classification Histogram')
    plt.xlabel('Quality')
    plt.ylabel('Count')
    plt.savefig(os.path.join(output_dir, 'post_sam_classification_histogram.png'))
    plt.close()

def plot_hq_histogram(variants, output_dir):
    hq_variants = [v for v in variants if v.quality_category == 'high_quality']
    if hq_variants:
        qualities = ['High'] * len(hq_variants)  # 단순히 High-quality 수 표시
        counts = [len(hq_variants)]
        plt.bar(qualities, counts)
        plt.title('High-Quality Histogram')
        plt.xlabel('Quality')
        plt.ylabel('Count')
        plt.savefig(os.path.join(output_dir, 'hq_classification_histogram.png'))
        plt.close()
    else:
        print("No High-quality variants to plot.")