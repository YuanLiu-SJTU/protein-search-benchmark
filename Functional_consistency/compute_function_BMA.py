#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import random
import argparse
from collections import defaultdict
from tqdm import tqdm
import subprocess
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm

from goatools.base import get_godag
from goatools.semsim.termwise.wang import SsWang

random.seed(42)
topk = 10
COLORS = cm.Oranges(0.5)

# -------------------------------
# Load GO DAG and annotations
# -------------------------------
with open('./GO.json', 'r') as f:
    go_anno = json.load(f)

godag = get_godag("./go-basic.obo", optional_attrs={'relationship'})


def simscore(goids, go_a, go_b, godag):
    """Compute Wang semantic similarity between two GO terms"""
    relationships = {'part_of'}
    wang_r1 = SsWang(goids, godag, relationships)
    return wang_r1.get_sim(go_a, go_b)


def read_search_result_tsv(file_path):
    """Parse search result TSV into dict: query -> list of targets"""
    results = defaultdict(list)
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            q, t = parts[:2]
            q = q[:-4] if q.endswith('.pdb') else q
            t = t[:-4] if t.endswith('.pdb') else t
            results[q].append(t)
    return results


def random_padding(result_dict, query_list, target_list):
    """Ensure each query has topk targets by adding random targets if necessary"""
    padding_dict = {}
    for q in query_list:
        if q not in result_dict:
            padding_dict[q] = random.sample(target_list, k=topk)
        elif len(result_dict[q]) < topk:
            available = list(set(target_list) - set(result_dict[q]))
            padding_dict[q] = result_dict[q] + random.sample(available, k=topk - len(result_dict[q]))
        else:
            padding_dict[q] = result_dict[q][:topk]
    return padding_dict


def get_BMA_result(search_result, query_list, target_list):
    """Compute GO semantic similarity using Best Match Average (BMA) with topk hits"""
    padded_results = random_padding(search_result, query_list, target_list)
    scores = []

    for q, t_list in tqdm(padded_results.items(), desc='Compute BMA'):
        if q not in go_anno:
            continue
        go_q = go_anno[q]
        for t in t_list:
            if t not in go_anno:
                continue
            go_t = go_anno[t]
            # BMA: query -> target
            bma_qt = np.mean([
                max([simscore([go1, go2], go1, go2, godag) for go2 in go_t])
                for go1 in go_q
            ])
            # BMA: target -> query
            bma_tq = np.mean([
                max([simscore([go1, go2], go1, go2, godag) for go1 in go_q])
                for go2 in go_t
            ])
            scores.append((q, t, (bma_qt + bma_tq) / 2))
    return scores


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--search_result', required=True, help='Search result TSV file')
    parser.add_argument('--output_dir', default='figures', help='Directory to save boxplot')
    args = parser.parse_args()

    search_result_file = args.search_result
    search_result = read_search_result_tsv(search_result_file)

    # Query list and target list
    query_list = list(search_result.keys())
    # Target list is all keys in GO annotation
    target_list = list(go_anno.keys())

    print(f"Number of queries: {len(query_list)}")
    print(f"Number of targets: {len(target_list)}")
    print(f"Processing top-{topk} hits with random padding...")

    # Compute BMA scores
    scores_data = get_BMA_result(search_result, query_list, target_list)
    df = pd.DataFrame(scores_data, columns=['Query', 'Target', 'GO Similarity'])

    # Boxplot
    os.makedirs(args.output_dir, exist_ok=True)
    plt.figure(figsize=(8, 5), dpi=300)
    sns.boxplot(y='GO Similarity', data=df, color=COLORS)
    plt.title(f'GO Semantic Similarity (top-{topk} hits)')
    plt.ylabel('BMA Similarity Score')
    plt.tight_layout()
    save_path = os.path.join(args.output_dir, f"{os.path.basename(search_result_file).replace('.tsv','')}.png")
    plt.savefig(save_path)
    plt.close()

    print(f"Boxplot saved to {save_path}")


if __name__ == '__main__':
    main()
