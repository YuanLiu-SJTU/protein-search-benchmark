#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
from collections import defaultdict
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import pickle

LEVELS = ['Class', 'Architecture', 'Topology', 'Homologous superfamily']

def read_search_result_tsv(file_path):
    """Read search results from a TSV file"""
    search_result = defaultdict(list)
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            query, target = parts[:2]
            search_result[query].append(target)
    return search_result

def read_cath_labels(cath_file):
    """Read CATH labels"""
    label_dict = {}
    with open(cath_file,'r') as f:
        for line in f:
            parts = line.strip().split()
            entry = parts[0]
            label_dict[entry] = parts[1:5]  # Class, Architecture, Topology, Homologous superfamily
    return label_dict

def compute_tp_count(label_dict, save_file=None):
    """Compute total TP for each query at each level (excluding self)"""
    if save_file and os.path.exists(save_file):
        with open(save_file,'rb') as p:
            return pickle.load(p)

    tp_dict = {}
    print("Computing total TP for each query...")
    for pdb in tqdm(label_dict.keys(), desc='TP count'):
        tp_list = [0]*4
        for pdb2 in label_dict.keys():
            if pdb2 == pdb:
                continue
            for l in range(4):
                if label_dict[pdb2][:l+1] == label_dict[pdb][:l+1]:
                    tp_list[l] += 1
        tp_dict[pdb] = tp_list

    if save_file:
        with open(save_file,'wb') as p:
            pickle.dump(tp_dict,p)
    return tp_dict

def get_fraction_points(query_list, search_result, level, tp_dict, label_dict):
    """Compute fraction points (Sensitivity up to 1st false positive) for a given level.

    Args:
        query_list (list): List of query IDs.
        search_result (dict): {query: [target1, target2, ...]} search results.
        level (str): One of LEVELS.
        tp_dict (dict): {query: [TP_class, TP_archi, TP_topo, TP_homo]} total TP counts.
        label_dict (dict): {entry: [class, archi, topo, homo]} CATH labels.

    Returns:
        fraction (np.ndarray): Array of query fractions (0-1) for plotting.
        point_result (list): Sensitivity at each fraction point.
        scatter (list): Sensitivity of all queries (used for AUC computation).
    """
    level_idx = LEVELS.index(level)
    sensitivities = []

    print(f"Processing {level} level for {len(query_list)} queries...")

    for query in tqdm(query_list, desc=f'Processing {level}'):
        total_tp = tp_dict.get(query, [0]*4)[level_idx]
        if total_tp == 0:
            continue  # skip queries with no TP

        tp = 0
        for target in search_result.get(query, []):
            if target == query:
                continue  # skip self-hit
            target_label = label_dict.get(target)
            if target_label is None:
                continue  # skip unknown target
            if target_label[:level_idx+1] == label_dict[query][:level_idx+1]:
                tp += 1
            else:
                sensitivities.append(tp / total_tp)
                break
        else:
            # No FP encountered, use final TP/total TP
            sensitivities.append(tp / total_tp)

    # Sort sensitivities descending
    sensitivities.sort(reverse=True)
    scatter = sensitivities.copy()

    # Compute fraction points for plotting
    fraction = np.linspace(0, 1, 21)
    if len(sensitivities) == 0:
        point_result = [0.0] * len(fraction)
    else:
        cut_off_indices = np.floor(len(sensitivities) * fraction).astype(int)
        cut_off_indices[0] = 0
        cut_off_indices[-1] = len(sensitivities) - 1
        point_result = [sensitivities[i] for i in cut_off_indices]

    return fraction, point_result, scatter


def plot_fraction_curve(fraction, point_result, scatter, method_name, level, save_dir='figures'):
    """Plot fraction curve and calculate AUC"""
    os.makedirs(save_dir, exist_ok=True)
    plt.figure(figsize=(6,4), dpi=300)
    plt.plot(fraction, point_result, label=method_name, alpha=0.8)
    plt.xlabel('Query fraction')
    plt.ylabel('Sensitivity up to the 1st FPs')
    plt.title(f'{method_name} - {level}')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    save_path = os.path.join(save_dir,f'{method_name}_{level}.png')
    plt.savefig(save_path)
    plt.close()
    auc = np.trapz(np.array(scatter), np.arange(len(scatter))/len(scatter))
    return auc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--search_result', required=True, help='Search result TSV file')
    parser.add_argument('--cath_labels', default='cath-domain-list-S20-v4_4_0.txt', help='CATH labels file')
    parser.add_argument('--save_dir', default='figures', help='Directory to save figures')
    args = parser.parse_args()

    print("Reading search results...")
    search_result = read_search_result_tsv(args.search_result)

    print("Reading CATH labels...")
    label_dict = read_cath_labels(args.cath_labels)

    print("Preparing query list from CATH labels...")
    query_list = list(label_dict.keys())

    print("Computing TP counts...")
    tp_dict = compute_tp_count(label_dict)

    method_name = os.path.basename(args.search_result).replace('.tsv','')
    for level in LEVELS:
        fraction, point_result, scatter = get_fraction_points(query_list, search_result, level, tp_dict, label_dict)
        auc = plot_fraction_curve(fraction, point_result, scatter, method_name, level, save_dir=args.save_dir)
        print(f'{level} AUC: {auc:.4f}')

if __name__=="__main__":
    main()
