# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Code was adjusted by Martin Steinegger to support alignment comparison


"""lDDT protein distance score."""
import json
import numpy as np
import argparse
from tqdm import tqdm

def read_ca_from_pdb(pdb_file):
    coords = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array([coords])

def lddt(target_points, query_points, aligned_pairs, cutoff=15., per_residue=False):
    assert len(target_points.shape) == 3
    assert target_points.shape[-1] == 3
    assert len(query_points.shape) == 3

    dmat_query = np.sqrt(np.sum((query_points[:, :, None] - query_points[:, None, :])**2, axis=-1))
    dists_to_score = ((dmat_query < cutoff).astype(np.float32) * (1. - np.eye(dmat_query.shape[1])))

    aligned_query_points = query_points[0, aligned_pairs[:, :, 0]]
    aligned_target_points = target_points[0, aligned_pairs[:, :, 1]]
    dmat_aligned_query = np.sqrt(np.sum((aligned_query_points[:, :, None] - aligned_query_points[:, None, :])**2, axis=-1))
    dmat_aligned_target = np.sqrt(np.sum((aligned_target_points[:, :, None] - aligned_target_points[:, None, :])**2, axis=-1))
    aligned_dists_to_score = ((dmat_aligned_query < cutoff).astype(np.float32) * (1. - np.eye(dmat_aligned_query.shape[1])))

    dist_l1 = np.abs(dmat_aligned_query - dmat_aligned_target)
    score = 0.25 * ((dist_l1 < 0.5).astype(np.float32) +
                    (dist_l1 < 1.0).astype(np.float32) +
                    (dist_l1 < 2.0).astype(np.float32) +
                    (dist_l1 < 4.0).astype(np.float32))
    score = score * (1. - np.eye(score.shape[1]))

    reduce_axes = (-1,) if per_residue else (-2, -1)
    norm = 1. / (np.sum(dists_to_score, axis=reduce_axes))
    norm_aligned = norm[0, aligned_pairs[:, :, 0]]
    score = norm_aligned * (np.sum(score * aligned_dists_to_score, axis=reduce_axes))
    return score

if __name__ == '__main__':
    BaseDir = '/data/liuyuan/workspace/Benchmark'
    parser = argparse.ArgumentParser()
    parser.add_argument('--aligned_json', type=str, help='aligned pairs JSON file')
    parser.add_argument('--query_folder', type=str, help='query PDB folder')
    parser.add_argument('--target_folder', type=str, help='target PDB folder')
    parser.add_argument('--output', type=str, help='output file')
    args = parser.parse_args()

    with open(args.aligned_json, 'r') as f:
        all_aligned = json.load(f)

    results = []
    for i, aln in tqdm(enumerate(all_aligned)):
        query_coords = read_ca_from_pdb(f"{args.query_folder}/{aln['query']}")
        target_coords = read_ca_from_pdb(f"{args.target_folder}/{aln['target']}")
        aligned_pairs = np.array([aln['pairs']], dtype=np.int32)

        per_residue_lddt = lddt(target_coords, query_coords, aligned_pairs, per_residue=True)
        avg_lddt = np.mean(per_residue_lddt[0])
        results.append(f"{aln['query']}\t{aln['target']}\t{avg_lddt:.4f}\n")

    if args.output:
        with open(args.output, 'w') as f:
            f.writelines(results)
