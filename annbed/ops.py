#!/usr/bin/env python3
"""
Author: Ken Chen (https://github.com/chenkenbio)
Date: 2025-01-13
"""

import argparse
import os
import sys
import math
import numpy as np
from tqdm import tqdm
import pandas as pd
from typing import Any, Dict, Iterable, List, Literal, Optional, Tuple, Union
from biock2 import random_string

TMPDIR = os.path.expanduser(os.environ.get("TMPDIR", "~/tmp"))
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)

## bedtools
def count_intersection(bed1: str, bed2: str, output: str, stranded: bool, *, min_overlap: int=None, max_dist: Union[int, Tuple[int, int]]=None):
    pass

def peak_enrichment_test_in_regions(bed1, bed2, regions, stranded: bool, prefix: str=None):
    if prefix is None:
        prefix = os.path.join(TMPDIR, random_string(10))
    else:
        prefix = prefix + "." + random_string(10)
    bedtools_intersect(bed1, regions, f"{prefix}.1.intersect", s=stranded)
    bedtools_intersect(bed2, regions, f"{prefix}.2.intersect", s=stranded)
    n_total_1 = count_peak_in_bed(bed1, stranded)
    n_total_2 = count_peak_in_bed(bed2, stranded)
    n_overlap_1 = count_peak_in_bed(f"{prefix}.1.intersect", stranded)
    n_overlap_2 = count_peak_in_bed(f"{prefix}.2.intersect", stranded)
    from scipy.stats import fisher_exact
    table = np.asarray([[n_overlap_1, n_overlap_2], [n_total_1 - n_overlap_1, n_total_2 - n_overlap_2]])
    oddsratio, pvalue = fisher_exact(table)
    oddsratio = float(oddsratio)
    pvalue = float(pvalue)
    table = pd.DataFrame(
        table,
        index=["in_regions", "not_in_regions"],
        columns=["peak1", "peak2"]
    )

    # clean up
    os.system(f"rm {prefix}.1.intersect {prefix}.2.intersect")

    # return (n_total_1, n_overlap_1), (n_total_2, n_overlap_2), pvalue, oddsratio
    from collections import namedtuple
    Overlap = namedtuple("Overlap", ["total", "overlap", "fraction"])
    Fisher_test = namedtuple("Fisher_test", ["pvalue", "oddsratio"])
    return Overlap(n_total_1, n_overlap_1, n_overlap_1/n_total_1), Overlap(n_total_2, n_overlap_2, n_overlap_2/n_total_2), Fisher_test(pvalue, oddsratio), table


def count_peak_in_bed(bed: str, stranded: bool, merge_overlap: bool=False):
    if merge_overlap:
        raise NotImplementedError("merge_overlap is not implemented yet")
    else:
        peaks = set()
        with open(bed) as f:
            for line in f:
                if stranded:
                    chrom, start, end, _, _, strand = line.strip().split("\t")[:6]
                    key = (chrom, start, end, strand)
                else:
                    chrom, start, end = line.strip().split("\t")[:3]
                    key = (chrom, start, end)
                    peaks.add(key)
    return len(peaks)
                

## Python wrapper for bedtools
def bedtools_intersect(bed1: str, bed2: str, output: str, *, s: bool):
    r"""
    Run bedtools intersect
    Args:
        bed1 (str): bed file 1
        bed2 (str): bed file 2
        output (str): output file
    """
    cmd = ' '.join([
        "bedtools intersect",
        "-a", bed1,
        "-b", bed2,
        "-wo",
        "-s" if s else "",
        f"> {output}"
    ])
    print("Run shell command: ", cmd, file=sys.stderr)
    os.system(cmd)

def bedtools_closest(
        bed1: str, bed2: str, 
        output: str, 
        *,
        s: bool, 
        d: bool=True,  # Report distance to closest feature by default
        D: Literal["a", "b"],
        d_cutoff: Union[int, Tuple[int, int]] = None,
    ):
    r"""
    Run bedtools closest
    Args:
        bed1 (str): bed file 1
        bed2 (str): bed file 2
        output (str): output file
        s (bool): Require same strandedness
        d (bool): Report distance to closest feature
        D (str): Report distance to closest feature in file D (a or b)
        d_cutoff (int, Tuple[int, int]): Report closest feature that is within this, or these, distance(s)
    """
    if D is not None:
        if d:
            d = False
        assert D in ["a", "b"], "D should be either 'a' or 'b'"
    if d_cutoff is not None:
        if isinstance(d_cutoff, int):
            assert D is None, "D should be None when d_cutoff is int"
            awk_cmd = f"awk '$NF < {d_cutoff}'"
        elif isinstance(d_cutoff, tuple):
            assert D is not None, "D should be either 'a' or 'b' when d_cutoff is tuple"
            assert len(d_cutoff) == 2, "d_cutoff should be a tuple of two integers"
            awk_cmd = f"awk '$NF > {d_cutoff[0]} && $NF < {d_cutoff[1]}'"
    else:
        awk_cmd = None
        
    cmd = [
        "bedtools closest",
        "-a", bed1,
        "-b", bed2,
        "-s" if s else "",
        "-d" if d else "",
        f"-D {D}", D if D is not None else "",
    ]
    if awk_cmd is not None:
        cmd.append(f"| {awk_cmd}")
            
    cmd = ' '.join(cmd) + f" > {output}"
    os.system(cmd)

   