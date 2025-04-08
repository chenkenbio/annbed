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
import warnings
_bioframe_inabled = False
try:
    import bioframe as bf
    _bioframe_inabled = True
except ImportError:
    warnings.warn("bioframe is not installed, use pandas instead. Some functionalities may not work.")
from typing import Any, Dict, Iterable, List, Literal, Optional, Tuple, Union
from biock2 import random_string, auto_open

TMPDIR = os.path.expanduser(os.environ.get("TMPDIR", "~/tmp"))
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)

## bedtools
def count_intersection(bed1: str, bed2: str, output: str, stranded: bool, *, min_overlap: int=None, max_dist: Union[int, Tuple[int, int]]=None):
    pass

def peak_enrichment_test_in_regions(bed1, bed2, regions, stranded: bool, max_dist: int=0, prefix: str=None, suppress_warning: bool=False):
    if prefix is None:
        prefix = os.path.join(TMPDIR, random_string(10))
    else:
        prefix = prefix + "." + random_string(10)
    if max_dist > 0:
        bedtools_window(bed1, regions, f"{prefix}.1.intersect", up=max_dist, down=max_dist, s=stranded, suppress_warning=suppress_warning)
        bedtools_window(bed2, regions, f"{prefix}.2.intersect", up=max_dist, down=max_dist, s=stranded, suppress_warning=suppress_warning)
    else:
        bedtools_intersect(bed1, regions, f"{prefix}.1.intersect", s=stranded, suppress_warning=suppress_warning)
        bedtools_intersect(bed2, regions, f"{prefix}.2.intersect", s=stranded, suppress_warning=suppress_warning)
    print(prefix)
    n_total_1 = count_peaks(bed1, stranded)
    n_total_2 = count_peaks(bed2, stranded)
    print(bed1, n_total_1, bed2, n_total_2)
    n_overlap_1 = count_peaks(f"{prefix}.1.intersect", stranded)
    n_overlap_2 = count_peaks(f"{prefix}.2.intersect", stranded)
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
    # os.system(f"rm {prefix}.1.intersect {prefix}.2.intersect")

    # return (n_total_1, n_overlap_1), (n_total_2, n_overlap_2), pvalue, oddsratio
    from collections import namedtuple
    Overlap = namedtuple("Overlap", ["total", "overlap", "fraction"])
    Fisher_test = namedtuple("Fisher_test", ["pvalue", "oddsratio"])
    return Overlap(n_total_1, n_overlap_1, n_overlap_1/n_total_1), Overlap(n_total_2, n_overlap_2, n_overlap_2/n_total_2), Fisher_test(pvalue, oddsratio), table


def count_peaks(bed: str, stranded: bool, merge_overlap: bool=False, start_column: int=0):
    r"""
    Count number of peaks in a bed file
    """
    if merge_overlap:
        raise NotImplementedError("merge_overlap is not implemented yet")
    else:
        peaks = set()
        with open(bed) as f:
            for line in f:
                if line.startswith("#") or line.startswith("track"):
                    continue
                if stranded:
                    chrom, start, end, _, _, strand = line.strip().split("\t")[start_column:start_column+6]
                    key = (chrom, start, end, strand)
                else:
                    chrom, start, end = line.strip().split("\t")[start_column:start_column+3]
                    key = (chrom, start, end)
                peaks.add(key)
    return len(peaks)

   

def peak_to_gene(bed, gene_bed, *, output, stranded, genome, d_cutoff=Tuple[int, int], peak_naming: Literal["loci", "name"]="loci"):
    pass

def adjust_bed(
        bed, 
        shift1: int, # minus: upstream, plus: downstream
        shift2: int, 
        output: str,
        anchor1: Literal["start", "end", "center"]="start", 
        anchor2: Literal["start", "end", "center"]="end",
        chromsize: Optional[str]=None,
        in_memory: bool=False,
    ):
    if chromsize is not None:
        chromsizes = pd.read_csv(chromsize, sep="\t", header=None, index_col=0)[1].to_dict()
    else:
        chromsizes = dict()
    
    def _adjust(start, end, shift, anchor, chromsize):
        if anchor == "start":
            pos = start
        elif anchor == "end":
            pos = end
        elif anchor == "center":
            pos = (start + end) // 2
        else:
            raise ValueError(f"Unknown anchor: {anchor}")
        pos += shift
        if pos < 0:
            pos = 0
        if chromsize is not None:
            if pos > chromsize:
                pos = chromsize
        return pos
    
    with auto_open(bed, 'rt') as infile, auto_open(output, 'wt') as outfile:
        for l in infile:
            if l.startswith("#") or l.startswith("track"):
                print(l, end="", file=outfile)
                continue
            fields = l.rstrip().split("\t")
            start, end = int(fields[1]), int(fields[2])
            chrom = fields[0]
            left = _adjust(start, end, shift1, anchor1, chromsizes.get(chrom, None))
            right = _adjust(start, end, shift2, anchor2, chromsizes.get(chrom, None))
            fields[1] = str(left)
            fields[2] = str(right)
            print("\t".join(fields), file=outfile)
    if in_memory:
        if _bioframe_inabled:
            return bf.read_table(output, schema="bed")
        else:
            return pd.read_csv(output, sep="\t", header=None)
    else:
        return output



## Python wrapper for bedtools
def bedtools_intersect(bed1: str, bed2: str, output: str, *, s: bool, fa: float=None, fb: float=None, header: bool=False, suppress_warning: bool=False):
    r"""
    Run bedtools intersect
    Args:
        bed1 (str): bed file 1
        bed2 (str): bed file 2
        output (str): output file
    """
    cmd = [
        "bedtools intersect",
        "-wo",
        "-s" if s else "",
    ]
    if header:
        cmd.append("-header")
    if fa is not None:
        cmd.append(f"-f {fa}")
    if fb is not None:
        cmd.append(f"-f {fb}")
    cmd.extend(["-a", bed1, "-b", bed2])
    if suppress_warning:
        cmd.append("2> /dev/null")
    cmd = ' '.join(cmd) + f" > {output}"
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
        suppress_warning: bool=False,
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
    if suppress_warning:
        cmd.append("2> /dev/null")
    if awk_cmd is not None:
        cmd.append(f"| {awk_cmd}")
            
    cmd = ' '.join(cmd) + f" > {output}"
    os.system(cmd)


def bedtools_window(
        bed1: str, bed2: str,
        output: str,
        *,
        up: int=0, down: int=0,
        s: bool=False,
        header: bool=False,
        to_dataframe: bool=False,
        suppress_warning: bool=False,
    ):
    r"""

    """
    cmd = [
        "bedtools window",
        "-a", bed1,
        "-b", bed2,
        "-l", str(up),
        "-r", str(down),
        "-sw -sm" if s else "",
        "-header" if header else ""
    ]
    if suppress_warning:
        cmd.append("2> /dev/null")
    cmd = ' '.join(cmd) + f" > {output}"
    print("Run shell command: ", cmd, file=sys.stderr)
    os.system(cmd)
    assert os.path.exists(output)
    return output


