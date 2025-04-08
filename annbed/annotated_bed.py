#!/usr/bin/env python3
"""
Author: Ken Chen (https://github.com/chenkenbio)
Date: 2024-12-31

Features:
- Annnotate BED to genomic regions/repeat elements
- plot distribution of genomic regions/repeat elements
- map BED to genes 
"""

import argparse
import os
import sys
import math
import glob
import numpy as np
from collections import defaultdict, OrderedDict
from typing import Any, Dict, Iterable, List, Literal, Optional, Tuple, Union, Callable
import time
from tqdm import tqdm
from multiprocessing import Process
import random, string

try:
    from .ops import bedtools_intersect, bedtools_closest, adjust_bed
except ImportError:
    from ops import bedtools_intersect, bedtools_closest, adjust_bed


"""
Desired funtionalities:
- Annotate BED to genomic regions/repeat elements
- plot numeric distribution of genomic regions/repeat elements
- enrichment analysis of genomic regions/repeat elements
- take bed subsets based on genomic regions/repeat elements

"""

HOME = os.environ.get("HOME", os.path.expanduser("~"))
REPEAT_DB = {
    "hg38_all": os.path.join(
        f"{HOME}/db/repeatmasker/hg38_repeatmasker.v0.bed.gz"
    ),
    "mm10_all": os.path.join(
        f"{HOME}/db/repeatmasker/mm10_repeatmasker.v0.bed.gz"
    )
}

def hash_string(s: str) -> int:
    return hash(s) % 10**16

TMPDIR = os.path.expanduser(os.environ.get("TMPDIR", "~/tmp"))
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)
def auto_open(file, mode='rt'):
    if file.endswith('.gz'):
        import gzip
        return gzip.open(file, mode)
    else:
        return open(file, mode)

def random_string(n: int=10) -> str:
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=n))

class AnnotatedBED(object):
    def __init__(
            self, 
            bed: str, 
            stranded: bool=True,
            key: Literal["loci", "name"]="loci",
            **kwargs):
        # bed file should have same columns in each row
        self.bed = bed
        self.std_bed = self.standardize_bed
        self.stranded = stranded
        assert key in ["loci", "name"], "key should be either 'loci' or 'name'"
        self.use_loci = key == "loci"

        self._prefix = os.path.join(TMPDIR, f"{os.path.basename(bed)}" + "." + random_string())
        assert len(glob.glob(f"{self._prefix}*")) == 0, "prefix collision with existing files, please try again"

        while len(glob.glob(f"{self._prefix}*")) > 0:
            self._prefix = os.path.join(TMPDIR, f"{os.path.basename(bed)}" + "." + random_string())
        self._peak_name2id = dict()
        self._peak_id2name = list()
        self.n_cols, self.header = self._prepare_bed(bed)
    
    def annotate(self):
        annotate_bed(
            self.bed, 
            gene_db=GR_DB["hg38"], 
            repeat_db=REPEAT_DB["hg38"],
            repeat_stranded=False,
            output=self._prefix + ".annotated.tsv",
            ignore_ids=False
        )
    
    @property
    def standardize_bed(self) -> str:
        # return BED-6 format
        std_bed = f"{self._prefix}.standard.bed"
        if not os.path.exists(std_bed):
            with auto_open(self.bed, 'rt') as infile, open(std_bed, 'wt') as outfile:
                for l in infile:
                    fields = l.rstrip().split('\t')
                    peak = get_peak_name(fields, self.use_loci, self.stranded)
                    peak_id = self.peak2id(peak)
                    fields[3] = str(peak_id)
                    if self.stranded:
                        fields = fields[:6]
                    else:
                        fields = fields[:4]
                    outfile.write("\t".join(fields) + "\n")
        return std_bed
    
    def _prepare_bed(self, bed) -> Tuple[int, List[str]]:
        r"""
        Check the number of columns in the bed file.
        Args:
            bed (str): bed file
            nlines (int): number of lines to check
        Returns:
            n_cols (int): number of columns
            header (List[str]): header of the bed file
        """
        self.header = None
        self.n_cols = 0
        cnt = 0
        with auto_open(bed, 'rt') as infile:
            for l in infile:
                fields = l.rstrip().split('\t')
                if l.startswith("#"):
                    if cnt == 0:
                        self.header = fields
                        self.n_cols = len(fields)
                    continue
                if self.n_cols > 0:
                    assert len(fields) == self.n_cols, f"Number of columns in line {cnt} is not consistent with header"
                else:
                    self.n_cols = len(fields)
                cnt += 1
                if self.use_loci:
                    if self.stranded:
                        peak = f"{fields[0]}_{fields[1]}_{fields[2]}_{fields[5]}"
                    else:
                        peak = f"{fields[0]}_{fields[1]}_{fields[2]}"
                else:
                    peak = fields[3]
                _ = self.peak2id(peak)

        if self.header is None:
            self.header = [f"col_{i}" for i in range(self.n_cols)]
        self.header[0] = "chrom"
        self.header[1] = "start"
        self.header[2] = "end"
        if self.n_cols > 5:
            self.header[5] = "strand"
        return self.n_cols, self.header

    def id2peak(self, id: int):
        assert isinstance(id, int), "id should be an integer"
        assert id < len(self._peak_id2name), "id should be less than the number of peaks"
        return self._peak_id2name[id]

    def peak2id(self, name):
        if name not in self._peak_name2id:
            self._peak_name2id[name] = len(self._peak_id2name)
            self._peak_id2name.append(name)
        return self._peak_name2id[name]
 
    def cleanup(self):
        # remove tmp files with prefix
        os.system(f"rm -f {self._prefix}*")

def get_peak_name(fields: List[str], use_loci: bool, stranded: bool) -> str:  
    if use_loci:
        if stranded:
            peak = f"{fields[0]}_{fields[1]}_{fields[2]}_{fields[5]}"
        else:
            peak = f"{fields[0]}_{fields[1]}_{fields[2]}"
    else:
        peak = fields[3]
    return peak

class BED2Repeats(object):
    def __init__(self, bed: str, repeat_bed: str, prefix: str=None, key: Literal["loci", "name"]="loci"):
        assert len(glob(prefix + "*")) == 0, "prefix should not exist"
        self.bed = bed
        self.repeat_bed = repeat_bed
        assert os.path.exists(bed), f"{bed} does not exist"
        assert os.path.exists(repeat_bed), f"{repeat_bed} does not exist"
        assert key in ["loci", "name"], "key should be either 'loci' or 'name'"

        self.use_loci = key == "loci"
        self.repeat_info = dict()
        self._prefix = prefix
    
    def _annotate_bed(self, bed: str, repeat_db: str) -> str:
        h = hash(os.path.realpath(bed) + os.path.realpath(repeat_db)) % 10**16
        output = f"{self._prefix}.to_repeat.{h}.intersect"
        bedtools_intersect(
            bed1=bed,
            bed2=repeat_db,
            output=output,
        )
        return output

    @staticmethod
    def _extract_rep_names(name: str) -> Tuple[str, str, str]:
        rep_type = name
        rep_family = name.split(',')[0]
        rep_class = rep_family.split('/')[0]
        return rep_type, rep_family, rep_class
    
    def parse_repeat_intersect(self, intersect: str):
        repclass = defaultdict(int)
        with open(intersect) as infile:
            for l in infile:
                if l.startswith("#") or l.startswith("track"):
                    continue
                fields = l.rstrip().split('\t')
                peak = get_peak_name(fields, self.use_loci, stranded=False)
                rep_id, rep_name = fields[-4], fields[-3]
                reptype, repfamily, repclass = self._extract_rep_names(rep_name)

# def annotate_bed(
#         in_bed, db_bed, *, 
#         stranded: bool, 
#         up: int=0, down: int=0, 
#         chromsize=None, 
#         prefix: str=None,
#         ):
#     if prefix is None:
#         prefix = os.path.join(TMPDIR, random_string(10))
#     else:
#         prefix = prefix + "." + random_string(10)
#     assert len(glob.glob(f"{prefix}*")) == 0, "prefix collision with existing files, please try again"
#     cmd = list()
#     if up > 0 or down > 0:
#         assert chromsize is not None and os.path.exists(chromsize), "chromsize should be provided"
#         cmd.append(f"bedtools slop -i {in_bed} -g {chromsize} -l {up} -r {down}")
#         if stranded:
#             cmd.append("-s |")
#         else:
#             cmd.append("|")
#     else:
#         cmd.append(f"cat {in_bed} | ")
#     cmd.append(f"bedtools intersect -wo -a stdin -b {db_bed}")
#     if stranded:
#         cmd.append("-s")
#     cmd.append(f"> {prefix}.intersect")
#     cmd = " ".join(cmd)
#     os.system(cmd)
#     return prefix + ".intersect"
    
        
import pandas as pd
class BED2GenomicRegions(object):
    def __init__(
            self, 
            bed: str, 
            genomic_regions: str, 
            prefix: str=None
        ):
        assert len(glob(prefix + "*")) == 0, "prefix should not exist"

def get_bed_info(bed: str, max_line: int=1000) -> Tuple[int, List[str]]:
    ncols = None
    header = None
    with auto_open(bed, 'rt') as infile:
        for i, l in enumerate(infile):
            if max_line is not None and i > max_line:
                break
            if l.startswith("#") or l.startswith("track"):
                if l.startswith("#"):
                    header = l.lstrip("#").rstrip().split("\t")
                continue
            fields = l.rstrip().split("\t")
            if ncols is None:
                ncols = len(fields)
            else:
                assert len(fields) == ncols, f"Number of columns in line {i} is not consistent with header"
    if header and len(header) != ncols:
        header = None
    return ncols, header

from biock2 import interval_distance

def _pair_peak_to_feature(bed, anno_bed, key_value: str, stranded: bool, up: int, down: int, *, key_extractor, type_extractor):
    assert key_value in ["loci", "name"], "key should be either 'loci' or 'name'"
    use_loci = key_value == "loci"
    tmp_tsv = os.path.join(TMPDIR, f"{os.path.basename(bed)}_annotated" + "." + random_string(10) + ".tsv")
    cmd = list()
    cmd.append(f"bedtools window -l {up} -r {down}")
    if stranded:
        cmd.append("-sw -sm")

    cmd.append(f"-b {bed} -a {anno_bed}")
    cmd.append(" > " + tmp_tsv)
    cmd = " ".join(cmd)
    # print("Running command:", cmd, file=sys.stderr)
    os.system(cmd)
    bed_cols, _ = get_bed_info(anno_bed)
    n_cols = None
    pairs = dict()
    with open(tmp_tsv) as infile:
        for l in infile:
            if l.startswith("#") or l.startswith("track"):
                continue
            fields = l.rstrip('\n').split("\t")
            if n_cols is None:
                n_cols = len(fields)
            else:
                assert len(fields) == n_cols, "Number of columns is not consistent"
            peak_key = get_peak_name(fields[bed_cols:], use_loci, stranded=True)
            feature_key = key_extractor(fields)
            feature_type = type_extractor(fields)
            # calculate  distance
            f_start, f_end = int(fields[1]), int(fields[2])
            p_start, p_end = int(fields[bed_cols + 1]), int(fields[bed_cols + 2])
            d = interval_distance((p_start, p_end), (f_start, f_end), ref="right")
            if not stranded:
                d = abs(d)
            elif fields[5] == "-":
                d = -d
            if peak_key not in pairs:
                pairs[peak_key] = list()
            pairs[peak_key].append((feature_key, feature_type, d))
    return pairs


HOME = os.environ.get("HOME", os.path.expanduser("~"))
HG38_GENE_PATH = f"{HOME}/db/gencode/GRCh38/v46/regions"
GR_DB = {
    "hg38": {
        "gene": [os.path.join(HG38_GENE_PATH, "gencode.v46.gene.bed.gz"), (0, 0), True], # path, (up, down), stranded
        "TSS": [os.path.join(HG38_GENE_PATH, "gencode.v46.TSS.bed.gz"), (3000, 0), True],
        "anti-TSS": [os.path.join(HG38_GENE_PATH, "gencode.v46.anti-TSS.bed.gz"), (0, 3000), True],
        "CDS": [os.path.join(HG38_GENE_PATH, "gencode.v46.CDS.bed.gz"), (0, 0), True],
        "exon": [os.path.join(HG38_GENE_PATH, "gencode.v46.exon.bed.gz"), (0, 0), True],
        "UTR3": [os.path.join(HG38_GENE_PATH, "gencode.v46.UTR3.bed.gz"), (0, 0), True],
        "UTR5": [os.path.join(HG38_GENE_PATH, "gencode.v46.UTR5.bed.gz"), (0, 0), True],
        "intron": [os.path.join(HG38_GENE_PATH, "gencode.v46.intron.bed.gz"), (0, 0), True],
        "TES": [os.path.join(HG38_GENE_PATH, "gencode.v46.TES.bed.gz"), (0, 3000), True],
    },
    "mm10": {
        "gene": ["/home/kenchen/db/gencode/GRCm38/vM25/gencode.vM25.gene.bed.gz", (0, 0), True],
        "TSS": ["/home/kenchen/db/gencode/GRCm38/vM25/gencode.vM25.TSS.bed.gz", (3000, 0), True],
        "anti-TSS": ["/home/kenchen/db/gencode/GRCm38/vM25/gencode.vM25.anti-TSS.bed.gz", (0, 3000), True],
        "CDS": ["/home/kenchen/db/gencode/GRCm38/vM25/gencode.vM25.CDS.bed.gz", (0, 0), True],
        "exon": ["/home/kenchen/db/gencode/GRCm38/vM25/gencode.vM25.exon.bed.gz", (0, 0), True],
        "UTR3": ["/home/kenchen/db/gencode/GRCm38/vM25/gencode.vM25.UTR3.bed.gz", (0, 0), True],
        "UTR5": ["/home/kenchen/db/gencode/GRCm38/vM25/gencode.vM25.UTR5.bed.gz", (0, 0), True],
        "intron": ["/home/kenchen/db/gencode/GRCm38/vM25/gencode.vM25.intron.bed.gz", (0, 0), True],
        "TES": ["/home/kenchen/db/gencode/GRCm38/vM25/gencode.vM25.TES.bed.gz", (0, 3000), True],
    }
}
REPEAT_DB = {
    "hg38": {
        "repeat": [
            os.path.join(f"{HOME}/db/repeatmasker/hg38_repeatmasker.v0.bed.gz"),
            (0, 0),
            False
        ]
    },
    "mm10": {
        "repeat": [
            os.path.join(f"{HOME}/db/repeatmasker/mm10_repeatmasker.v0.bed.gz"),
            (0, 0),
            False
            ]
        }
}


def split_repeat_type_name(name: str) -> Tuple[str, str, str]:
    """
    Return:
    - repeat type
    - repeat family
    - repeat class
    """
    rep_type = name
    rep_family = name.split(',')[1]
    rep_class = rep_family.split('/')[0]
    return rep_type, rep_family, rep_class

def annotate_bed(
        bed, 
        *, 
        gene_db: Dict[str, Tuple[str, Tuple[int, int]]]=GR_DB["hg38"], 
        repeat_db: Dict[str, Tuple[str, Tuple[int, int]]]=REPEAT_DB["hg38"],
        repeat_stranded: bool=False,
        output, 
        ignore_ids: bool=False,
        in_memory: bool=True):
    db = dict()
    for k, v in gene_db.items():
        db[k] = v
    for k, v in repeat_db.items():
        if repeat_stranded:
            db[k] = (v[0], v[1], True)
        else:
            db[k] = (v[0], v[1], False)
    links = OrderedDict()
    is_stranded = False
    for feature_type, (db_bed, (up, down), stranded) in db.items():
        if stranded:
            is_stranded = True
        if feature_type in {"gene"}:
            key_extractor = lambda x: (feature_type + "|" + '/'.join(x[3].split('|')[0:3]))
            type_extractor = lambda x: x[4]
        elif feature_type in {"TSS", "anti-TSS", "CDS", "exon", "UTR3", "UTR5", "intron", "TES"}:
            key_extractor = lambda x: (feature_type + "|" + '/'.join(x[3].split('|')[0:3]))
            type_extractor = lambda x: x[4]
        elif feature_type == "repeat":
            key_extractor = lambda x: x[3]
            type_extractor = lambda x: split_repeat_type_name(x[4])
        else:
            raise ValueError(f"key {feature_type} not recognized")

        pairs = _pair_peak_to_feature(
            bed, db_bed, key_value="loci", stranded=stranded, up=up, down=down,
            key_extractor=key_extractor, type_extractor=type_extractor
        )
        links[feature_type] = pairs
    
    n_cols, header = get_bed_info(bed)
    if header is None:
        header = [f"col_{i}" for i in range(n_cols)]
        header[0] = "chrom"
        header[1] = "start"
        header[2] = "end"
        if is_stranded:
            header[5] = "strand"

    for feature_type in links.keys():
        if feature_type == "repeat":
            header.append("repeat_class\trepeat_family\trepeat_type")
            if not ignore_ids:
                header.append(f"{feature_type}_ids")
        elif feature_type == "gene":
            header.append("gene")
            if not ignore_ids:
                header.append(f"{feature_type}_ids")
        else:
            header.append(f"{feature_type}")
            if not ignore_ids:
                header.append(f"{feature_type}_ids")
    
    print(bed)
    with auto_open(bed, 'rt') as infile, auto_open(output, 'wt') as outfile:
        print("\t".join(header), file=outfile)
        for l in infile:
            if l.startswith("#") or l.startswith("track"):
                continue
            fields = l.rstrip().split("\t")
            peak_key = get_peak_name(fields, use_loci=True, stranded=is_stranded)
            for group, pairs in links.items():
                feature_keys, feature_types = set(), set()
                if peak_key in pairs:
                    for feature_key, feature_type, d in pairs[peak_key]:
                        feature_keys.add(f"{feature_key}|{d}")
                        feature_types.add(feature_type)
                    feature_keys = ";".join(sorted(feature_keys))
                    if group == "repeat":
                        rep_types, rep_families, rep_classes = set(), set(), set()
                        for (rtype, rfamily, rclass) in feature_types:
                            rep_types.add(rtype)
                            rep_families.add(rfamily)
                            rep_classes.add(rclass)
                        rep_types = ";".join(sorted(rep_types))
                        rep_families = ";".join(sorted(rep_families))
                        rep_classes = ";".join(sorted(rep_classes))
                        feature_types = f"{rep_classes}\t{rep_families}\t{rep_types}"
                    else:
                        feature_types = ";".join(sorted(feature_types))
                else:
                    if group == "repeat":
                        feature_keys = "."
                        feature_types = ".\t.\t."
                    else:
                        feature_keys, feature_types = ".", "."
                fields.append(feature_types)
                if not ignore_ids:
                    fields.append(feature_keys)
                # fields.append(feature_keys)
            print("\t".join(fields), file=outfile)
    return links

def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("bed", help="Input bed file")
    p.add_argument("--output", "-o", help="Output file", required=True)
    p.add_argument("--gene-db", help="Gene database", default="hg38", choices=["hg38", "mm10"], required=True)
    p.add_argument("--repeat-db", help="Repeat database", default="hg38", choices=["hg38", "mm10"], required=True)
    p.add_argument("--repeat-stranded", action="store_true", help="Repeat database is stranded")
    p.add_argument("--ignore-ids", action="store_true", help="Ignore IDs")
    p.add_argument("--region-from", choices=["bed", "index"], default="bed", help="Genomic region to annotate")
    p.add_argument("--header", choices=("Y", "N"), help="Header in the input bed file")
    args = p.parse_args()
    if args.region_from == "index":
        assert args.header is not None, "Header should be provided if region_from is index"
        tmp_bed = os.path.join(TMPDIR, f"{os.path.basename(args.bed)}" + "." + random_string())
        with open(args.bed, 'rt') as infile, open(tmp_bed, 'wt') as outfile:
            for nr, l in enumerate(infile):
                fields = l.rstrip().split("\t")
                if nr == 0 and args.header == "Y":
                    print("#chrom\tstart\tend\tname\tname2\tstrand\t" + '\t'.join(fields[1:]), file=outfile)
                else:
                    chrom, start, end, strand = fields[0].split("_")
                    name = fields[0]
                    print(f"{chrom}\t{start}\t{end}\t{name}\t{name}\t{strand}\t" + '\t'.join(fields[1:]), file=outfile)
        print(f"Using {tmp_bed} as input", file=sys.stderr)
        args.bed = tmp_bed
    else:
        tmp_bed = None

    annotate_bed(
        args.bed, 
        gene_db=GR_DB[args.gene_db], 
        repeat_db=REPEAT_DB[args.repeat_db],
        repeat_stranded=args.repeat_stranded,
        output=args.output,
        ignore_ids=args.ignore_ids
    )
    # if tmp_bed:
    #     os.system(f"rm -f {tmp_bed}")

if __name__ == "__main__":
    main()
