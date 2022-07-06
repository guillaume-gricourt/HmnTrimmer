#!/usr/bin/env python
# coding: utf8

import argparse
import logging
import os
import random
import sys

from Bio import SeqIO

VERSION = "0.5.0"

# Init.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%d-%m-%Y %H:%M",
)
random.seed(0)

description = "Compute test file for HmnTrimmer"
parser = argparse.ArgumentParser(description=description)
parser.add_argument("--input-forward", required=True, help="Input forward")
parser.add_argument("--input-reverse", required=True, help="Input forward")
parser.add_argument("--indices", nargs="+", help="Indices to keep, 1-based")
# Output.
parser.add_argument("--outdir", required=True, help="Outdir")
parser.add_argument("--basename", required=True, help="Basename")

if len(sys.argv) < 1:
    print("No args provided")
    sys.exit()
else:
    args = parser.parse_args()

input_forward = args.input_forward
input_reverse = args.input_reverse
indices = args.indices

basename = args.basename
outdir = args.outdir

recs_f = [x for x in SeqIO.parse(input_forward, "fastq")]
recs_r = [x for x in SeqIO.parse(input_reverse, "fastq")]

recs_f_new, recs_r_new = [], []
for i in indices:
    i = int(i) - 1
    recs_f_new.append(recs_f[i])
    recs_r_new.append(recs_r[i])

recs_f = recs_f_new
recs_r = recs_r_new

recs_inter = []
for rec_f, rec_r in zip(recs_f, recs_r):
    rec_f.id = rec_f.id + "\\1"
    rec_f.description = ""
    rec_f.name = ""
    recs_inter.append(rec_f)
    rec_r.id = rec_r.id + "\\2"
    rec_r.description = ""
    rec_r.name = ""
    recs_inter.append(rec_r)

with open(os.path.join(outdir, basename + ".R1.fastq"), "w") as fod:
    SeqIO.write(recs_f, fod, "fastq")
with open(os.path.join(outdir, basename + ".R2.fastq"), "w") as fod:
    SeqIO.write(recs_r, fod, "fastq")
with open(os.path.join(outdir, basename + ".Interleaved.fastq"), "w") as fod:
    SeqIO.write(recs_inter, fod, "fastq")
