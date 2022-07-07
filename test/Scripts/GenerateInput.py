#!/usr/bin/env python
# coding: utf8

import argparse
import logging
import os
import random
import sys
from typing import Any, Dict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
parser.add_argument(
    "--mode",
    choices=["big", "lengthmin", "qualsld", "qualtail", "infodust"],
    default="length",
    help="What kind of test file is to test",
)
# Output.
parser.add_argument("--outdir", help="Outdir")
parser.add_argument("--version", help="Show version and exit")

if len(sys.argv) < 1:
    print("No args provided")
    sys.exit()
else:
    args = parser.parse_args()

version = args.version

if version:
    print("Version " + VERSION)
    sys.exit(0)

mode = args.mode

outdir = os.getcwd()
if args.outdir:
    outdir = args.outdir

# Var.
INSTRUMENT = "M99999"
RUN = "100"
FLOWCELL = "000000000-BL3BP"
LANE = "1"
TILE = (1101, 2110)
COORDX = (2000, 29000)
COORDY = (1, 9999)
COORDYFORMAT = "{0:04d}"
QUAL = (0, 41)
QUALBAD = (1, 3)
QUALMEAN = (18, 21)
QUALGOOD = (25, 30)

DNA4 = ["A", "T", "C", "G"]
DNA5 = DNA4 + ["N"]


# FUNCTIONS
#  Global.
def generate_record_qual(generator, title, length, quals):
    seqs, qual_tmp = generate_seqs_quals(generator, length)
    rec = SeqRecord(Seq(seqs), id=title, name="", description="")
    rec.letter_annotations["phred_quality"] = quals
    return rec


def generate_record(title, seq, qual):
    rec = SeqRecord(Seq(seq), id=title, name="", description="")
    rec.letter_annotations["phred_quality"] = qual
    return rec


def generate_seq(generator, length):
    seqs = []
    for i in range(length):
        seqs.append(generator.choice(DNA4))
    return "".join(seqs)


def generate_qual(generator, length):
    quals = []
    for i in range(length):
        quals.append(generator.randrange(QUALGOOD[0], QUALGOOD[1]))
    return quals


def generate_seqs_quals(generator, length):
    return (generate_seq(generator, length), generate_qual(generator, length))


def generate_title(generator):
    tile = generator.randrange(TILE[0], TILE[1], 2)
    coord_x = generator.randrange(COORDX[0], COORDX[1], step=2)
    coord_y = generator.randrange(COORDY[0], COORDY[1], step=2)
    coord_y = COORDYFORMAT.format(coord_y)
    return ":".join([INSTRUMENT, RUN, FLOWCELL, LANE, str(tile), str(coord_x), coord_y])


#  Length.
def generate_record_length(generator, title, length):
    seqs, quals = generate_seqs_quals(generator, length)
    rec = SeqRecord(Seq(seqs), id=title, name="", description="")
    rec.letter_annotations["phred_quality"] = quals
    return rec


def generate_records_length(generator, length_forward, length_reverse):
    title = generate_title(generator)
    rec_forward = generate_record_length(generator, title, length_forward)
    rec_reverse = generate_record_length(generator, title, length_reverse)
    return rec_forward, rec_reverse


#  QualSld.
def generate_quality_qualsld(
    generator, length, qual_good, qual_bad, bin_start, bin_stop
):
    quals = []
    for i in range(bin_start):
        quals.append(generator.randrange(qual_good[0], qual_good[1]))
    for i in range(bin_start, bin_stop):
        quals.append(generator.randrange(qual_bad[0], qual_bad[1]))
    for i in range(bin_stop, length):
        quals.append(generator.randrange(qual_good[0], qual_good[1]))
    assert length == len(quals), "QualSld, len qual different of length"
    return quals


def generate_records_qualsld(
    generator, length, qual_good, qual_bad, bin_start, bin_stop, side
):
    title = generate_title(generator)
    rec_forward, rec_reverse = None, None
    quals = generate_quality_qualsld(
        generator, length, qual_good, qual_bad, bin_start, bin_stop
    )
    if side == "same":
        rec_forward = generate_record_qual(generator, title, length, quals)
        rec_reverse = generate_record_qual(generator, title, length, quals)
    elif side == "forward":
        rec_forward = generate_record_qual(generator, title, length, quals)
        rec_reverse = generate_record_length(generator, title, length)
    elif side == "reverse":
        rec_forward = generate_record_length(generator, title, length)
        rec_reverse = generate_record_qual(generator, title, length, quals)
    else:
        raise ValueError("QualSld, Indicated side")

    return rec_forward, rec_reverse


#  QualTail.
def generate_quality_qualtail(generator, length, qual_bad, bin_start):
    quals = []
    for i in range(bin_start):
        quals.append(generator.randrange(QUALGOOD[0], QUALGOOD[1]))
    for i in range(bin_start, length):
        quals.append(qual_bad)
    assert length == len(quals), "QualTail, len qual different of length"
    return quals


def generate_records_qualtail(generator, length, qual_bad, bin_start, side):
    title = generate_title(generator)
    rec_forward, rec_reverse = None, None
    quals = generate_quality_qualtail(generator, length, qual_bad, bin_start)
    if side == "same":
        rec_forward = generate_record_qual(generator, title, length, quals)
        rec_reverse = generate_record_qual(generator, title, length, quals)
    elif side == "forward":
        rec_forward = generate_record_qual(generator, title, length, quals)
        rec_reverse = generate_record_length(generator, title, length)
    elif side == "reverse":
        rec_forward = generate_record_length(generator, title, length)
        rec_reverse = generate_record_qual(generator, title, length, quals)
    else:
        raise ValueError("QualSld, Indicated side")
    return rec_forward, rec_reverse


#  Info.
def generate_records_info(generator, seqf, seqr):
    title = generate_title(generator)
    qualf = generate_qual(generator, len(seqf))
    qualr = generate_qual(generator, len(seqr))
    rec_forward = generate_record(title, seqf, qualf)
    rec_reverse = generate_record(title, seqr, qualr)
    return rec_forward, rec_reverse


records: Dict[Any, Any] = dict(forward=[], reverse=[])

if mode == "big":
    #  Number, length
    LENGTHSEQ = [
        (1500, 30),
        (2000, 40),
        (3000, 50),
        (1500, 75),
        (2000, 100),
    ]
    for ix in range(len(LENGTHSEQ)):
        params = random.choice(LENGTHSEQ)
        nb, length = params
        for x in range(nb):
            #  Rec.
            rec_f, rec_r = generate_records_length(random, length, length)
            records["forward"].append(rec_f)
            records["reverse"].append(rec_r)
        #  Update.
        LENGTHSEQ.remove(params)
elif mode == "lengthmin":
    # Length R1, Length R2
    LENGTHSEQ = [
        (20, 19),
        (20, 20),
        (20, 21),
        (49, 50),
        (51, 50),
        (50, 50),
        (100, 99),
        (100, 100),
        (101, 101),
    ]
    for ix in range(len(LENGTHSEQ)):
        lengths = random.choice(LENGTHSEQ)
        # Rec.
        rec_f, rec_r = generate_records_length(random, lengths[0], lengths[1])
        records["forward"].append(rec_f)
        records["reverse"].append(rec_r)

        # Update.
        LENGTHSEQ.remove(lengths)
elif mode == "qualsld":
    LENGTH = 100
    # Rec 1.
    rec_f, rec_r = generate_records_length(random, 5, 15)
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)
    # Rec 2.
    rec_f, rec_r = generate_records_length(random, 15, 5)
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)
    # Rec 3.
    rec_f, rec_r = generate_records_qualsld(
        random,
        LENGTH,
        QUALGOOD,
        QUALMEAN,
        75,
        80,
        "forward",
    )
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)
    # Rec 4.
    rec_f, rec_r = generate_records_qualsld(
        random,
        LENGTH,
        QUALGOOD,
        QUALMEAN,
        75,
        80,
        "reverse",
    )
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)
    # Rec 5.
    rec_f, rec_r = generate_records_qualsld(
        random,
        LENGTH,
        QUALGOOD,
        QUALMEAN,
        75,
        80,
        "same",
    )
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)
    # Rec 6.
    rec_f, rec_r = generate_records_qualsld(
        random,
        LENGTH,
        QUALGOOD,
        QUALMEAN,
        50,
        70,
        "forward",
    )
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)
    # Rec 7.
    rec_f, rec_r = generate_records_qualsld(
        random,
        LENGTH,
        QUALMEAN,
        QUALBAD,
        70,
        80,
        "reverse",
    )
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)
    # Rec 8.
    rec_f, rec_r = generate_records_qualsld(
        random,
        LENGTH,
        QUALMEAN,
        QUALBAD,
        3,
        6,
        "forward",
    )
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)
elif mode == "qualtail":
    LENGTH = 150
    # Rec 1.
    rec_f, rec_r = generate_records_length(random, 5, 15)
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)
    # Rec 2.
    rec_f, rec_r = generate_records_length(random, 15, 5)
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)
    # Rec 3.
    rec_f, rec_r = generate_records_qualtail(random, LENGTH, 5, 100, "forward")
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)
    # Rec 4.
    rec_f, rec_r = generate_records_qualtail(random, LENGTH, 5, 100, "reverse")
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)
    # Rec 5.
    rec_f, rec_r = generate_records_qualtail(random, LENGTH, 5, 100, "same")
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)
    # Rec 6.
    rec_f, rec_r = generate_records_qualtail(random, LENGTH, 10, 120, "forward")
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)
    # Rec 7.
    rec_f, rec_r = generate_records_qualsld(
        random, LENGTH, QUALMEAN, QUALBAD, 70, 80, "same"
    )
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)

elif mode == "infodust":
    # Rec 1. dust-1.06
    seqf = "GGCTGCAACAATCTTCCTTGTGTTAGCGTCATAAGAAATCAGGATCGTGT"
    seqr = "GGCTGCAACAATCTTCCTTGTGTTAGCGTCATAAGAAATCAGGATCGTGT"
    rec_f, rec_r = generate_records_info(random, seqf, seqr)
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)

    # Rec 2. dust-1.06 1.42
    seqf = "GGCTGCAACAATCTTCCTTGTGTTAGCGTCATAAGAAATCAGGATCGTGT"
    seqr = "TTCACATGGGAATGGCCGAAGAAAAGTAACGCGCCCGGGCGACCGTGTCG"
    rec_f, rec_r = generate_records_info(random, seqf, seqr)
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)

    # Rec 2. dust-1.98 2.40
    seqf = "TCCAACGCAACACTGGCCGCCTGAACACGTCGCACGGTATTGAGTAGAAA"
    seqf += "CCGCCTCCTTCCAACGTAACCCCCGTGATTCCACGACAAGAGCGCGCGT"
    seqf += "TCCCTTTACCGCTGACGGTTTGTTTTGTCTTTGTCAGTTTCGGAGTGAGAG"
    seqr = "TATGGTTGACTTGATAAATAAGATTTGCCATTGCAGCACACCGACCAGCC"
    seqr += "CCATTTTGCCCTGATCCGCAGGACCAGCACGGCAGCACTCCACCCCATT"
    seqr += "GTCCTATAAGCAGGTCCCGATATTGGCAACAGGCTGCATGTTTTGGGCGTG"
    rec_f, rec_r = generate_records_info(random, seqf, seqr)
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)

    # Rec 2. dust-2.97 4.45
    seqf = "TGTATCCCTGGTTTCAATTGTACGAAAAATTTTACGGGTAAAAGTTGATT"
    seqf += "TTTTCTTGTCGAAATCCAATTCACAATTTTGAATTTTTAAAGATTTCTC"
    seqf += "GCATCCTATAGGTAATGGATGGCGATTGGTAATATGACTTGCTAACCGT"
    seqf += "ATGAACGTGAAAAGGAATAAATATTTATCATTTAATGATAAATTATAACCTG"
    seqr = "GTATCTAAAAAGATCAACAATTTTAAGTATACCTAAAATACATATATAAA"
    seqr += "ACACAGTATAAACTATGAATCCAATGGACATTAATAAAAAATGATCTAA"
    seqr += "TTAATTACAAATCATACTTACTATACTAATTATACTCATTAATTCGTAA"
    seqr += "TACGAATTTATATATTATAAATACCAATATACTCTTTTAAACTAAATTA"
    seqr += "TATCTATAATAGATTACCTTAATGTTTATTTGAAGTCCTATCTATTCTT"
    seqr += "TTCCGATTATATATCTAGAATAACTCTATCTTATATATATTTTAACTAGCAAAC"
    rec_f, rec_r = generate_records_info(random, seqf, seqr)
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)

    # Rec 2. dust-6.03 2.40
    seqf = "AGTGAAAAGGGTTTTAAAAGGTCACTAATAATAATATTACACATATTTTA"
    seqf += "GTCACGTATAAAAAATATATATTACTATAATAGATTATAATATAATGAG"
    seqf += "TTAATATTATTACATAAAATAAACTTATAAATATGTAAACTATATAATA"
    seqf += "TAATTAGGACATATTTTATAATATTTTATTTAGCAAAAGCTCTTTAATG"
    seqf += "CAAACATACCAATTTACATAATAGCTATTTATATTTAAAAATTTTTATT"
    seqf += "ATTTATAATAATTTAATTAAAATACTTTGCAAAAAATAAATATCAATAAATGTT"
    seqr = "TATGGTTGACTTGATAAATAAGATTTGCCATTGCAGCACACCGACCAGCC"
    seqr += "CCATTTTGCCCTGATCCGCAGGACCAGCACGGCAGCACTCCACCCCATT"
    seqr += "GTCCTATAAGCAGGTCCCGATATTGGCAACAGGCTGCATGTTTTGGGCGTG"
    rec_f, rec_r = generate_records_info(random, seqf, seqr)
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)

    # Rec 2. dust-3.97 4.45
    seqf = "GCTAAATTAAAAGACTAGTAGCCTATATATACAGTTTTTGGTATTTAAAT"
    seqf += "TACAGCTTGTGTTTAGAAAGAAACTTCTTTTTACTTATATTAAAATTTT"
    seqf += "TTTTTATATTTAACCGCTCATATAATAATTATGTAATAATTTTTGGAAA"
    seqf += "AATGATATTATACCTAGTAACAAACACTTCAACTATTTTCTTATCCTATCTA"
    seqr = "GTATCTAAAAAGATCAACAATTTTAAGTATACCTAAAATACATATATAAA"
    seqr += "ACACAGTATAAACTATGAATCCAATGGACATTAATAAAAAATGATCTAA"
    seqr += "TTAATTACAAATCATACTTACTATACTAATTATACTCATTAATTCGTAA"
    seqr += "TACGAATTTATATATTATAAATACCAATATACTCTTTTAAACTAAATTA"
    seqr += "TATCTATAATAGATTACCTTAATGTTTATTTGAAGTCCTATCTATTCTT"
    seqr += "TTCCGATTATATATCTAGAATAACTCTATCTTATATATATTTTAACTAGCAAAC"
    rec_f, rec_r = generate_records_info(random, seqf, seqr)
    records["forward"].append(rec_f)
    records["reverse"].append(rec_r)
else:
    raise ValueError("Mode not implemented")

# Write
name_output = mode.upper()

# R1.
fforward = os.path.join(outdir, ".".join([name_output, "R1", "fastq"]))
with open(fforward, "w") as fod:
    SeqIO.write(records["forward"], fod, "fastq")
# R2.
freverse = os.path.join(outdir, ".".join([name_output, "R2", "fastq"]))
with open(freverse, "w") as fod:
    SeqIO.write(records["reverse"], fod, "fastq")

# Interleaved.
recs_inter = []
for rec_f, rec_r in zip(records["forward"], records["reverse"]):
    rec_f.id = rec_f.id + "\\1"
    rec_f.description = ""
    rec_f.name = ""
    recs_inter.append(rec_f)
    rec_r.id = rec_r.id + "\\2"
    rec_r.description = ""
    rec_r.name = ""
    recs_inter.append(rec_r)
finterleaved = os.path.join(outdir, ".".join([name_output, "Interleaved", "fastq"]))
with open(finterleaved, "w") as fod:
    SeqIO.write(recs_inter, fod, "fastq")
