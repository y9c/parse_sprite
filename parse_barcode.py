#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 Ye Chang yech1990@gmail.com
# Distributed under terms of the MIT license.
#
# Created: 2022-03-13 23:00

"""parse SPRITE barcode.

R1: seq
R2: barcode
"""

import logging
import sys

import dnaio
import regex

from utils import get_min_edit_distance, overlap_pairend, reverse_complement

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(message)s")


def read_barcodes(barcode_file):
    barcode_info = {}
    with open(barcode_file, "r") as f:
        for line in f.readlines():
            n, b = line.strip().split("\t")
            g = n.split("-")[0]
            if g not in barcode_info:
                barcode_info[g] = {n: b}
            else:
                barcode_info[g][n] = b
    return barcode_info


def parse_barcode(record, barcode_info):
    r2 = record[1].sequence
    q2 = record[1].qualities
    m = regex.match(
        r"^(?P<TERM>[ACTGN]{8,11})(?:TGACTTG){s<=1}(?P<EVEN>[ACTGN]{14,18})(?:TGACAACT){s<=1}(?P<ODD>[ACTGN]{14,18})(?:TTGACTTG){s<=1}(?P<EVEN>[ACTGN]{14,18})(?:TGACAACT){s<=1}(?P<ODD>[ACTGN]{14,18})TT(?:GACTTGTCATGTCTTCCGAT){e<=2}CT(?P<DPM>[ACTGN]{7,9})[ACGTN]T(?P<SEQ>[ACTGN]{20,24})$",
        r2,
    )
    if m is None:
        return None

    barcode_list = [
        get_min_edit_distance(b, barcode_info[n])
        for n, bs in m.capturesdict().items()
        if n != "SEQ"
        for b in bs
    ]
    if "NA" in barcode_list:
        return None

    # cut read1 that overlap with read2
    adapter = (
        reverse_complement(m.group("SEQ"))
        + "AN"
        + reverse_complement(m.group("DPM"))
        + "AGATCGGAAGACATGACAAGTCAA"
    )
    term_len = len(m.group("TERM"))
    r1_cut = overlap_pairend(
        record[0].sequence,
        adapter,
        min_overlap=3,
        left_shift=term_len,
        right_shift=len(m.group("SEQ")) + 2,
    )

    new_name = record[0].name.split()[0] + " CB:" + ",".join(barcode_list)
    # edit read1
    record[0].name = new_name
    record[0].sequence = record[0].sequence[term_len:r1_cut]
    record[0].qualities = record[0].qualities[term_len:r1_cut]
    # edit read2
    record[1].name = new_name
    record[1].sequence = m.group("SEQ")
    record[1].qualities = q2[m.start("SEQ") : m.end("SEQ")]
    return record


if __name__ == "__main__":
    input_file_r1 = sys.argv[1]
    input_file_r2 = sys.argv[2]
    output_file_r1 = sys.argv[3]
    output_file_r2 = sys.argv[4]
    barcode_info = read_barcodes(sys.argv[5])

    fin = dnaio.open(input_file_r1, file2=input_file_r2, mode="r")
    fout = dnaio.open(
        output_file_r1, file2=output_file_r2, fileformat="fastq", mode="w"
    )
    n_input = 0
    n_failed = 0
    for record in fin:
        n_input += 1
        parsed_record = parse_barcode(record, barcode_info)
        if parsed_record is None:
            n_failed += 1
        else:
            fout.write(parsed_record[0], parsed_record[1])
            if n_input % 10_000 == 0:
                logging.debug(
                    f"{n_input:,} reads processed. {1-n_failed/n_input:.2%} with barcode."
                )
    fin.close()
    fout.close()

    logging.info(f"Number of input pairs: {n_input}")
    logging.info(f"Number of pairs without barcode: {n_failed}")
