#!/usr/bin/env python

import os
import argparse


def main():
    parser = argparse.ArgumentParser(description="", usage="")
    parser.add_argument("vcf", metavar="<snvs.vcf>", type=str, help="vcf file with SNVs.")
    parser.add_argument("fai", metavar="<contigs.fasta.fai>", type=str, help="fai file with all reference sequences.")
    parser.add_argument("-c", metavar="STR", type=str, help="liftover contig")
    parser.add_argument("-f", metavar="INT", type=int, help="liftover offset")

    args = parser.parse_args()

    vcf_file = os.path.abspath(args.vcf)
    fai_file = os.path.abspath(args.fai)

    # Load the ctgs
    ctg_lens = dict()
    with open(fai_file, "r") as f:
        for line in f:
            h, l, x, y, z = line.split("\t")
            ctg_lens[h] = int(l)

    # Get the liftover contig
    lo_ctg = args.c
    if lo_ctg not in ctg_lens:
        raise ValueError("Liftover contig must be in the provided fai file")

    # Get the liftover positional offset
    lo_offset = args.f
    if lo_offset < 0:
        raise ValueError("Liftover offset must be >= 0")

    # Iterate over the VCF file
    with open(vcf_file, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("#"):  # header lines
                if line.startswith("##contig="):
                    for i in sorted(list(ctg_lens.keys())):
                        print("##contig=<ID={}>".format(i))
                else:
                    print(line)

            else:  # variant lines
                fields = line.split("\t")
                fields[0] = lo_ctg
                fields[1] = str(int(fields[1]) + lo_offset)
                print("\t".join(fields))