#!/usr/bin/env python

import os
import sys
import argparse

import pysam

from patchpolish_utilities.utilities import log, run_oe


def make_coord_name(h, s, e):
    return "{}_{}_{}".format(h, s, e)


def main():
    parser = argparse.ArgumentParser(description='Reference-guided misassembly correction',
                                     usage="ragtag.py correct <reference.fa> <query.fa>")

    parser.add_argument("reference", metavar="<patched_chromosomes.fa>", nargs='?', default="", type=str, help="reference fasta file containing raw patches. must not be gzipped.")
    parser.add_argument("reads", metavar="<ont_reads.fa>", nargs='?', default="", type=str, help="reads. must not be gzipped.")
    parser.add_argument("-o", metavar="STR", type=str, default="patchpolish_output", help="output directory [patchpolish_output]")

    args = parser.parse_args()

    if not args.reference or not args.reads:
        parser.print_help()
        sys.exit()

    reference_file = os.path.abspath(args.reference)
    reads_file = os.path.abspath(args.reads)

    # Setup the output directory
    output_path = args.o.replace("/", "").replace(".", "")
    cwd = os.getcwd()
    output_path = cwd + "/" + output_path + "/"
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    supp_reads = [
        "06210b04-c852-422b-b1e6-adeb3da63ab1",
        "223bfa25-d77a-4369-b960-9189d76fe220",
        "242b72a0-23b8-4a09-9feb-02a3fd885c1c",
        "28f594fc-2c81-45d4-a737-4be17277d509",
        "2a6b7a57-1a1e-461a-93de-4bbfec509d54",
        "34408519-ed07-4913-800a-40ee56af74ef",
        "49492f04-e697-46ca-8958-b1a30783d8db",
        "65f8fdc2-7ff0-472f-a978-47c2f2b12881",
        "73bc2a53-8e0a-4862-9c02-73030c90bac3",
        "874a8282-2dca-4698-b188-c8e6a78d419a",
        "8c51c819-20c5-46ce-9f9f-75ecdd6bb106",
        "959b8cb8-a3a7-4265-b8c1-426088657fc0",
        "d396acbe-63c4-4ddf-96db-d19ed77f9147",
        "da3e36bc-f8d3-41cf-abcd-b0e84e2e736b",
        "e1856fd3-dbee-4103-ae2a-3647d93404d3",
        "e2106675-aabd-4768-be89-eb898e46094a"
    ]

    patch_chr = "chr17"
    patch_start = 24831095
    patch_end = 24835982
    log("polishing patch {}:{}-{}".format(patch_chr, patch_start, patch_end))

    # Cut out the patch + flanking sequence
    log("extracting patch")
    #  TODO change to 100k flank sizes
    left_flank_size = 50000
    right_flank_size = 50000

    ff_ref = pysam.FastaFile(reference_file)
    patch_chr_len = ff_ref.get_reference_length(patch_chr)

    patch_start_flank = patch_start - left_flank_size
    if patch_start_flank < 0:
        patch_start_flank = 0
        left_flank_size = patch_start

    patch_end_flank = patch_end + right_flank_size
    if patch_end_flank > patch_chr_len:
        patch_end_flank = patch_chr_len
        right_flank_size = patch_chr_len - patch_end

    basename = make_coord_name(patch_chr, patch_start_flank, patch_end_flank)
    with open(output_path + basename + ".fasta", "w") as f:
        f.write(">" + make_coord_name(patch_chr, patch_start_flank, patch_end_flank) + "\n")
        f.write(ff_ref.fetch(patch_chr, patch_start_flank, patch_end_flank) + "\n")

    # Write the supporting reads to a fasta file
    log("extracting supporting reads")
    ff_reads = pysam.FastaFile(reads_file)
    with open(output_path + basename + ".reads.fasta", "w") as f:
        for i in supp_reads:
            f.write(">" + i + "\n")
            f.write(ff_reads.fetch(i) + "\n")

    # Run Medaka
    cmd = [
        "medaka_consensus",
        "-v",
        "-i",
        output_path + basename + ".reads.fasta",
        "-d",
        output_path + basename + ".fasta",
        "-o",
        basename + "_medaka"
    ]
    os.chdir(output_path)
    run_oe(cmd, output_path + basename + ".medaka.out", output_path + basename + ".medaka.err")
    os.chdir(cwd)

    # Filter the vcf
    cmd = [
        "tabix",
        output_path + basename + "_medaka/variants.vcf.gz",
        basename + ":" + str(left_flank_size) + "-" + str(left_flank_size + (patch_end-patch_start))
    ]
    run_oe(cmd, output_path + basename + "_medaka/variants.patch.vcf", output_path + basename + "_medaka/tabix.err")


if __name__ == "__main__":
    main()