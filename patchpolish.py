#!/usr/bin/env python

import os
import sys
import argparse

import pysam

from patchpolish_utilities.utilities import log, run_oe


def make_coord_name(h, s, e):
    return "{}_{}_{}".format(h, s, e)


def main():
    parser = argparse.ArgumentParser(description="Polish specified patch regions with supporting reads", usage="patchpolish.py <reference.fa> <reads.fa> <patches.bed>")
    parser.add_argument("reference", metavar="<patched_chromosomes.fa>", nargs='?', default="", type=str, help="reference fasta file containing raw patches. must not be gzipped.")
    parser.add_argument("reads", metavar="<ont_reads.fa>", nargs='?', default="", type=str, help="reads. must not be gzipped.")
    parser.add_argument("patches", metavar="<patches.bed>", nargs='?', default="", type=str, help="patches in bed format. fourth column is a CSV list of supporting reads.")
    parser.add_argument("-o", metavar="STR", type=str, default="patchpolish_output", help="output directory [patchpolish_output]")

    args = parser.parse_args()

    if not args.reference or not args.reads or not args.patches:
        parser.print_help()
        sys.exit()

    reference_file = os.path.abspath(args.reference)
    reads_file = os.path.abspath(args.reads)
    patch_file = os.path.abspath(args.patches)

    # Setup the output directory
    output_path = args.o.replace("/", "").replace(".", "")
    cwd = os.getcwd()
    output_path = cwd + "/" + output_path + "/"
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # Iterate over the patches
    with open(patch_file, "r") as fh:
        for line in fh:
            patch_chr, patch_start, patch_end, read_csv = line.rstrip().split("\t")
            patch_start, patch_end = int(patch_start), int(patch_end)
            supp_reads = read_csv.split(",")

            log("polishing patch {}:{}-{}".format(patch_chr, patch_start, patch_end))

            # Cut out the patch + flanking sequence
            log("extracting patch")
            left_flank_size = 100000
            right_flank_size = 100000

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
                "-h",
                output_path + basename + "_medaka/variants.vcf.gz",
                basename + ":" + str(left_flank_size + 1) + "-" + str(left_flank_size + (patch_end-patch_start))
            ]
            run_oe(cmd, output_path + basename + "_medaka/variants.patch.vcf", output_path + basename + "_medaka/tabix.err")

            cmd = [
                "liftover_vcf.py",
                output_path + basename + "_medaka/variants.patch.vcf",
                reference_file + ".fai",
                "-c",
                patch_chr,
                "-f",
                str(patch_start_flank)
            ]
            run_oe(cmd, output_path + basename + ".variants.patchpolish.vcf", output_path + "liftover.err")


if __name__ == "__main__":
    main()
