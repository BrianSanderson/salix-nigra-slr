#!/usr/bin/env python

"""
A script that takes a list of candidate sequences for targeted sequence capture
design, and a VCF created by mpileup as input, and quantifies the polymorphism
in the candidate regions

note: "import vcf" refers to the module pyvcf
"""

import csv
import argparse
import vcf
from tqdm import tqdm

__author__ = "Brian J. Sanderson <brian@biologicallyrelevant.com>"


def probe_locs(candidate_file, num_cand):
    """Defines the positions on each chromosome covered by candidate regions

    Parameters
    ----------
    candidate_file : str
        The path to a text file with the following fields
            0: the chromosome
            1: the start coordinate
            2: the stop coordinate
            3: the length of the region
            4: the name of the locus

    num_cand : int
        (Optional) Number of candidate regions. Will provide better progress
        output.

    Returns
    -------
    probeSet : dict
        Dictionary of sets for each chromsome (or contig or scaffold)
        that includes the range of positions for which the candidate probe
        sequences cover.
    """
    print("Processing candidate regions")
    with tqdm(total=num_cand) as pbar:
        with open(candidate_file) as f:
            probeSet = {}
            for line in f:
                chrom = line.split()[0]
                start = int(line.split()[1])
                stop = int(line.split()[2])
                if chrom in probeSet:
                    probeSet[chrom].update(set(range(start, stop + 1)))
                else:
                    probeSet[chrom] = set(range(start, stop + 1))
                pbar.update(1)
    print("")
    return probeSet


def quantify_poly(probeSet, vcf_file, num_vcf):
    """Quantify SNPs and indels for each position in the candidate probes

    Parameters
    ----------
    probeSet : dict
        The dictionary of sets from probe_locs()
    vcf_file : str
        The path to the VCF file containing the variant calls
    num_vcf : int
        (Optional) Number of candidate regions. Will provide better progress
        output.
    Returns
    -------
    snps : dict
        A dictionary of sets of SNP locations for each chromosome
    indels: dict
        A dictionary of sets of indel locations for each chromosome
    """
    vcf_reader = vcf.Reader(filename=vcf_file)
    indels = {}
    snps = {}
    print("Quantifying SNPs and indels in VCF")
    with tqdm(total=num_vcf) as pbar:
        for record in vcf_reader:
            try:
                record.POS in probeSet[record.CHROM]
            except ValueError:
                pass
            else:
                if record.POS in probeSet[record.CHROM]:
                    if record.is_monomorphic:
                        continue
                    elif record.is_snp:
                        if record.CHROM in snps:
                            snps[record.CHROM].update(set([record.POS]))
                        else:
                            snps[record.CHROM] = set([record.POS])
                    elif record.is_indel:
                        if record.CHROM in indels:
                            indels[record.CHROM].update(set([record.POS]))
                        else:
                            indels[record.CHROM] = set([record.POS])
                    else:
                        continue
                else:
                    continue
            finally:
                pbar.update(1)
    print("")
    return snps, indels


def summarize_poly(snps, indels, candidate_file, out_file, num_cand):
    """Write out a tab-delimited file that summarizes polymorphism per region

    Parameters
    ----------
    snps : dict
        Dictionary of sets for each chromosome with the positions of SNPs
    indels : dict
        Dictionary of sets for each chromosome with the positions of indels
    candidate_file : str
        The path to a text file with the following fields
            0: the chromosome
            1: the start coordinate
            2: the stop coordinate
            3: the length of the region
            4: the name of the region
    out_file : str
        The path to a tab-delimited text file to be written out
    num_cand : int
        (Optional) Number of candidate regions. Will provide better progress
        output.

    Notes
    -----
    This function writes out a tab-delimited text file with the following
    fields
        0: the chromosome
        1: the start coordinate
        2: the stop coordinate
        3: the length of the region
        4: the count of indels in region
        5: the count of SNPs in the region
        6: the name of the locus
    """
    with open(out_file, 'w', newline='\n') as ofile:
        writer = csv.writer(ofile, delimiter='\t',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['CHROM', 'start', 'stop', 'length', 'indelCount',
                         'snpCount', 'locusName'])
        with tqdm(total=num_cand) as pbar:
            print("Writing out summary for each candidate region")
            with open(candidate_file) as g:
                for line in g:
                    chrom = line.split()[0]
                    start = int(line.split()[1])
                    stop = int(line.split()[2])
                    indelCount = 0
                    snpCount = 0
                    try:
                        indels[chrom]
                    except KeyError:
                        pass
                    else:
                        for j in range(start, stop + 1):
                            if j in indels[chrom]:
                                indelCount += 1
                            else:
                                continue
                    try:
                        snps[chrom]
                    except KeyError:
                        pass
                    else:
                        for j in range(start, stop + 1):
                            if j in snps[chrom]:
                                snpCount += 1
                            else:
                                continue
                    writer.writerow([line.split()[0], line.split()[1],
                                     line.split()[2], line.split()[3],
                                     indelCount, snpCount,
                                     line.split()[4]])
                    pbar.update(1)
        print("")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    requiredNamed = parser.add_argument_group("required named arguments")
    requiredNamed.add_argument("-c", dest="candidate_file",
                               metavar="example.txt",
                               help="file containing region names and ranges",
                               required=True)
    requiredNamed.add_argument("-i", dest="vcf_file", metavar="input.vcf",
                               help="the input VCF file", required=True)
    requiredNamed.add_argument("-o", dest='out_file', metavar="output.txt",
                               help="the output polymorhpism file",
                               required=True)
    optionalNamed = parser.add_argument_group("optional named arguments")
    optionalNamed.add_argument("--num_cand", dest="num_cand",
                               metavar="num_cand",
                               help="number of candidate regions; if provided "
                               "will generate an informative progress bar",
                               required=False, default=None)
    optionalNamed.add_argument("--num_vcf", dest="num_vcf", metavar="num_vcf",
                               help="number of sites in VCF; if provided will "
                               "generate an informative progress bar",
                               required=False, default=None)
    args = parser.parse_args()

    if args.num_cand is None:
        num_cand = '-inf'
    else:
        num_cand = int(args.num_cand)
    if args.num_vcf is None:
        num_cand = '-inf'
    else:
        num_vcf = int(args.num_vcf)

    probeSet = probe_locs(args.candidate_file, num_cand)
    snps, indels = quantify_poly(probeSet, args.vcf_file, num_vcf)
    summarize_poly(snps, indels, args.candidate_file, args.out_file,
                   num_cand)
