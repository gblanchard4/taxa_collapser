#!/usr/bin/env python
import argparse
import sys
from collections import defaultdict
from operator import add, div, mul

__author__ = "Gene Blanchard"
__email__ = "me@geneblanchard.com"

'''
Collapse OTUs Based on Taxa Assignments
p.s. this makes me feel dirty :(
'''


def check_biom_format(biomfile):
    with open(biomfile, 'r') as biom:
        line = biom.readline()
        if line.startswith('#'):
            return True
        else:
            return False


def clean_taxa_string(taxa):
    prefixes = ["k__", "; p__", "; c__", "; o__", "; f__", "; g__", "; s__"]
    while any(taxa.endswith(prefix) for prefix in prefixes):
        taxa = taxa.replace(prefixes[-1], '')
        prefixes.pop()
    return taxa


def main():
    # Argument Parser
    parser = argparse.ArgumentParser(description='Collapses taxa in a biom file')

    # Input file
    parser.add_argument('-i', '--input', dest='input', required=True, help='The tab seperated biom file')
    # Filter percentage
    parser.add_argument('-f', '--filter', dest='filter', type=int, default=.01, help='The minumum percent abundance to keep out DEFAULT=.01 or 1%')
    # Output file
    parser.add_argument('-o', '--output', dest='output', required=True, help='The output file')

    # Parse arguments
    args = parser.parse_args()
    infile = args.input
    min_abund = args.filter
    outfile = args.output

    # Confirm biom is tsv
    if not check_biom_format(infile):
        print "File not a TSV!\nTry to convert using:\nbiom convert -i otus.biom --to-tsv -o otus.txt --header-key taxonomy"
        sys.exit()

    # Read the biom to a dictionaty of taxa
    taxa_dict = defaultdict(list)
    with open(infile, 'r') as biom:
        for line in biom:
            if line.startswith('#OTU ID\t'):
                headerline = line
                header = line.rstrip('\n').split('\t')
                otu = header[0]
                taxa = header[-1]
                samples = header[1:-1]
                sums = [0] * len(samples)
            if not line.startswith('#'):
                data = line.rstrip('\n').split('\t')
                otu = data[0]
                taxa = data[-1]
                counts = map(float, data[1:-1])
                sums = map(add, sums, counts)
                # Check if it's an empty row
                if not sum(counts) == 0:
                    taxa = clean_taxa_string(taxa)
                    if taxa_dict[taxa]:
                        taxa_dict[taxa] = map(add, taxa_dict[taxa], counts)
                    else:
                        taxa_dict[taxa] = counts

    # # Convert to relative abundance
    # for key in taxa_dict:
    #     taxa_dict[key] = map(div, taxa_dict[key], sums)
    # Format Output
    with open(outfile, 'w') as outhandle:
        outhandle.write(headerline)
        for index, key in enumerate(taxa_dict.keys()):
            outhandle.write("{}\t{}\t{}\n".format(index, '\t'.join(map(str, taxa_dict[key])), key))


if __name__ == '__main__':
    main()
