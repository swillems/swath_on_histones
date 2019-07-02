#!/usr/bin/python2.7


import csv
from collections import *


in_file_name = "/home/sander/Documents/UGent/Proteomics/Laura/Relative_abundance/hESC_206_validaded_opgekuist.csv"
out_file_name = "/home/sander/Documents/UGent/Proteomics/Laura/Relative_abundance/hESC_206_validaded_opgekuist_ra.csv"

sequence_col_name = "Sequence"
mod_col_name = "Variable modifications ([position] description)"
normalized_abundance_col_name_pre = "171006_LDC_hESC_Naive"


sequences_ra = defaultdict(lambda: defaultdict(Counter))
sequences_total = defaultdict(Counter)
with open(in_file_name, "rb") as raw_infile:
    infile = csv.reader(raw_infile)
    header1 = infile.next()
    header2 = infile.next()
    header3 = infile.next()
    seq_cols = header3.index(sequence_col_name)
    mod_cols = header3.index(mod_col_name)
    int_cols = {
        i: sample for i, sample in enumerate(header3) if sample.startswith(
            normalized_abundance_col_name_pre
        )
    }
    for line in infile:
        sequence = line[seq_cols]
        mods = line[mod_cols].split("|")
        for int_col, sample in int_cols.iteritems():
            try:
                intensity = float(line[int_col])
            except:
                continue
            if intensity <= 0:
                continue
            sequences_total[sequence][sample] += intensity
            for mod in mods:
                sequences_ra[sequence][mod][sample] += intensity

with open(out_file_name, "wb") as raw_outfile:
    outfile = csv.writer(raw_outfile)
    outfile.writerow(
        [
            "Sequence",
            "PTM"
        ] + sorted(int_cols.values())
    )
    for sequence, mods in sorted(sequences_ra.iteritems()):
        for mod, samples in sorted(mods.iteritems()):
            if "Propionyl" in mod:
                continue
            col = [sequence, mod]
            for sample, count in sorted(samples.iteritems()):
                total = sequences_total[sequence][sample]
                col.append(count / total)
            outfile.writerow(col)

