#!venv/bin/python


import csv
import numpy as np


base_path = "/home/sander/Documents/LabFBT/Colleague_work/180306_LdC_SWATH/Relative_abundance/190129_HDAC_BCS_SWATHpaper/"
inputs = [
    "HDAC",
    "HDAC_norm",
    "BCC_453",
    "BCC_468",
]
for input in inputs:
    if input == "HDAC":
        in_file_name = base_path + "190130_HDAC_RA_SW_limma.csv"
        out_file_name = base_path + "190130_HDAC_RA_SW_limma_matrices.csv"
    if input == "HDAC_norm":
        in_file_name = base_path + "181212_HDAC_norm_parse_edit_RA_forStats_SW_limma.csv"
        out_file_name = base_path + "181212_HDAC_norm_parse_edit_RA_forStats_SW_limma_matrices.csv"
    if input == "BCC_453":
        in_file_name = base_path + "181211_BCC_Norm_parse_edit_RA_forStats_453_SW_limma.csv"
        out_file_name = base_path + "181211_BCC_Norm_parse_edit_RA_forStats_453_SW_limma_matrices.csv"
    if input == "BCC_468":
        in_file_name = base_path + "181211_BCC_Norm_parse_edit_RA_forStats_468_SW_limma.csv"
        out_file_name = base_path + "181211_BCC_Norm_parse_edit_RA_forStats_468_SW_limma_matrices.csv"
    peps = []
    contrasts = []
    with open(in_file_name, "r") as raw_infile:
        infile = csv.reader(raw_infile)
        header = next(infile)
        for row in infile:
            if row[0].startswith("Time"):
                t1, t2 = row[0].split("-")
                t1 = int(t1[4:])
                t2 = int(t2[4:])
                contrasts.append((t1, t2))
            if row[0].startswith("adj.P.Val"):
                peps.append([float(i) for i in row[1:]])
    peps = np.array(peps)
    max_size = np.max(contrasts)
    with open(out_file_name, "w") as raw_outfile:
        outfile = csv.writer(raw_outfile)
        for pep_index, pep_name in enumerate(header[1:]):
            mat = np.zeros((max_size - 1, max_size - 1))
            for contrast_index, (time_a, time_b) in enumerate(contrasts):
                mat[time_b - 2, time_a - 1] = peps[contrast_index, pep_index]
            tmp = outfile.writerow([pep_name])
            for row in mat:
                tmp = outfile.writerow(row)
            tmp = outfile.writerow("")
