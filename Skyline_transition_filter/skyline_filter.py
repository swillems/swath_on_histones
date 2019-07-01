#!/usr/bin/python


# To run, copy paste the following command to a shell:
# "python skyline_filter.py -p 'Test/parameter_file.py'"
# with correct location of a parameter file


import csv
from collections import *
import re
import sys
import os
import json
from operator import attrgetter


def main():
    print("Start Analysis")
    if not os.path.exists(OUTPUT_PATH):
        os.makedirs(OUTPUT_PATH)
    filtered_transitions = []
    with open(FILTERED_TRANSITION_FILE_NAME, "wb") as out_file:
        csv_out_file = csv.writer(out_file, delimiter=DELIMITER)
        for filtered_transition in generateFilteredTransitions():
            csv_out_file.writerow(filtered_transition)
            filtered_transitions.append(filtered_transition)
    with open(UNIQUE_TRANSITION_FILE_NAME, "wb") as out_file:
        csv_out_file = csv.writer(out_file, delimiter=DELIMITER)
        for unique_transition in generateUniqueTransitions(
            filtered_transitions
        ):
            csv_out_file.writerow(unique_transition)
    print("Finished Analysis")


def generateFilteredTransitions():
    with open(TRANSITION_FILE_NAME, "rb") as in_file:
        csv_in_file = csv.reader(in_file, delimiter=DELIMITER)
        header = csv_in_file.next()
        sequence_column = getColumnIndex(header, "pep_mod_seq")
        precursor_mz_column = getColumnIndex(header, "pep_mz")
        product_mz_column = getColumnIndex(header, "prod_mz")
        ion_type_column = getColumnIndex(header, "prod_ion_type")
        protein_name_column = getColumnIndex(header, "protein")
        if CURATE_PTMS:
            # curated_protein_dict = createCuratedProteinDict()
            curated_protein_dict = createCuratedProteinDictFromJSON()
        for row in csv_in_file:
            ion_type = row[ion_type_column]
            skyline_sequence_string = row[sequence_column]
            # Exact match with database file!
            protein_name = row[protein_name_column].split("|")[-1]
            precursor_mz = row[precursor_mz_column]
            product_mz = row[product_mz_column]
            if REMOVE_PRECURSOR_TRANSITIONS and ion_type == "precursor":
                continue
            if CURATE_PTMS and not sequenceIsCurated(
                curated_protein_dict,
                protein_name,
                skyline_sequence_string
            ):
                continue
            transition = (
                skyline_sequence_string,
                precursor_mz,
                product_mz,
            )
            yield transition


def getColumnIndex(header, column_string):
    return header.index(
        TRANSITION_FILE_COLUMN_NAMES[column_string]
    )


def createCuratedProteinDictFromJSON():
    PTM_ABBREVIATIONS = {'N-acetylserine': 'Ac', 'Asymmetricdimethylarginine': 'Me2', 'Phosphoserine': 'Ph', 'N6,N6,N6-trimethyllysine': 'Me3', 'N-acetylvaline': 'Ac', 'Citrulline': 'Cit', 'Omega-N-methylatedarginine': 'Me', 'N6-crotonyllysine': 'Cr', 'Phosphothreonine': 'Ph', 'N-acetylmethionine': 'Ac', 'O-AMP-threonine': '', 'N-pyruvate2-iminyl-valine': '', 'N-acetylglycine': 'Ac', 'Omega-N-methylarginine': 'Me', 'Pyrrolidonecarboxylicacid': 'Py', 'O-AMP-tyrosine': '', 'N6,N6-dimethyllysine': 'Me2', 'N-acetylthreonine': 'Ac', 'N6-malonyllysine': 'Ma', 'Cysteinemethylester': 'Me', '3-hydroxyproline': 'OH', 'Blockedaminoend(Cys)': '', 'Blockedaminoend(Thr)': '', 'N6-(pyridoxalphosphate)lysine': '', 'Methionine(R)-sulfoxide': 'Ox', 'N6-succinyllysine': 'Su', 'Methioninesulfoxide': 'Ox', 'Symmetricdimethylarginine': 'Me2', 'N-acetylalanine': 'Ac', 'Cysteinepersulfide': 'Pe', 'N5-methylglutamine': 'Me', 'mM-': 'mM-', 'N-acetylaspartate': 'Ac', 'Hypusine': '', 'ADP-ribosylasparagine': 'Ar', 'S-nitrosocysteine': 'Ni', 'N-acetylproline': 'Ac', 'Deamidatedasparagine': 'Am-', 'Nitratedtyrosine': '', 'ADP-ribosylcysteine': 'Ar', 'Phosphotyrosine': 'Ph', '5-glutamyl': '', 'Tele-methylhistidine': 'Me', '(Z)-2,3-didehydrotyrosine': 'Hy2', '5-hydroxylysine': '', 'N6-methylatedlysine': 'Me', 'Cysteinederivative': '', 'N-acetylglutamate': 'Ac', 'Diphthamide': '', '(3S)-3-hydroxyhistidine': 'OH', 'Blockedaminoend(Ser)': '', 'PolyADP-ribosylglutamicacid': '', 'Dimethylatedarginine': 'Me2', 'Allysine': '', 'Pros-methylhistidine': '', 'S-glutathionylcysteine': '', '4-hydroxyproline': 'OH', 'N6-acetyllysine': 'Ac', 'N6-methyllysine': 'Me'}
    with open(DATABASE_FILE_NAME, "rb") as infile:
        tmp_curated_protein_dict = json.load(infile)
    curated_protein_dict = {}
    for protein_name, (protein_sequence, ptms) in tmp_curated_protein_dict.items():
        ptm_string = "".join(
            "({}|{})".format(
                loc,
                PTM_ABBREVIATIONS[
                    "".join(ptm.split())
                ]
            ) if ptm in PTM_ABBREVIATIONS else "({}|??)".format(loc) for loc, ptm in ptms
        )
        curated_protein_dict[protein_name] = processProtein(
            protein_name,
            protein_sequence,
            ptm_string,
        )
    return curated_protein_dict


def createCuratedProteinDict():
    curated_protein_dict = {}
    with open(DATABASE_FILE_NAME, "rb") as in_file:
        for line in in_file:
            if line.startswith(">"):
                try:
                    curated_protein_dict[protein_name] = processProtein(
                        protein_name,
                        protein_sequence,
                        ptm_string,
                    )
                except UnboundLocalError:
                    pass
                protein_sequence = ""
                ptm_string = ""
                for protein_property in line.split():
                    if protein_property.startswith("\\ID="):
                        protein_name = protein_property.split("=")[1]
                    elif protein_property.startswith("\\ModRes="):
                        ptm_string = protein_property.split("=")[1]
            else:
                protein_sequence += line.rstrip()
        curated_protein_dict[protein_name] = processProtein(
            protein_name,
            protein_sequence,
            ptm_string,
        )
    return curated_protein_dict


def processProtein(
    protein_name,
    protein_sequence,
    ptm_string,
):
    parsed_ptms = parsePTMFromProteinSequence(
        ptm_string,
        protein_sequence
    )
    protein = Protein(
        protein_name,
        protein_sequence,
        parsed_ptms
    )
    if USE_EXPLICIT_SNAPSHOT_VARIANTS:
        protein.generateExplicitVariants()
    return protein


def parsePTMFromProteinSequence(ptm_string, sequence):
    ptm_set = set()
    pattern = "[^()]+"
    raw_ptms = re.findall(pattern, ptm_string)
    for raw_ptm in raw_ptms:
        loc, mod = raw_ptm.split("|")
        ptm = PTM(mod, sequence[int(loc)], loc)
        ptm_set.add(ptm)
    return ptm_set


def sequenceIsCurated(
    curated_protein_dict,
    protein_name,
    skyline_sequence_string
):
    protein = curated_protein_dict[protein_name]
    sequence, ptms = parsePTMFromSkylineSequenceString(
        skyline_sequence_string
    )
    if protein.hasCuratedSequence(sequence, ptms):
        return True
    if USE_EXPLICIT_SNAPSHOT_VARIANTS:
        for protein_name, protein in curated_protein_dict.iteritems():
            # Slow manner of checking snapshot variants
            if protein_name.startswith("Snapshot"):
                if protein.hasCuratedSequenceVariant(sequence, ptms):
                    return True
    return False


def parsePTMFromSkylineSequenceString(skyline_sequence_string):
    sequence = ""
    elements = re.split("[\[,\]]", skyline_sequence_string)
    parsed_mods = []
    for element in elements:
        try:
            mass = float(element)
        except ValueError:
            sequence += element
            continue
        aa = sequence[-1]
        # NOTE: Skyline n-term is considered as a ptm on
        #       the first amino acid!
        loc = len(sequence)
        mod = MASS_TRANSLATIONS[round(mass,4)]
        if loc == 1:
            if mod in CHEMICAL_PTMS["N-term"]:
                aa = "N-term"
        ptm = PTM(mod, aa, loc)
        if ptm.aa in CHEMICAL_PTMS:
            if ptm.mod in CHEMICAL_PTMS[ptm.aa]:
                ptm.biological = False
                if CHEMICAL_PTMS[ptm.aa][ptm.mod] == "Fixed":
                    ptm.variable = False
        if len(ptm.aa) != 1 or ptm.aa < "A" or ptm.aa > "Z":
            ptm.biological = False
        parsed_mods.append(ptm)
    return sequence, parsed_mods


def generateUniqueTransitions(filtered_transitions):
    windows = createSWATHWindows()
    windows = fillSWATHWindowsWithData(windows, filtered_transitions)
    for unique_transition in generateUniqueTransitionsPerWindow(windows):
        yield unique_transition


def createSWATHWindows():
    windows = {}
    with open(WINDOW_FILE_NAME, "rb") as in_file:
        csv_in_file = csv.reader(in_file, delimiter=DELIMITER)
        row = csv_in_file.next()
        lower_limit, upper_limit, window_size = map(float, row)
        windows[(0, lower_limit)] = defaultdict(list)
        windows[(lower_limit, upper_limit)] = defaultdict(list)
        for row in csv_in_file:
            lower_limit, upper_limit, window_size = map(float, row)
            windows[(lower_limit, upper_limit)] = defaultdict(list)
        windows[(upper_limit, sys.float_info.max)] = defaultdict(list)
    return windows


def fillSWATHWindowsWithData(windows, filtered_transitions):
    limits = sorted(windows)
    upper_ppm = (1 + PRECURSOR_PPM_ERROR / 1000000.)
    lower_ppm = (1 - PRECURSOR_PPM_ERROR / 1000000.)
    for transition in filtered_transitions:
        precursor_mass, product_mass = map(float, transition[1:])
        for lower_limit, upper_limit in limits:
            if upper_limit * upper_ppm < precursor_mass:
                continue
            if lower_limit * lower_ppm > precursor_mass:
                break
            window = (lower_limit, upper_limit)
            windows[window][product_mass].append(transition)
    return windows


def generateUniqueTransitionsPerWindow(windows):
    upper_ppm = (1 + PRODUCT_PPM_ERROR / 1000000.)
    for transition_dict in windows.itervalues():
        transition_list = sorted(transition_dict.iterkeys())
        unique = True
        for i, transition in enumerate(transition_list):
            try:
                if transition * upper_ppm > transition_list[i + 1]:
                    unique = False
                    continue
            except IndexError:
                pass
            if unique:
                rows = transition_dict[transition]
                if len(rows) == 1:
                    yield rows[0]
            unique = True


class PTM(object):

    def __init__(self, mod, aa, loc):
        self.mod = mod
        self.loc = int(loc)
        self.aa = aa
        self.biological = True
        self.variable = True
        # self.id = "{0}@{1}{2}".format(mod, aa, loc)
        self.id = "{0}+{1}({2})".format(self.mod, self.aa, self.loc)

    def __str__(self):
        return "{0}".format(self.id)

    def __eq__(self, other):
        if self.id == other.id:
            return True
        return False

    def __hash__(self):
        return hash(self.id)

    def isSNP(self):
        return self.mod.startswith("m") and self.mod[-1] != "-"

    def shiftLocation(self, offset):
        return PTM(self.mod, self.aa, self.loc + offset)

    def changeButyrylToMethyl(self):
        if self.mod == "Bu":
            return PTM("Me", self.aa, self.loc)
        return self


class Protein(object):

    def __init__(self, accesion, sequence, ptms):
        self.accesion = accesion
        self.sequence = sequence
        self.ptms = ptms
        self.id = accesion

    def __str__(self):
        return "{}".format(self.id)

    def generateExplicitVariants(self):
        self.sequence_variants = set([self.sequence])
        snps = [ptm for ptm in self.ptms if ptm.isSNP()]
        for snp in sorted(snps, key=attrgetter('loc')):
            loc = snp.loc
            aa = snp.aa
            seq_end = self.sequence[loc:]
            new_seqs = set()
            for seq in self.sequence_variants:
                seq_begin = seq[:loc - 1]
                new_seqs.add(
                    "".join(
                        [seq_begin, aa, seq_end]
                    )
                )
            self.sequence_variants |= new_seqs

    def hasCuratedSequence(self, sequence, ptms):
        if self.hasCuratedSubSequence(sequence, ptms):
            return True
        return False

    def hasCuratedSequenceVariant(self, sequence, ptms):
        for ref_sequence in self.sequence_variants:
            if sequence in ref_sequence:
                if self.hasCuratedSubSequence(sequence, ptms):
                    return True
        return False

    def hasCuratedSubSequence(self, sequence, ptms):
        indices = [m.start() for m in re.finditer(sequence, self.sequence)]
        for index_offset in indices:
            for ptm in ptms:
                if not ptm.biological:
                    continue
                reindexed_ptm = ptm.shiftLocation(index_offset - 1)
                if reindexed_ptm not in self.ptms:
                    if reindexed_ptm.changeButyrylToMethyl() in self.ptms:
                        continue
                    break
            else:
                return True


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Filter a Skyline file.'
    )
    parser.add_argument(
        "-p",
        "--parameter_file",
        help="The parameter file (required, string)",
        required=True,
    )
    parameter_file_name = parser.parse_args().parameter_file
    execfile(parameter_file_name)
    main()
