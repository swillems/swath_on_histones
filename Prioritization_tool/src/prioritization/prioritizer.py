#!src/venv/bin/python


import re
import bisect
import numpy as np
import csv
import json
import itertools
import os
from zipfile import ZipFile, ZIP_DEFLATED
from shutil import rmtree
from time import asctime
from collections import defaultdict, Counter


class Prioritization(object):

    def __init__(self):
        self.PARAMETERS = {}
        self.AA_DICT = {}
        self.PTM_DICT = {}
        self.PROTEIN_DICT = {}
        self.ISOFORM_LIST = []
        self.PRECURSOR_LIST = []
        self.ANNOTATION_BOUNDARY_LIST = []
        self.PTMC_WEIGHT_COUNTER = Counter()
        self.PRIORITY_LIST = []

    def loadParameters(self, parameter_file_name):
        self.PARAMETERS = self.__load(
            file_name=parameter_file_name,
            text="parameters",
            count_top_level=False
        )

    def loadAADict(self):
        self.AA_DICT = self.__load(
            file_name=self.PARAMETERS["SAMPLE_PREP"]["AA_FILE_NAME"],
            text="amino acids"
        )

    def loadPTMDict(self):
        self.PTM_DICT = self.__load(
            file_name=self.PARAMETERS["SAMPLE_PREP"]["PTM_FILE_NAME"],
            text="PTMs"
        )

    def loadProteinDict(self):
        self.PROTEIN_DICT = self.__load(
            file_name=self.PARAMETERS["SAMPLE_PREP"]["PROTEIN_FILE_NAME"],
            text="proteins"
        )

    def loadPrecursorList(self):
        self.PRECURSOR_LIST = self.__load(
            file_name=self.PARAMETERS["MS"]["PRECURSOR_FILE_NAME"],
            text="precursors"
        )

    def __load(self, file_name, text=None, count_top_level=True):
        if text is not None:
            print(
                "{} > Loading {}".format(
                    asctime(),
                    text
                )
            )
        with open(file_name, 'rb') as infile:
            result = json.load(infile)
        if text is not None:
            if count_top_level:
                count = len(result)
            else:
                count = sum(len(i) for i in result.itervalues())
            print(
                "{} >    Found {} {}".format(
                    asctime(),
                    count,
                    text
                )
            )
        return result

    def runAll(self, save_as=None):
        print("{} > Start analysis".format(asctime()))
        self.loadAADict()
        self.loadPTMDict()
        self.loadProteinDict()
        self.loadPrecursorList()
        self.generateIsoforms()
        self.annotatePrecursors()
        self.calculatePTMCWeights()
        self.prioritizeAllPTMs()
        self.writePrioritizedPTMsToCSV()
        if save_as is not None:
            self.save(save_as)
        print("{} > Completed analysis".format(asctime()))

    def generateIsoforms(self):
        print("{} > Generating peptide isoforms".format(asctime()))
        n_ptms = self.getBasePTMs("n")
        c_ptms = self.getBasePTMs("c")
        non_selectable_ptms = set(
            ["null"] + [
                ptm for ptm in self.PARAMETERS["SAMPLE_PREP"]["PTMS"]["fixed"].values()
            ]
        )
        isoform_list = []
        for protein_name, (sequence, ptms) in self.PROTEIN_DICT.iteritems():
            full_ptms = self.getFullPTMList(sequence, ptms)
            digestion_points = self.digestionPoints(protein_name, sequence)
            for miscleavage in xrange(
                self.PARAMETERS["SAMPLE_PREP"]["MAX_MISCLEAVAGES"] + 1
            ):
                for peptide_start, peptide_end in zip(
                    digestion_points[:-(miscleavage + 1)],
                    digestion_points[miscleavage + 1:]
                ):
                    if peptide_end - peptide_start > self.PARAMETERS["SAMPLE_PREP"]["MAX_PEPTIDE_LENGTH"]:
                        continue
                    pep_seq = "n{}c".format(sequence[peptide_start: peptide_end])
                    try:
                        pep_mass = sum(self.AA_DICT[aa] for aa in pep_seq)
                    except KeyError:
                        continue
                    pep_ptms = [n_ptms] + full_ptms[
                        peptide_start: peptide_end
                    ] + [c_ptms]
                    pep_ptmcs = self.createIsoforms(pep_ptms)
                    for ptmc, count in pep_ptmcs.iteritems():
                        try:
                            ptmc_mass = sum(self.PTM_DICT[ptm]["MM"] for ptm in ptmc)
                        except TypeError:
                            continue
                        except KeyError:
                            continue
                        total_mass = pep_mass + ptmc_mass
                        ptmc_tuple = tuple(
                            set(
                                ptm for ptm in ptmc if ptm not in non_selectable_ptms
                            )
                        )
                        ptmc_representation = tuple(
                            "[{}] {}".format(
                                self.PTM_DICT[ptm]["TG"],
                                self.PTM_DICT[ptm]["CF"]
                            ) for ptm in ptmc_tuple
                        )
                        isoform = (
                            total_mass,
                            protein_name,
                            pep_seq,
                            ptmc,
                            ptmc_tuple,
                            ptmc_representation,
                            count
                        )
                        isoform_list.append(isoform)
        print(
            "{} >    Generated {} peptide isoforms".format(
                asctime(), len(isoform_list)
            )
        )
        self.ISOFORM_LIST = sorted(isoform_list)

    def getFullPTMList(self, sequence, ptms):
        ptm_list = [self.getBasePTMs(aa) for aa in sequence]
        for loc, ptm in ptms:
            ptm_list[loc].append(ptm)
        return ptm_list

    def getBasePTMs(self, aa):
        ptms = []
        if aa in self.PARAMETERS["SAMPLE_PREP"]["PTMS"]["fixed"]:
            ptms.append(self.PARAMETERS["SAMPLE_PREP"]["PTMS"]["fixed"][aa])
        else:
            ptms.append("null")
        if aa in self.PARAMETERS["SAMPLE_PREP"]["PTMS"]["variable"]:
            ptms += self.PARAMETERS["SAMPLE_PREP"]["PTMS"]["variable"][aa]
        return ptms

    def digestionPoints(self, protein_name, sequence):
        digestion_points = [0] + [
            i.start() + 1 for i in re.finditer(
                self.PARAMETERS["SAMPLE_PREP"]["CLEAVAGE_AAS"],
                sequence
            )
        ]
        if sequence[-1] not in self.PARAMETERS["SAMPLE_PREP"]["CLEAVAGE_AAS"]:
            digestion_points.append(len(sequence))
        return digestion_points

    def createIsoforms(self, pep_ptms):
        pep_ptms = sorted(pep_ptms, key=lambda s: len(s))
        previous_ptm_tuples = Counter()
        for ptm in pep_ptms[0]:
            previous_ptm_tuples[(ptm,)] = 1
        for ptms in pep_ptms[1:]:
            new_ptm_tuples = Counter()
            for previous_ptm_tuple, count in previous_ptm_tuples.items():
                for ptm in ptms:
                    ptm_list = list(previous_ptm_tuple)
                    bisect.insort(ptm_list, ptm)
                    new_ptm_tuples[tuple(ptm_list)] += count
            previous_ptm_tuples = new_ptm_tuples
        return new_ptm_tuples

    def annotatePrecursors(self):
        print("{} > Matching precursor candidates".format(asctime()))
        lower_limit = [0] * (self.PARAMETERS["MS"]["MAX_C13"] + 1)
        upper_limit = [0] * (self.PARAMETERS["MS"]["MAX_C13"] + 1)
        annotation_boundaries = []
        for mass in self.PRECURSOR_LIST:
            boundaries = []
            for isotope in xrange(self.PARAMETERS["MS"]["MAX_C13"] + 1):
                upper_mass = (
                    mass - self.PARAMETERS["CONSTANTS"]["DELTA_C_MASS"] * isotope
                ) * (
                    1 + self.PARAMETERS["MS"]["PPM"] / 1000000.
                )
                lower_mass = (
                    mass - self.PARAMETERS["CONSTANTS"]["DELTA_C_MASS"] * isotope
                ) * (
                    1 - self.PARAMETERS["MS"]["PPM"] / 1000000.
                )
                try:
                    while lower_mass > self.ISOFORM_LIST[lower_limit[isotope]][0]:
                        lower_limit[isotope] += 1
                except IndexError:
                    lower_limit[isotope] = len(self.ISOFORM_LIST)
                try:
                    while upper_mass > self.ISOFORM_LIST[upper_limit[isotope]][0]:
                        upper_limit[isotope] += 1
                except IndexError:
                    upper_limit[isotope] = len(self.ISOFORM_LIST)
                boundaries.append(
                    (lower_limit[isotope], upper_limit[isotope])
                )
            annotation_boundaries.append(boundaries)
        self.ANNOTATION_BOUNDARY_LIST = annotation_boundaries
        print(
            "{} >    Matched {} precursor candidates".format(
                asctime(),
                sum(1 for boundaries in annotation_boundaries if sum(
                    u - l for l, u in boundaries
                ) > 0)
            )
        )

    def calculatePTMCWeights(self):
        print("{} > Calculating PTMC weights".format(asctime()))
        ptmc_weights = Counter()
        precursor_weight = 1. / len(self.ANNOTATION_BOUNDARY_LIST)
        for boundaries in self.ANNOTATION_BOUNDARY_LIST:
            sequences_dict = defaultdict(Counter)
            for lower, upper in boundaries:
                for (
                    total_mass,
                    protein_name,
                    seq,
                    ptmc,
                    ptmc_tuple,
                    ptmc_representation,
                    count
                ) in self.ISOFORM_LIST[lower: upper]:
                    if len(ptmc_representation) <= self.PARAMETERS["PRIORITIZATION"]["MAX_PTMC_SIZE"]:
                        sequences_dict[seq][ptmc_representation] += 1
            try:
                seq_weight = precursor_weight / len(sequences_dict)
            except ZeroDivisionError:
                continue
            for seq, ptmc_counter in sequences_dict.iteritems():
                ptmc_weight = seq_weight / len(ptmc_counter)
                for ptmc in ptmc_counter.iterkeys():
                    ptmc_weights[ptmc] += ptmc_weight
        self.PTMC_WEIGHT_COUNTER = ptmc_weights

    def prioritizeAllPTMs(self):
        print("{} > Prioritizing PTMs".format(asctime()))
        sorted_ptms = sorted(
            set(
                ptm for ptmc in self.PTMC_WEIGHT_COUNTER.keys() for ptm in ptmc
            )
        )
        sorted_ptms_inverse = {ptm: i for i, ptm in enumerate(sorted_ptms)}
        sorted_ptmcs = []
        sorted_ptm_weights = np.zeros(len(self.PTMC_WEIGHT_COUNTER), dtype=float)
        M = np.zeros(
            (len(self.PTMC_WEIGHT_COUNTER), len(sorted_ptms))
        )
        for i, (ptmc, ptmc_weight) in enumerate(
            sorted(self.PTMC_WEIGHT_COUNTER.items())
        ):
            sorted_ptmcs.append(ptmc)
            sorted_ptm_weights[i] = ptmc_weight
            ptm_indices = [sorted_ptms_inverse[ptm] for ptm in ptmc]
            M[i, ptm_indices] = ptmc_weight
        print(
            "{} > \tFound {} PTMs in {} PTMCs".format(
                asctime(),
                len(sorted_ptms),
                len(self.PTMC_WEIGHT_COUNTER),
            )
        )
        unannotatables = (
            ["Unannotatable"],
            1 - np.sum(sorted_ptm_weights)
        )
        priorities = [unannotatables]
        try:
            unmodified = sorted_ptmcs.index(())
            unmodifieds = (
                ["No variable PTMs"],
                sorted_ptm_weights[unmodified]
            )
            priorities.append(unmodifieds)
        except ValueError:
            pass
        # np_sorted_ptms = np.array(sorted_ptms)
        # np_sorted_ptmcs = np.array(sorted_ptmcs)
        for i in xrange(
            self.PARAMETERS["PRIORITIZATION"]["SEQUENTIAL_SEARCHES"]
        ):
            next_ptms, next_ptmcs = self.prioritizeSinglePTMC(
                M,
                len(sorted_ptms),
                sorted_ptm_weights
            )
            if len(next_ptms) == 0:
                break
            M[next_ptmcs] = 0
            priorities.append(
                (
                    [sorted_ptms[ptm] for ptm in next_ptms],
                    np.sum(sorted_ptm_weights[next_ptmcs])
                )
            )
        self.PRIORITY_LIST = priorities

    # def prioritizeSinglePTMC(self, M):
    #     N = M.copy()
    #     N_sum = np.sum(N, axis=0)
    #     N_todo = (N_sum > 0).nonzero()[0]
    #     while len(N_todo) > self.PARAMETERS["PRIORITIZATION"]["MAX_PTMC_SIZE"]:
    #         worst_ptm = N_todo[np.argmin(N_sum[N_todo])]
    #         ptmcs_to_update = N[:, worst_ptm] > 0
    #         N[ptmcs_to_update] = 0
    #         N_sum = np.sum(N, axis=0)
    #         N_todo = (N_sum > 0).nonzero()[0]
    #     selected_ptmcs = np.any(N, axis=1).nonzero()[0]
    #     return N_todo, selected_ptmcs

    def prioritizeSinglePTMC(self, M, ptm_count, sorted_ptm_weights):
        max_results = [0, None, None]
        for combo in itertools.combinations(
            xrange(ptm_count),
            ptm_count - self.PARAMETERS["PRIORITIZATION"]["MAX_PTMC_SIZE"]
        ):
            full_ptmcs = np.logical_not(
                np.any(
                    M[:, combo],
                    axis=1
                )
            ).nonzero()[0]
            total_sum = np.sum(sorted_ptm_weights[full_ptmcs])
            if total_sum > max_results[0]:
                selected_ptms = np.ones(ptm_count, dtype=bool)
                selected_ptms[list(combo)] = False
                selected_ptms = selected_ptms.nonzero()[0]
                max_results[0] = total_sum
                max_results[1] = selected_ptms
                max_results[2] = full_ptmcs
        return max_results[1:]

    def writePrioritizedPTMsToCSV(self):
        print("{} > Saving priorities".format(asctime()))
        with open(
            self.PARAMETERS["PRIORITIZATION"]["PRIORITIES_FILE_NAME"], "wb"
        ) as raw_outfile:
            outfile = csv.writer(raw_outfile)
            weight = 0
            for ptmc, ptmc_weight in self.PRIORITY_LIST:
                weight += ptmc_weight
                outfile.writerow(
                    [
                        "; ".join(ptmc),
                        ptmc_weight,
                        weight
                    ]
                )

    def save(
        self,
        save_as_file_name,
        item_to_save=None,
    ):
        if item_to_save is not None:
            print("{} > Saving".format(asctime()))
            to_save_list = [(save_as_file_name, item_to_save)]
        else:
            print("{} > Saving analysis".format(asctime()))
            tmp_dir = ".tmp"
            if os.path.exists(tmp_dir):
                print("WARNING: .tmp dir exists!")
                print("Please check it for valuable information")
                raw_input("Press ctrl-C if you do not want to continue...")
                rmtree(tmp_dir)
            os.makedirs(tmp_dir)
            to_save_list = [
                (os.path.join(tmp_dir, "parameters.json"), self.PARAMETERS),
                (os.path.join(tmp_dir, "amino_acids.json"), self.AA_DICT),
                (os.path.join(tmp_dir, "ptms.json"), self.PTM_DICT),
                (os.path.join(tmp_dir, "proteins.json"), self.PROTEIN_DICT),
                (os.path.join(tmp_dir, "precursors.json"), self.PRECURSOR_LIST),
                (os.path.join(tmp_dir, "priorities.json"), self.PRIORITY_LIST),
            ]
        for file_name, item in to_save_list:
            with open(file_name, 'wb') as outfile:
                json.dump(item, outfile, indent=4, sort_keys=True)
        if item_to_save is None:
            with ZipFile(save_as_file_name, 'w', ZIP_DEFLATED) as zip_file:
                for root, dirs, files in os.walk(tmp_dir):
                    for file_name in files:
                        rooted_file_name = os.path.join(root, file_name)
                        zip_file.write(rooted_file_name, file_name)
            rmtree(tmp_dir)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Prioritize PTMs.'
    )
    parser.add_argument(
        "-p",
        "--parameter_file",
        help="The parameter file (required, string)",
        required=True,
    )
    parameter_file_name = parser.parse_args().parameter_file
    prioritization = Prioritization()
    prioritization.loadParameters(parameter_file_name)
    # prioritization.runAll(save_as="analysis.zip")
    prioritization.runAll()

# from src.prioritizer import ptm_prioritization
# parameter_file_name = "lib/default/parameters.json"
# prioritization = ptm_prioritization.Prioritization()
# prioritization.loadParameters(parameter_file_name)
# prioritization.runAll()
