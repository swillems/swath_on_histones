#!src/venv/bin/python


import json
import re
import urllib
import urllib2
from time import asctime


PROTON_MASS = 1.007276
HYDROGEN_MASS = 1.007825
AMINO_ACID_ABBREVIATIONS = {
    "Alanine": "A",
    "Arginine": "R",
    "Asparagine": "N",
    "Aspartate": "D",
    "Cysteine": "C",
    "Glutamate": "E",
    "Glutamine": "Q",
    "Glycine": "G",
    "Histidine": "H",
    "Isoleucine": "I",
    "Leucine": "L",
    "Lysine": "K",
    "Methionine": "M",
    "Phenylalanine": "F",
    "Proline": "P",
    "Serine": "S",
    "Threonine": "T",
    "Seleno-cysteine": "U",
    "Tryptophan": "W",
    "Tyrosine": "Y",
    "Valine": "V",
    "Any": "X",
    "None": "-",
    "N-terminal": "n",
    "C-terminal": "c",
}


def loadProteinDict(protein_in_file_name):
    print("{} > Loading proteins".format(asctime()))
    proteins = {}
    if protein_in_file_name.endswith(".fasta"):
        line_generator = loadFastaFile(protein_in_file_name)
    else:
        line_generator = loadUniprotFile(protein_in_file_name)
    for line in line_generator:
        if line.startswith("ID"):
            protein_name = line.split()[1]
            proteins[protein_name] = [None, []]
            cleave_init_met = False
            while not line.startswith("//"):
                if line.startswith("FT   "):
                    if "MOD_RES" in line:
                        location = int(line.split()[2])
                        ptm_description = line[34:]
                        ptm_description = ptm_description
                        line = line_generator.next()
                        while line.startswith("FT          "):
                            ptm_description += line[34:]
                            line = line_generator.next()
                        # ptm_name = re.split("[;.{]", ptm_description)[0]
                        # proteins[protein_name][1].append(
                        #     (location - 1, " ".join(ptm_name.split()))
                        # )
                        ptm_name = " ".join(ptm_description.split())
                        ptm_name = re.split(";|\. ", ptm_name)[0]
                        proteins[protein_name][1].append(
                            (location - 1, ptm_name)
                        )
                    elif "INIT_MET" in line:
                        # TODO, leave as variable?
                        cleave_init_met = True
                        line = line_generator.next()
                    else:
                        line = line_generator.next()
                else:
                    line = line_generator.next()
                # raw_input(line)
                if line.startswith("SQ   "):
                    line = line_generator.next()
                    sequence = ""
                    while line.startswith("  "):
                        data = line.split()
                        for d in data:
                            sequence += d
                        line = line_generator.next()
                    if cleave_init_met:
                        proteins[protein_name][0] = sequence[1:]
                        proteins[protein_name][1] = [
                            (
                                loc - 1, ptm
                            ) for loc, ptm in proteins[protein_name][1]
                        ]
                    else:
                        proteins[protein_name][0] = sequence
    print("{} > \tLoaded {} proteins".format(asctime(), len(proteins)))
    return proteins


def loadUniprotFile(protein_in_file_name):
    with open(protein_in_file_name, "rb") as infile:
        for line in infile:
            yield line


def loadFastaFile(fasta_file_name):
    print("{} > \tDownloading Uniprot accessions from fasta file".format(asctime()))
    proteins = []
    with open(fasta_file_name, "rb") as infile:
        for line in infile:
            if line.startswith(">"):
                raw_protein = line.split()[0]
                protein = raw_protein.split("|")[-1]
                proteins.append(protein)
    url = 'http://www.uniprot.org/uploadlists/'
    params = {
        'from': 'ACC',
        'to': 'ACC',
        'format': 'txt',
        'query': " ".join(proteins)
    }
    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    response = urllib2.urlopen(request)
    page = response.read()
    for line in page.splitlines():
        yield line


def loadPrecursorList(precursor_in_file_name):
    print("{} > Loading precursor list".format(asctime()))
    pep_mass = 0
    charge = 0
    precursor_masses = []
    with open(precursor_in_file_name, "rb") as mgf_file:
        for line in mgf_file:
            if line.startswith("PEPMASS"):
                pep_mass = float((line.split("=")[-1]).split()[-1])
            elif line.startswith("CHARGE"):
                charge = int((line.split("=")[-1]).split()[-1][0])
            if pep_mass and charge:
                mass = (pep_mass - PROTON_MASS) * charge
                precursor_masses.append(mass)
                pep_mass = 0
                charge = 0
    print(
        "{} > \tLoaded {} precursors".format(
            asctime(), len(precursor_masses)
        )
    )
    return sorted(precursor_masses)


def loadPTMDict(
    ptm_in_file_name,
):
    print("{} > Loading Uniprot PTM dict".format(asctime()))
    ptm_dict = {
        "null": {
            "MM": 0.,
            "TG": None,
            "CF": None,
            "KW": None,
            "DR": None,
            "UD": False,
        }
    }
    with open(ptm_in_file_name, "rb") as infile:
        for line in infile:
            if line.startswith("ID   "):
                name = line[5:].rstrip()
                on_residue_side = True
                terminal = ""
                ptm = {
                    "MM": None,
                    "TG": None,
                    "CF": None,
                    "KW": None,
                    "DR": None,
                    "UD": False,
                }
                for line in infile:
                    try:
                        key, value = line.split("   ")
                    except ValueError:
                        ptm = validatePTM(ptm, name, on_residue_side, terminal)
                        ptm_dict[name] = ptm
                        break
                    value = value.rstrip()
                    if key in ptm:
                        if key == "DR":
                            if "RESID" in value:
                                tmp, value = value.split("; ")
                            else:
                                continue
                        if key == "MM":
                            value = float(value)
                        elif key in ["TG", "KW", "DR"]:
                            value = value[:-1]
                        ptm[key] = value
                    elif key == "PA":
                        if value == "Amino acid backbone.":
                            on_residue_side = True
                    elif key == "PP":
                        if value == "C-terminal.":
                            terminal = "C-terminal"
                        elif value == "N-terminal.":
                            terminal = "N-terminal"
                    elif key == "FT" and value != "MOD_RES":
                        break
    print("{} > \tLoaded {} Uniprot PTMs".format(asctime(), len(ptm_dict)))
    return ptm_dict


def validatePTM(ptm, name, on_residue_side, terminal):
    if ptm["MM"] is None:
        print(
            "{} > \t{} has no mass definition".format(
                asctime(),
                name,
            )
        )
    if on_residue_side and terminal != "":
        print(
            "{} > \t{} target set as {}".format(
                asctime(),
                name,
                terminal,
            )
        )
        ptm["TG"] = terminal
    if ptm["CF"] is None:
        print(
            "{} > \t{} has no chemical function".format(
                asctime(),
                name,
            )
        )
    elif hasFalseParity(ptm):
        # try:
        if True:
            ptm["MM"] -= HYDROGEN_MASS
            composition = ptm["CF"].split()
            for i, element in enumerate(composition):
                if element.startswith("H"):
                    atom, count = element[0], int(element[1:])
                    count -= 1
                    new_element = "".join([atom, str(count)])
                    composition[i] = new_element
                    break
            ptm["CF"] = " ".join(composition)
            print(
                "{} > \t{} was mass corrected (negative parity of [{}])".format(
                    asctime(),
                    name,
                    ptm["CF"],
                )
            )
        # except TypeError:
        #     pass
    try:
        ptm["TG"] = AMINO_ACID_ABBREVIATIONS[ptm["TG"]]
    except KeyError:
        print(
            "{} > \t{} has unclear target {}".format(
                asctime(),
                name,
                ptm["TG"],
            )
        )
    return ptm


def hasFalseParity(ptm):
    pattern = re.compile("([a-zA-Z]+)([0-9,\-,\+]+)")
    parity_atoms = set(["H", "N", "P", "Cl", "Br", "F", "I"])
    atomcounts = ptm["CF"].split()
    even_count = True
    for atomcount in atomcounts:
        atom, count = pattern.match(atomcount).groups()
        if atom in parity_atoms:
            if int(count) % 2 == 1:
                even_count = not even_count
    if not even_count:
        return True
    else:
        return False


def save(out_file_name, out_dict):
    with open(out_file_name, 'wb') as outfile:
        json.dump(out_dict, outfile, indent=4, sort_keys=True)


if __name__ == "__main__":
    if False:
        ##############################
        FASTA_FILE_NAME = "files/M/Bovine_Histone_CoEx_Crap_IS_180321_1.fasta"
        if True:
            PROTEIN_IN_FILE_NAME = "files/test2/bovin_histones_uniprot.txt"
        else:
            PROTEIN_IN_FILE_NAME = FASTA_FILE_NAME
        PROTEIN_OUT_FILE_NAME = "files/test2/proteins.json"
        protein_dict = loadProteinDict(PROTEIN_IN_FILE_NAME)
        save(PROTEIN_OUT_FILE_NAME, protein_dict)
        ##############################
        PRECURSOR_IN_FILE_NAME = "files/test2/spectra.mgf"
        PRECURSOR_OUT_FILE_NAME = "files/test2/precursors.json"
        precursor_dict = loadPrecursorList(PRECURSOR_IN_FILE_NAME)
        save(PRECURSOR_OUT_FILE_NAME, precursor_dict)
        ##############################
        PTM_IN_FILE_NAME = "libs/uniprot/ptmlist.txt"
        PTM_OUT_FILE_NAME = "ptms.json"
        ptm_dict = loadPTMDict(PTM_IN_FILE_NAME)
        save(PTM_OUT_FILE_NAME, ptm_dict)
        ##############################
    if True:
        PTM_IN_FILE_NAME = "libs/uniprot/ptmlist.txt"
        PTM_OUT_FILE_NAME = "ptms.json"
        ptm_dict = loadPTMDict(PTM_IN_FILE_NAME)
        save(PTM_OUT_FILE_NAME, ptm_dict)
