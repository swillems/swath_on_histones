#!/usr/bin/python2.7


from collections import *


def main():
    print("Starting Analysis")
    sequences = getSequencesFromUniprotFile()
    writeToPeff(sequences)
    print("Analysis Finished")


def getSequencesFromUniprotFile():
    print("Reading Uniprot File")
    sequences = {}
    with open(UNIPROT_FILE, "rb") as infile:
        for line in infile:
            if line.startswith("ID"):
                seq_id = line.split()[1]
                sequences[seq_id] = {
                    "modifications": [],
                    "sequence": ""
                }
                while not line.startswith("//"):
                    line = infile.next()
                    if line.startswith("FT"):
                        if "MOD_RES" in line:
                            data = line.split()
                            ptm_description = line[34:]
                            ptm_description = ptm_description.split(".")[0]
                            ptm_name = ptm_description.split(";")[0]
                            sequences[seq_id]["modifications"].append(
                                (int(data[2]) - 1, ptm_name)
                            )
                        if "INIT_MET" in line:
                            sequences[seq_id]["modifications"].append((0, "mM-"))
                    elif line.startswith("SQ"):
                        line = infile.next()
                        sequence = ""
                        while line.startswith("  "):
                            data = line.split()
                            for d in data:
                                sequence += d
                            line = infile.next()
                        sequences[seq_id]["sequence"] = sequence
    return sequences


def writeToPeff(sequences):
    print("Writing results")
    pxfa_data = []
    for protein in sequences:
        pxfa_data.append("> " + protein + " 7 None" "\n")
        seq = sequences[protein]["sequence"]
        temp_mods = sequences[protein]["modifications"]
        mods = defaultdict(lambda: "")
        for x in temp_mods:
            try:
                mods[x[0]] += PTM_TRANSLATIONS[x[1]] + " "
            except KeyError:
                PTM_TRANSLATIONS[x[1]] = ""
                print("PTM not in parsing dict: " + x[1])
                mods[x[0]] += " "
        for i in range(len(seq)):
            pxfa_data.append(seq[i] + "\t")
            if i in mods:
                pxfa_data.append(mods[i])
            pxfa_data.append("\n")
    infile = "".join(pxfa_data).split("\n")
    with open(PEFF_FILE, "wb") as outfile:
        if INCLUDE_SNAPSHOT:
            exclude_list = writeSnapshot(outfile)
        else:
            exclude_list = []
        protein_sequence = ""
        protein_name = ""
        protein_max_ptm = 0
        ptm_list = []
        modifiable = False
        for line in infile:
            if not line:
                continue
            if line[0] == ">":
                if protein_name != "":
                    writeProtein(
                        outfile,
                        protein_name,
                        protein_sequence,
                        protein_max_ptm,
                        modifiable,
                        ptm_list,
                        exclude_list
                    )
                protein_sequence = ""
                data = line.split()
                protein_name = data[1]
                protein_max_ptm = data[2]
                ptm_list = []
                modifiable = False
            else:
                try:
                    data = line.split()
                    protein_sequence += data[0]
                    if len(data) > 1:
                        ptm_list.append(data[1:])
                        modifiable = True
                    else:
                        ptm_list.append([])
                except IndexError:
                    pass
        writeProtein(
            outfile,
            protein_name,
            protein_sequence,
            protein_max_ptm,
            modifiable,
            ptm_list,
            exclude_list
        )


def writeSnapshot(outfile):
    outfile.write(">Snapshot_Histone_H3.1 \ID=Snapshot_Histone_H3.1 \MaxPTMS=7 \ModRes=(0|mM-)(2|Me)(2|Me2)(2|Cit)(3|Ph)(4|Me)(4|Me2)(4|Me3)(4|Ac)(4|Cr)(4|Hib)(6|Ph)(8|Me)(8|Me2)(8|Cit)(9|Me)(9|Me2)(9|Me3)(9|Ac)(9|Bu)(9|Cr)(9|Hib)(10|Ac)(10|Ph)(11|Ph)(14|Me)(14|Me2)(14|Me3)(14|Ac)(14|Bu)(14|Hib)(14|Su)(14|Ub)(17|Me)(17|Me2)(17|Cit)(18|Me)(18|Me2)(18|Me3)(18|Ac)(18|Bu)(18|Cr)(18|Hib)(18|Su)(18|Fo)(18|Ub)(22|Ac)(23|Me)(23|Me2)(23|Me3)(23|Ac)(23|Pr)(23|Bu)(23|Cr)(23|Hib)(23|Su)(23|Fo)(23|Ub)(24|mAV)(26|Me)(26|Me2)(26|Cit)(27|Me)(27|Me2)(27|Me3)(27|Ac)(27|Bu)(27|Cr)(27|Hib)(27|Su)(27|Ub)(27|Ar)(28|Ac)(28|Ph)(31|mAS)(36|Me)(36|Me2)(36|Me3)(36|Ac)(36|Hib)(36|Ub)(37|Me)(37|Me2)(37|Me3)(37|Ac)(37|Ar)(40|Me)(40|Me2)(41|Ph)(42|Me)(42|Me2)(45|Ph)(56|Me)(56|Me2)(56|Me3)(56|Ac)(56|Pr)(56|Cr)(56|Hib)(56|Su)(56|Ma)(56|Fo)(56|Ub)(57|Ph)(63|Me)(63|Me2)(64|Me)(64|Me2)(64|Me3)(64|Ac)(64|Hib)(64|Fo)(71|mVM)(79|Me)(79|Me2)(79|Me3)(79|Ac)(79|Bu)(79|Cr)(79|Hib)(79|Su)(79|Fo)(79|Ub)(80|Ac)(80|Ph)(83|Me)(83|Me2)(86|Ph)(87|mSA)(89|mVI)(90|mMG)(96|mCS)(98|mAS)(107|Ph)(111|mAV)(115|Ac)(115|Bu)(122|Me)(122|Me2)(122|Me3)(122|Ac)(122|Bu)(122|Cr)(122|Hib)(122|Su)(122|Fo)(122|Ub)(128|Me)(128|Me2) \Length=136\n")
    outfile.write("MARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTE\n")
    outfile.write("LLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEACEAYLVGLFEDTNLCAIHAKRVTI\n")
    outfile.write("MPKDIQLARRIRGERA\n")
    outfile.write(">Snapshot_Histone_H4 \ID=Snapshot_Histone_H4 \MaxPTMS=7 \ModRes=(0|mM-)(1|Ph)(3|Me)(3|Me2)(3|Cit)(5|Me)(5|Me2)(5|Me3)(5|Ac)(5|Pr)(5|Bu)(5|Cr)(5|Hib)(8|Ac)(8|Pr)(8|Bu)(8|Cr)(8|Hib)(12|Me)(12|Me2)(12|Me3)(12|Ac)(12|Pr)(12|Bu)(12|Cr)(12|Hib)(12|Su)(12|Fo)(16|Me)(16|Me2)(16|Me3)(16|Ac)(16|Pr)(16|Bu)(16|Cr)(16|Hib)(16|Ar)(17|Me)(17|Me2)(17|Cit)(18|Ph)(19|Me)(19|Me2)(19|Cit)(20|Me)(20|Me2)(20|Me3)(20|Ac)(23|Me)(23|Me2)(31|Me)(31|Me2)(31|Me3)(31|Ac)(31|Pr)(31|Bu)(31|Hib)(31|Su)(31|Fo)(31|Ub)(35|Me)(35|Me2)(44|Pr)(44|Bu)(44|Hib)(47|Ph)(51|Ph)(51|OH)(55|Me)(55|Me2)(59|Me)(59|Me2)(59|Me3)(59|Cr)(59|Hib)(59|Fo)(59|Ub)(67|Me)(67|Me2)(72|Ph)(75|Ph)(77|Me)(77|Me2)(77|Me3)(77|Ac)(77|Pr)(77|Bu)(77|Cr)(77|Hib)(77|Su)(77|Fo)(77|Ub)(78|Me)(78|Me2)(79|Me)(79|Me2)(79|Me3)(79|Ac)(79|Pr)(79|Bu)(79|Hib)(79|Su)(79|Fo)(88|Ph)(88|OH)(91|Ac)(91|Pr)(91|Bu)(91|Cr)(91|Hib)(91|Su)(91|Fo)(91|Ub)(92|Me)(92|Me2) \Length=103\n")
    outfile.write("MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLK\n")
    outfile.write("VFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG\n")
    outfile.write(">Snapshot_Histone_H2A1H \ID=Snapshot_Histone_H2A1H \MaxPTMS=7 \ModRes=(0|mM-)(1|Ph)(3|Me)(3|Me2)(3|Cit)(5|Ac)(5|Bu)(5|Hib)(9|Me)(9|Me2)(9|Me3)(9|Ac)(9|Hib)(9|Su)(11|Me)(11|Me2)(13|Ac)(13|Ub)(13|Ar)(15|Ac)(15|Ub)(16|mTS)(29|Me)(29|Me2)(36|Ac)(36|Cr)(36|Hib)(36|Su)(36|Fo)(36|Ub)(39|OH)(40|mAS)(42|Me)(42|Me2)(50|Ph)(51|mLM)(59|Ph)(71|Me)(71|Me2)(74|Me)(74|Me2)(74|Me3)(74|Ac)(74|Hib)(75|Hib)(88|Me)(88|Me2)(95|Me)(95|Me2)(95|Me3)(95|Ac)(95|Pr)(95|Bu)(95|Cr)(95|Hib)(95|Su)(95|Fo)(95|Ub)(99|Me)(99|Me2)(99|Me3)(99|mKR)(101|Ph)(118|Me)(118|Me2)(118|Me3)(118|Ac)(118|Cr)(118|Hib)(118|Fo)(118|Ub)(119|Cr)(119|Ub)(120|Ph)(122|Ph)(125|Me)(125|Me2)(125|Me3)(125|Pr)(125|Cr)(125|Ub)(126|mAT)(127|Ac)(128|mG-)(129|Ac)(129|mK-) \Length=130\n")
    outfile.write("MSGRGKQGGKARAKAKTRSSRAGLQFPVGRVHRLLRKGNYAERVGAGAPVYLAAVLEYLT\n")
    outfile.write("AEILELAGNAARDNKKTRIIPRHLQLAIRNDEELNKLLGKVTIAQGGVLPNIQAVLLPKK\n")
    outfile.write("TESHHKAKGK\n")
    outfile.write(">Snapshot_Histone_H2B1C/E/F/G/I \ID=Snapshot_Histone_H2B1C/E/F/G/I \MaxPTMS=7 \ModRes=(0|mM-)(2|Ar)(2|mED)(4|mAS)(4|mAT)(5|Me)(5|Me2)(5|Me3)(5|Ac)(5|Bu)(5|Cr)(5|Hib)(5|Su)(5|Fo)(6|Ph)(11|Ac)(11|Bu)(11|Cr)(12|Me)(12|Me2)(12|Me3)(12|Ac)(12|Bu)(12|Cr)(12|Hib)(14|Ph)(15|Me)(15|Me2)(15|Me3)(15|Ac)(15|Bu)(15|Cr)(16|Ac)(16|Bu)(16|Cr)(18|mVI)(19|Ac)(20|Me)(20|Me2)(20|Me3)(20|Ac)(20|Bu)(20|Cr)(20|Hib)(20|Ub)(23|Me)(23|Me2)(23|Me3)(23|Ac)(23|Bu)(23|Cr)(23|Hib)(24|Ac)(24|Hib)(30|Ar)(32|Ph)(34|Me)(34|Me2)(34|Me3)(34|Cr)(34|Hib)(34|Su)(34|Fo)(34|Ub)(36|Ph)(37|OH)(39|mVI)(43|Me)(43|Me2)(43|Me3)(43|Ac)(43|Hib)(43|Su)(43|Fo)(46|Me)(46|Me2)(46|Me3)(46|Ac)(46|Bu)(46|Hib)(46|Su)(46|Fo)(46|Ub)(52|Ph)(56|Ph)(57|Me)(57|Me2)(57|Me3)(57|Ac)(57|Hib)(57|Ub)(64|Ph)(75|Ph)(79|Me)(79|Me2)(83|OH)(85|Me)(85|Me2)(85|Me3)(85|Ac)(85|Hib)(85|Su)(86|Me)(86|Me2)(87|Ph)(88|Ph)(91|Ph)(92|Me)(92|Me2)(99|Me)(99|Me2)(108|Me)(108|Me2)(108|Me3)(108|Ac)(108|Cr)(108|Hib)(108|Su)(108|Fo)(108|Ub)(115|Ph)(116|Me)(116|Me2)(116|Me3)(116|Ac)(116|Cr)(116|Hib)(116|Ma)(116|Su)(116|Fo)(116|Ub)(119|Ph)(120|Ac)(120|Hib)(120|Su)(120|Fo)(120|Ub)(124|mSA)(125|Ac) \Length=126\n")
    outfile.write("MPEPAKSAPAPKKGSKKAVTKAQKKDGKKRKRSRKESYSVYVYKVLKQVHPDTGISSKAM\n")
    outfile.write("GIMNSFVNDIFERIASEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVT\n")
    outfile.write("KYTSSK\n")
    outfile.write(">Snapshot_Histone_H1.2 \ID=Snapshot_Histone_H1.2 \MaxPTMS=7 \ModRes=(0|mM-)(1|Ph)(2|Ar)(3|Ph)(16|Me)(16|Me2)(16|Me3)(16|Ac)(16|Fo)(21|Me)(21|Me2)(21|Me3)(21|Ac)(22|Hib)(25|Me)(25|Me2)(25|Me3)(25|Ac)(25|Hib)(26|Me)(26|Me2)(26|Me3)(26|Hib)(33|Me)(33|Me2)(33|Me3)(33|Ac)(33|Cr)(33|Hib)(33|Su)(33|Fo)(33|Ub)(35|Ac)(35|Ph)(40|Ph)(45|Ac)(45|Hib)(45|Su)(45|Fo)(45|Ub)(50|Ac)(51|Me)(51|Me2)(51|Me3)(51|Ac)(51|Hib)(51|Ub)(53|Cit)(54|Ph)(62|Me)(62|Me2)(62|Me3)(62|Ac)(62|Hib)(62|Su)(62|Fo)(63|Me)(63|Me2)(63|Me3)(63|Ac)(63|Cr)(63|Hib)(63|Su)(63|Fo)(63|Ub)(70|Ph)(70|OH)(74|Ac)(74|Hib)(74|Fo)(74|Ub)(80|Hib)(80|Fo)(84|Ac)(84|Cr)(84|Hib)(84|Su)(84|Fo)(89|Me)(89|Me2)(89|Me3)(89|Ac)(89|Cr)(89|Hib)(89|Su)(89|Fo)(89|Ub)(96|Me)(96|Me2)(96|Me3)(96|Ac)(96|Cr)(96|Hib)(96|Su)(96|Fo)(96|Ub)(103|Ph)(105|Me)(105|Me2)(105|Me3)(105|Ac)(105|Su)(105|Ub)(109|Me)(109|Me2)(109|Me3)(109|Hib)(109|Ub)(112|Ac)(116|Hib)(116|Ub)(118|Me)(118|Me2)(118|Me3)(120|Hib)(120|Su)(126|Ub)(128|Me)(128|Me2)(128|Me3)(128|Hib)(135|Hib)(139|Fo)(139|Ub)(145|Ph)(147|Me)(147|Me2)(147|Me3)(147|Hib)(158|Cr)(158|Hib)(159|Fo)(159|Ub)(164|Ph)(167|Me)(167|Me2)(167|Me3)(167|Ac)(167|Cr)(167|Hib)(167|Ub)(168|Me)(168|Me2)(168|Me3)(168|Ac)(172|Ph)(179|Ph)(179|mTA)(186|Me)(186|Me2)(186|Me3)(190|Ac)(212|Hib) \Length=213\n")
    outfile.write("MSETAPAAPAAAPPAEKAPVKKKAAKKAGGTPRKASGPPVSELITKAVAASKERSGVSLA\n")
    outfile.write("ALKKALAAAGYDVEKNNSRIKLGLKSLVSKGTLVQTKGTGASGSFKLNKKAASGEAKPKV\n")
    outfile.write("KKAGGTKPKKPVGAAKKPKKAAGGATPKKSAKKTPKKAKKPAAATVTKKVAKSPKKAKVT\n")
    outfile.write("KPKKAAKSAAKAVKPKAAKPKVVKPKKAAPKKK\n")
    return [
        "H12_HUMAN",
        "H2A1B_HUMAN",
        "H2A1C_HUMAN",
        "H2A1D_HUMAN",
        "H2A1H_HUMAN",
        "H2A1J_HUMAN",
        "H2A2A_HUMAN",
        "H2B1B_HUMAN",
        "H2B1C_HUMAN",
        "H2B1D_HUMAN",
        "H2B1H_HUMAN",
        "H2B1J_HUMAN",
        "H2B1K_HUMAN",
        "H2B1N_HUMAN",
        "H2B1O_HUMAN",
        "H2B2E_HUMAN",
        "H31_HUMAN",
        "H31T_HUMAN",
        "H33_HUMAN",
        "H4_HUMAN",
    ]


def writeProtein(
    outfile,
    protein_name,
    protein_sequence,
    protein_max_ptm,
    modifiable,
    ptm_list,
    exclude_list
):
    if protein_name in exclude_list:
        return
    outfile.write(
        ">{} \ID={} \MaxPTMS={}".format(
            protein_name,
            protein_name,
            protein_max_ptm
        )
    )
    if modifiable:
        outfile.write(" \ModRes=")
        for i, ptms in enumerate(ptm_list):
            for ptm in ptms:
                outfile.write("({}|{})".format(i, ptm))
    outfile.write(" \Length=" + str(len(protein_sequence)) + "\n")
    while protein_sequence:
        outfile.write(protein_sequence[:60] + "\n")
        protein_sequence = protein_sequence[60:]


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Create a PEFF file from Uniprot.'
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

