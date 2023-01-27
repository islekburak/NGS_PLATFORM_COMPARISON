

CHROMOSOME_MAPPING = {
    "NC_012920.1": "M",
    "NC_000001.10": "1",
    "NC_000002.11": "2",
    "NC_000003.11": "3",
    "NC_000004.11": "4",
    "NC_000005.9": "5",
    "NC_000006.11": "6",
    "NC_000007.13": "7",
    "NC_000008.10": "8",
    "NC_000009.11": "9",
    "NC_000010.10": "10",
    "NC_000011.9": "11",
    "NC_000012.11": "12",
    "NC_000013.10": "13",
    "NC_000014.8": "14",
    "NC_000015.9": "15",
    "NC_000016.9": "16",
    "NC_000017.10": "17",
    "NC_000018.9": "18",
    "NC_000019.9": "19",
    "NC_000020.10": "20",
    "NC_000021.8": "21",
    "NC_000022.10": "22",
    "NC_000023.10": "X",
    "NC_000024.9": "Y"
}

EVIDENCE_CODE_MAPPING = {
    'PVS1': 'PVS1',
    'PVS1-': 'PVS1_Strong',
    'PVS1--': 'PVS1_Moderate',
    'PVS1---': 'PVS1_Supporting',
    'PP5+++': 'PP5_Very_Strong',
    'PS1': 'PS1',
    'PS2': 'PS2',
    'PS3': 'PS3',
    'PS4': 'PS4',
    'PP5++': 'PP5_Strong',
    'PM5+': 'PM5_Strong',
    'PM1': 'PM1',
    'PS1-': 'PS1_Moderate',
    'PM2': 'PM2',
    'PM3': 'PM3',
    'PM4': 'PM4',
    'PM5': 'PM5',
    'PM6': 'PM6',
    'PP5+': 'PP5_Moderate',
    'PP1': 'PP1',
    'PP2': 'PP2',
    'PP3': 'PP3',
    'PP4': 'PP4',
    'PP5': 'PP5',
    'PM5-': 'PM5_Supporting',
    'PS1--': 'PS1_Supporting',
    'BA1': 'BA1',
    'BA1-': "BA1_Ignore",
    'BP6++': 'BP6_Stand_Alone',
    'BS1': 'BS1',
    'BS2': 'BS2',
    'BS3': 'BS3',
    'BS4': 'BS4',
    'BP6+': 'BP6_Moderate',
    'BP1': 'BP1',
    'BP2': 'BP2',
    'BP3': 'BP3',
    'BP4': 'BP4',
    'BP5': 'BP5',
    'BP6': 'BP6',
    'BP7': 'BP7',
}

EVIDENCE_CODE_REVERSE_MAPPING = {
    'PVS1': 'verystrong_pathogenic',
    'PVS1_Strong': 'strong_pathogenic',
    'PVS1_Moderate': 'moderate_pathogenic',
    'PVS1_Supporting': 'supporting_pathogenic',
    'PP5+++': 'verystrong_pathogenic',
    'PS1': 'strong_pathogenic',
    'PS1_Moderate': 'moderate_pathogenic',
    'PS1_Supporting': 'supporting_pathogenic',
    'PS2': 'strong_pathogenic',
    'PS3': 'strong_pathogenic',
    'PS4': 'strong_pathogenic',
    'PP5++': 'strong_pathogenic',
    'PM5+': 'strong_pathogenic',
    'PM1_Strong': 'strong_pathogenic',
    'PM1': 'moderate_pathogenic',
    'PM1_Supporting': 'supporting_pathogenic',
    'PS1-': 'moderate_pathogenic',
    'PM2': 'moderate_pathogenic',
    'PM2_Supporting': 'supporting_pathogenic',
    'PM3': 'moderate_pathogenic',
    'PM3_Moderate': 'moderate_pathogenic',
    'PM3_Strong': 'strong_pathogenic',
    'PM3_Supporting': 'strong_pathogenic',
    'PM3_Very Strong': 'verystrong_pathogenic',
    'PM4': 'moderate_pathogenic',
    'PM4_Supporting': 'supporting_pathogenic',
    'PM5_Supporting': 'supporting_pathogenic',
    'PM5': 'moderate_pathogenic',
    'PM5_Strong': 'strong_pathogenic',
    'PM6': 'moderate_pathogenic',
    'PP5+': 'moderate_pathogenic',
    'PP1': 'supporting_pathogenic',
    'PP2': 'supporting_pathogenic',
    'PP3': 'supporting_pathogenic',
    'PP3_Moderate': 'moderate_pathogenic',
    'PP4': 'supporting_pathogenic',
    'PP4_Moderate': 'moderate_pathogenic',
    'PP4_Strong': 'strong_pathogenic',
    'PP5': 'supporting_pathogenic',
    'PM5-': 'supporting_pathogenic',
    'PS1--': 'supporting_pathogenic',
    'BA1': 'standalone_benign',
    'BA1-': "ignore",
    'BP6++': 'standalone_benign',
    'BS1_Stand Alone': 'standalone_benign',
    'BS1': 'strong_benign',
    'BS1_Supporting': 'supportive_benign',
    'BS2': 'strong_benign',
    'BS2_Stand Alone': 'standalone_benign',
    'BS2_Supporting': 'supportive_benign',
    'BS3': 'strong_benign',
    'BS4': 'strong_benign',
    'BP6+': 'strong_benign',
    'BP1': 'supportive_benign',
    'BP2': 'supportive_benign',
    'BP3': 'supportive_benign',
    'BP4': 'supportive_benign',
    'BP5': 'supportive_benign',
    'BP6': 'supportive_benign',
    'BP7': 'supportive_benign',
    'PVSX': 'verystrong_pathogenic',
    'PSX': 'strong_pathogenic',
    'PMX': 'moderate_pathogenic',
    'PPX': 'supporting_pathogenic',
    'BAX': 'standalone_benign',
    'BSX': 'strong_benign',
    'BPX': 'supportive_benign',
    'VAX': 'standalone_vus',
    'IGX': 'ignore'
}

DISEASE_SPECIFIC_CODES = {
    "PP1",
    "PS2",
    "PS3",
    "PS4",  # GWAS???
    "PM3",
    "PM6",
    "PP4",
    "BS3",
    "BS4",
    "BP2",
    "BP5",
    "BP5",
    ""
}

pathogenicity_mapping = {
    "pathogenic": "P",
    "likely pathogenic": "LP",
    "benign": "B",
    "likely benign": "LB",
    "uncertain significance": "VUS",
}


def get_variant_identifer(identifiers):
    for ident in identifiers:
        chromosome_id, value = ident.split(":")
        if chromosome_id not in CHROMOSOME_MAPPING:
            continue
        expression = value.split("g.")[-1]
        substitution = re.split("[0-9]+", expression)[-1]
        position = expression.replace(substitution, "")
        if ">" not in substitution:
            continue
        ref, alt = substitution.split(">")
        return CHROMOSOME_MAPPING[chromosome_id], position, ref, alt
    return None


def read_clingen_json_data(indata, number=None):
    clingen_dict = {}
    counter = 0
    for entry in indata["data"]:
        identifier = get_variant_identifer(entry["hgvs_ids"])
        evidence_codes = entry["evidence_codes"]
        pathogenicity = pathogenicity_mapping[entry["pathogenicity"].lower()]
        if identifier is None:
            continue
        counter += 1
        chromosome, position, ref, alt = identifier
        clingen_dict["-".join([chromosome, position, ref, alt])] = {
            "hgvs_ids": entry["hgvs_ids"],
            "pathogenicity": pathogenicity,
            "evidence_codes": evidence_codes,
        }
        if number is not None and counter >= number:
            break
    return clingen_dict


def read_franklin_data(infile):
    PATHOGENICITY = {"PATHOGENIC":"P", "LIKELY_PATHOGENIC": "LP","POSSIBLY_PATHOGENIC_MODERATE": "LP", "UNCERTAIN_SIGNIFICANCE":"VUS",
                     "POSSIBLY_PATHOGENIC_LOW": "LP","POSSIBLY_BENIGN":"LB","BENIGN":"B", "LIKELY_BENIGN":"LB"}
    data = {}
    counter = 0
    with open(infile) as fh:
        for line in fh:
            counter += 1
            if counter == 1:
                continue
            ls = line.split("\t")
            chromosome, pos, ref, alt = ls[2], ls[4], ls[5], ls[6]
            pathogenicity, evidences = PATHOGENICITY[ls[63]], ls[64]
            data[(chromosome, pos, ref, alt)] = (pathogenicity,evidences)
    return data




def create_comparison_file(clingen_json, annotation_json, franklin_txt, outfile):
    outfh = open(outfile, "w")
    print(
        "#identifier", "hash_id", "gene", "transcript_id", "hgvsc", "seq-pathogenicity", "clingen-pathogenicity",
        "franklin_pathogenicity"
        "seq-evidence_codes", "clingen-evidence-codes", "franklin_evidence_codes", "matched", "clingen_hgvs_ids",
        sep="\t", file=outfh
    )
    clingen_dict = read_clingen_json_data(clingen_json)
    annotation_data = load_json(annotation_json)
    franklin_data= read_franklin_data(franklin_txt)
    missing_in_franklin = 0
    for annotation in annotation_data:
        variant_info = annotation["variant"]
        chromosome, position, ref, alt = (
            variant_info["chromosome"],
            variant_info["position"],
            variant_info["ref"],
            variant_info["alt"],
        )
        franklin_pat = franklin_data[(chromosome, str(position), ref, alt)][0] if (chromosome, str(position), ref,
                                                                                   alt) in franklin_data else ""
        franklin_evidences = franklin_data[(chromosome, str(position), ref, alt)][1] if (chromosome, str(position), ref,
                                                                                         alt) in franklin_data else ""
        missing_in_franklin+=int(franklin_pat!="")
        hash_id = str(get_variant_uuid(chromosome, int(position), ref, alt, "hg19"))
        identifier = f"{chromosome}-{position}-{ref}-{alt}"
        clingen_data = clingen_dict[identifier]
        pathogenicity_summaries = annotation["annotations"]["summary"]
        found = False
        for sub_annotation in pathogenicity_summaries:
            gene_symbol, transcript_id, hgvsc, hgvsp, pathogenicity = (
                sub_annotation["gene_symbol"],
                sub_annotation["transcript_id"],
                sub_annotation["hgvsc"],
                sub_annotation["hgvsp"],
                sub_annotation["pathogenicity"]
            )
            if hgvsc not in clingen_data["hgvs_ids"] and hgvsp not in clingen_data["hgvs_ids"]:
                continue
            found = True

            print(
                identifier,
                hash_id,
                gene_symbol,
                transcript_id,
                hgvsc,
                pathogenicity["acmg_class"],
                clingen_data["pathogenicity"],
                franklin_pat,
                ",".join(pathogenicity["evidence_codes"]),
                ",".join(clingen_data["evidence_codes"]),
                franklin_evidences,
                "1",
                ",".join(clingen_data["hgvs_ids"]),
                sep="\t", file=outfh
            )
            break
        if not found:
            print(
                identifier,
                hash_id,
                pathogenicity_summaries[0]["gene_symbol"],
                pathogenicity_summaries[0]["transcript_id"],
                pathogenicity_summaries[0]["hgvsc"],
                pathogenicity_summaries[0]["pathogenicity"]["acmg_class"],
                clingen_data["pathogenicity"],
                franklin_pat,
                ",".join(pathogenicity_summaries[0]["pathogenicity"]["evidence_codes"]),
                ",".join(clingen_data["evidence_codes"]),
                franklin_evidences,
                "0",
                ",".join(clingen_data["hgvs_ids"]),
                sep="\t", file=outfh
            )
    outfh.close()

def myfunction(mystring):
    print(mystring)


"""

if __name__ == '__main__':
    globals()[sys.argv[1]](sys.argv[2])

"""