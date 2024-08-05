from typing import Dict, Tuple
from CGRtools import SMILESRead
from CGRtools.files import RDFRead
import logging
import pickle


def parse_reactions(
    input_file: str, id_tag: str, do_basic_standardization: bool = True
) -> Dict:
    """
    Reads an RDF/SMILES file. Returns a dictionary where the key is the reaction ID, and the value is the
    ReactionContainer.
    :param do_basic_standardization: bool
    :param input_file: str
    :param id_tag: str
    :return: Dict
    """
    data = {}
    if input_file.endswith(".rdf"):
        with RDFRead(input_file) as f:
            while True:
                try:
                    r = next(f)
                    key = r.meta[id_tag]
                    if do_basic_standardization:
                        r.thiele()
                        r.standardize()
                    data[key] = r
                    logging.info(f"Reaction {key} passed..")
                except StopIteration:
                    break
    elif input_file.endswith(".smi") or input_file.endswith(".smiles"):
        with SMILESRead(
            input_file, ignore=True, store_log=True, remap=False, header=True
        ) as ifile, open(input_file) as meta_searcher:
            id_tag_position = meta_searcher.readline().strip().split().index(id_tag)
            if id_tag_position is None or id_tag_position == 0:
                logging.critical(f"No reaction ID tag was found in the header!")
                raise ValueError(f"No reaction ID tag was found in the header!")
            for reaction in ifile._data:
                if isinstance(reaction, tuple):
                    meta_searcher.seek(reaction.position)
                    line = meta_searcher.readline().strip().split()
                    if len(line) <= id_tag_position:
                        logging.critical(
                            f"No reaction ID tag was found in line {reaction.number}!"
                        )
                        raise ValueError(
                            f"No reaction ID tag was found in line {reaction.number}!"
                        )
                    r_id = line[id_tag_position]
                    logging.critical(
                        f"Reaction {r_id}: Parser has returned an error message\n{reaction.log}"
                    )
                    continue
                if do_basic_standardization:
                    reaction.thiele()
                    reaction.standardize()
                key = reaction.meta[id_tag]
                data[key] = reaction
                logging.info(f"Reaction {key} passed..")
    return data


def filter_data(referenced_aam: dict, generated_aam: dict) -> Tuple[Dict, Dict]:
    """
    To prevent any errors, only common reactions are selected, and others are removed.
    :param referenced_aam: dict
    :param generated_aam: dict
    :return: Tuple[Dict, Dict]
    """
    common_IDs = set(referenced_aam.keys()) & set(generated_aam.keys())
    new_ref, new_gen = {}, {}
    for key in common_IDs:
        new_ref[key] = referenced_aam[key]
        new_gen[key] = generated_aam[key]
    return new_ref, new_gen


def make_cgr(data: dict) -> Dict:
    """
    Transforms the reactions to CGRs.
    :param data: dict
    :return: Dict
    """
    result = {}
    for key, r in data.items():
        cgr = ~r
        result[key] = cgr

    return result


def __config_log(log_file):
    logging.basicConfig(
        filename=log_file,
        level=logging.DEBUG,
        filemode="w",
        format="%(asctime)s: %(message)s",
        datefmt="%d/%m/%Y %H:%M:%S",
    )


def main(
    reference_file: str,
    generated_file: str,
    log_file: str,
    id_tag: str,
    archive_file: str,
    ignore_basic_stand: bool,
):
    """
    Computes statistics on the comparison of the referenced and generated mappings.
    :param ignore_basic_stand: bool
    :param reference_file: str
    :param generated_file: str
    :param log_file: str
    :param archive_file: str
    :param id_tag: str
    """
    __config_log(log_file=log_file)

    print("Loading the reference dataset..")
    logging.info("Loading the reference dataset..")
    ref_mapping = parse_reactions(
        input_file=reference_file,
        id_tag=id_tag,
        do_basic_standardization=(not ignore_basic_stand),
    )
    print(f"{len(ref_mapping)} reactions were found..")

    print("Loading the generated dataset..")
    logging.info("Loading the generated dataset..")
    gen_mapping = parse_reactions(
        input_file=generated_file,
        id_tag=id_tag,
        do_basic_standardization=(not ignore_basic_stand),
    )
    print(f"{len(gen_mapping)} reactions were found..")

    print("Filtering data..")
    ref_mapping, gen_mapping = filter_data(
        referenced_aam=ref_mapping, generated_aam=gen_mapping
    )
    print(f"{len(ref_mapping)} reactions stayed..")

    print("Convert reference data to CGR..")
    ref_cgr = make_cgr(data=ref_mapping)

    print("Convert generated data to CGR..")
    gen_cgr = make_cgr(data=gen_mapping)

    statistics = dict(
        same_mapping=0, differentCD=0, differentCGRs=0, not_equal_reactions=0
    )
    reactions = dict(
        good_mapping=[],
        bad_mapping_diff_CD=[],
        bad_mapping_same_CD_but_diff_CGRs=[],
        not_equal_reactions=[],
    )

    for key in sorted(ref_cgr.keys()):
        print(key)
        if ref_mapping[key] != gen_mapping[key]:
            statistics["not_equal_reactions"] += 1
            reactions["not_equal_reactions"].append(key)
        elif ref_cgr[key] != gen_cgr[key]:
            chemical_distance = abs(
                len(ref_cgr[key].center_atoms)
                + len(ref_cgr[key].center_bonds)
                - len(gen_cgr[key].center_atoms)
                - len(gen_cgr[key].center_bonds)
            )
            if chemical_distance > 0:
                statistics["differentCD"] += 1
                reactions["bad_mapping_diff_CD"].append(key)
            else:
                statistics["differentCGRs"] += 1
                reactions["bad_mapping_same_CD_but_diff_CGRs"].append(key)
        else:
            statistics["same_mapping"] += 1
            reactions["good_mapping"].append(key)
    if archive_file:
        pickle.dump(reactions, open(archive_file, "wb"))

    print("----------------------------------------------------------------------")
    print("**********************************************************************")
    print(f'Reactions mapped identically: {statistics["same_mapping"]}')
    print(f'Reactions mapped differently (CD>0): {statistics["differentCD"]}')
    print(
        f'Reactions mapped differently (CD=0 but CGRs are different): {statistics["differentCGRs"]}'
    )
    print(f'Not equal reactions: {statistics["not_equal_reactions"]}')
    print("**********************************************************************")
    print("----------------------------------------------------------------------")
    return reactions


if __name__ == "__main__":
    from pathlib import Path
    import sys

    root_dir = Path(__file__).parents[2]
    sys.path.append(str(root_dir))
    from syntemp.SynUtils.utils import load_database, save_database
    import pandas as pd

    df = pd.read_csv(f"{root_dir}/Data/AAM/cgrtool_benchmark/USPTO_3K.csv")
    df["GroundTruth"] = df["LocalMapper"]
    df["ID"] = df.index
    df[["GroundTruth", "ID"]].to_csv(
        f"{root_dir}/Data/AAM/cgrtool_benchmark/ref.smiles", index=False, sep=" "
    )
    df[["RXNMapper", "ID"]].to_csv(
        f"{root_dir}/Data/AAM/cgrtool_benchmark/rxnmapper.smiles", index=False, sep=" "
    )
    df[["GraphMapper", "ID"]].to_csv(
        f"{root_dir}/Data/AAM/cgrtool_benchmark/graphmapper.smiles",
        index=False,
        sep=" ",
    )
    df[["LocalMapper", "ID"]].to_csv(
        f"{root_dir}/Data/AAM/cgrtool_benchmark/localmapper.smiles",
        index=False,
        sep=" ",
    )
    df["CGRTool_rxnmapper"] = False
    df["CGRTool_graphmapper"] = False
    df["CGRTool_localmapper"] = False

    rxnmapper_check = main(
        reference_file=f"{root_dir}/Data/AAM/cgrtool_benchmark/ref.smiles",
        generated_file=f"{root_dir}/Data/AAM/cgrtool_benchmark/rxnmapper.smiles",
        log_file=None,
        id_tag="ID",
        archive_file=None,
        ignore_basic_stand=True,
    )

    graphmapper_check = main(
        reference_file=f"{root_dir}/Data/AAM/cgrtool_benchmark/ref.smiles",
        generated_file=f"{root_dir}/Data/AAM/cgrtool_benchmark/graphmapper.smiles",
        log_file=None,
        id_tag="ID",
        archive_file=None,
        ignore_basic_stand=True,
    )

    localmapper_check = main(
        reference_file=f"{root_dir}/Data/AAM/cgrtool_benchmark/ref.smiles",
        generated_file=f"{root_dir}/Data/AAM/cgrtool_benchmark//localmapper.smiles",
        log_file=None,
        id_tag="ID",
        archive_file=None,
        ignore_basic_stand=True,
    )

    for i in rxnmapper_check["good_mapping"]:
        df.loc[int(i), "CGRTool_rxnmapper"] = True

    for i in graphmapper_check["good_mapping"]:
        df.loc[int(i), "CGRTool_graphmapper"] = True

    for i in localmapper_check["good_mapping"]:
        df.loc[int(i), "CGRTool_localmapper"] = True

    df.to_csv(f"{root_dir}/Data/AAM/cgrtool_benchmark/uspto_3k_cgrtool_old.csv")
