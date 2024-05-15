import os
import shutil
from uuid import uuid4
from chython import smiles
from rxnmapper import RXNMapper
from contextlib import contextmanager
from localmapper import localmapper
import re
from rdkit import Chem

# Initialize RXNMapper instance
rxn_mapper = RXNMapper()


def map_with_rxn_mapper(reaction_smiles: str, rxn_mapper: RXNMapper) -> str:
    """
    Maps a reaction using the RXNMapper.

    Parameters:
        reaction_smiles (str): The SMILES string of the reaction to be mapped.
        rxn_mapper (RXNMapper): The instance of RXNMapper used for mapping.

    Returns:
        str: The mapped reaction SMILES string.
    """
    try:
        # Map reaction using RXNMapper
        mapped_rxn = rxn_mapper.get_attention_guided_atom_maps(
            [reaction_smiles], canonicalize_rxns=False
        )[0]["mapped_rxn"]
        # Return mapped reaction
        return mapped_rxn.split(" ")[0] if " " in mapped_rxn else mapped_rxn
    except Exception as e:
        print(f"RXNMapper mapping failed: {e}")
        return reaction_smiles


def map_with_graphormer(reaction_smiles: str) -> str:
    """
    Maps a reaction using the Graphormer.

    Parameters:
        reaction_smiles (str): The SMILES string of the reaction to be mapped.

    Returns:
        str: The mapped reaction SMILES string.
    """
    try:
        # Initialize Graphormer with reaction SMILES
        myrxn = smiles(reaction_smiles)
        myrxn.reset_mapping()
        # Get mapped reaction
        mapped_rxn = format(myrxn, "m")
        # Return mapped reaction
        return mapped_rxn.split(" ")[0] if " " in mapped_rxn else mapped_rxn
    except Exception as e:
        print(f"Graphormer mapping failed: {e}")
        return reaction_smiles


def map_with_local_mapper(reaction_smiles: str, mapper=localmapper()) -> str:
    """
    Maps a reaction using the AtomMapper.

    Parameters:
        reaction_smiles (str): The SMILES string of the reaction to be mapped.
    Returns:
        str: The mapped reaction SMILES string.
    """

    try:
        # Map reaction using AtomMapper
        result = mapper.get_atom_map(reaction_smiles)
        # Return mapped reaction
        return result
    except Exception as e:
        print(f"AtomMapper mapping failed: {e}")
        return reaction_smiles


@contextmanager
def temporary_change_dir(target_directory):
    """
    Context manager to temporarily change the working directory.
    """
    original_directory = os.getcwd()  # Save the original working directory
    try:
        os.chdir(target_directory)
        yield
    finally:
        os.chdir(original_directory)  # Revert back to the original directory


def map_with_rdt(reaction_smiles: str, rdt_jar_path: str, working_dir: str) -> str:
    """
    Maps a reaction using the RDT (Retro-Data) tool.

    Parameters:
        reaction_smiles (str): The SMILES string of the reaction to be mapped.
        rdt_jar_path (str): The file path to the RDT JAR file.
        working_dir (str): The working directory to execute RDT.

    Returns:
        str: The mapped reaction SMILES string.
    """
    unique_id = uuid4()
    unique_dir = os.path.join(working_dir, f"RDT_{unique_id}")
    output_file = "ECBLAST_smiles_AAM.txt"

    try:
        os.makedirs(unique_dir)

        with temporary_change_dir(unique_dir):
            command = (
                f'java -jar "{rdt_jar_path}" -Q SMI -q "{reaction_smiles}"'
                + f" -c -j AAM -f TEXT > {output_file}"
            )
            os.system(command)

            # Read the output from the file
            with open(output_file, "r") as inputFile:
                RDTaam = inputFile.read().splitlines()[3]

        return RDTaam.split(" ")[0] if " " in RDTaam else RDTaam
    except Exception as e:
        print(f"RDT mapping failed: {e}")
        return reaction_smiles
    finally:
        shutil.rmtree(unique_dir, ignore_errors=True)


def remove_atom_mapping(smiles: str) -> str:
    pattern = re.compile(r":\d+")
    smiles = pattern.sub("", smiles)
    pattern = re.compile(r"\[(?P<atom>(B|C|N|O|P|S|F|Cl|Br|I){1,2})(?:H\d?)?\]")
    smiles = pattern.sub(r"\g<atom>", smiles)
    return smiles


def normalize_smiles(smiles: str) -> str:
    if ">>" in smiles:
        return ">>".join([normalize_smiles(t) for t in smiles.split(">>")])
    elif "." in smiles:
        token = sorted(
            smiles.split("."),
            key=lambda x: (sum(1 for c in x if c.isupper()), sum(ord(c) for c in x)),
            reverse=True,
        )
        return ".".join([normalize_smiles(t) for t in token])
    else:
        return Chem.CanonSmiles(remove_atom_mapping(smiles))
