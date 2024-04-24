import pandas as pd
from rdkit import Chem
from rdkit.Chem import MolStandardize
from joblib import Parallel, delayed
from typing import Tuple, Optional, Union
from typing import List

from SynTemp.SynUtils.chemutils import (
    normalize_molecule,
    canonicalize_tautomer,
    salts_remover,
    reionize_charges,
    uncharge_molecule,
    assign_stereochemistry,
    fragments_remover,
    remove_hydrogens_and_sanitize,
)


class Standardize:
    """
    The SMILESStandardizer class is designed for comprehensive standardization of chemical structures
    represented in SMILES (Simplified Molecular Input Line Entry System) format.
    This class utilizes various functionalities from the RDKit library to process and normalize chemical structures
    for consistency and comparability in cheminformatics applications.
    """

    def __init__(self):
        self.normalizer = MolStandardize.normalize.Normalizer()

    @staticmethod
    def standardize_mol(
        mol: Chem.Mol,
        normalize: bool = True,
        tautomerize: bool = False,
        remove_salts: bool = False,
        handle_charges: bool = False,
        uncharge: bool = True,
        handle_stereo: bool = True,
        remove_fragments: bool = False,
        largest_fragment_only: bool = False,
    ) -> Chem.Mol:
        """
        Apply comprehensive standardization to a molecule.

        Parameters:
        mol (Chem.Mol): RDKit Mol object to be standardized.
        normalize (bool, optional): Perform normalization (corrects functional groups and recharges).
        tautomerize (bool, optional): Canonicalize tautomers.
        remove_salts (bool, optional): Remove salt fragments from the molecule.
        handle_charges (bool, optional): Adjust molecule to its most likely ionic state using Reionizer.
        uncharge (bool, optional): Neutralize molecule by removing counter-ions using Uncharger.
        handle_stereo (bool, optional): Handle stereochemistry.
        remove_fragments (bool, optional): Remove small fragments, keeping only the largest one.
        largest_fragment_only (bool, optional): Keep only the largest fragment in the molecule.

        Returns:
        Chem.Mol: Standardized RDKit Mol object.
        """
        if mol is None:
            return None

        # Ensure ring information is computed
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)

        # Apply standardization steps based on the specified options
        if normalize:
            mol = normalize_molecule(mol)
        if tautomerize:
            try:
                mol = canonicalize_tautomer(mol)
            except ValueError:
                mol = mol
        if remove_salts:
            mol = salts_remover(mol)
        if handle_charges:
            mol = reionize_charges(mol)
        if uncharge:
            mol = uncharge_molecule(mol)
        if handle_stereo:
            assign_stereochemistry(mol, cleanIt=True, force=True)
        if remove_fragments or largest_fragment_only:
            mol = fragments_remover(mol)

        # Remove explicit hydrogens and sanitize the molecule
        mol = remove_hydrogens_and_sanitize(mol)

        return mol

    @staticmethod
    def standardize_smiles(smiles: str, **kwargs) -> Tuple[str, Optional[Chem.Mol]]:
        """
        Standardize a SMILES string.

        Parameters
        ----------
        smiles : str
            SMILES string to be standardized.

        Returns
        -------
        Tuple[str, Optional[Chem.Mol]]
            A tuple containing the standardized SMILES string and the standardized RDKit Mol object.
        """
        original_mol = Chem.MolFromSmiles(smiles)
        if not original_mol:
            return None, None
        try:
            standardized_mol = Standardize.standardize_mol(original_mol, **kwargs)
            standardized_smiles = Chem.MolToSmiles(standardized_mol)
            return standardized_smiles, standardized_mol
        except Chem.MolSanitizeException:
            return "Sanitization failed for SMILES: " + smiles, None

    def standardize_dict_smiles(
        self,
        data_input: Union[pd.DataFrame, List[dict]],
        keys: List[str],
        parallel: bool = True,
        n_jobs: int = 4,
        keep_mol: bool = True,
        **kwargs
    ) -> Union[pd.DataFrame, List[dict]]:
        """
        Standardize SMILES strings in a pandas DataFrame or a list of dictionaries for multiple keys.

        Args:
            data_input (Union[pd.DataFrame, List[dict]]): Data containing SMILES strings to be standardized.
            keys (List[str]): Keys or column names for SMILES strings in the data.
            visualize (bool, optional): If True, visualize the molecules during standardization.
            parallel (bool, optional): If True and data_input is a list of dicts, use parallel processing.
            n_jobs (int, optional): Number of jobs to run in parallel (if parallel is True).
            keep_mol (bool, optional): If True, keep the RDKit Mol objects in the output.

        Returns:
            Union[pd.DataFrame, List[dict]]: Data with standardized SMILES strings and, optionally,
            standardized RDKit Mol objects for each key.
        """
        if isinstance(data_input, pd.DataFrame):
            for key in keys:
                if parallel:
                    with Parallel(n_jobs=n_jobs, verbose=1) as parallel:
                        standardized_results = parallel(
                            delayed(self.standardize_smiles)(smiles, **kwargs)
                            for smiles in data_input[key]
                        )
                    data_input["standardized_" + key] = [
                        result[0] for result in standardized_results
                    ]
                    if keep_mol:
                        data_input["standardized_mol_" + key] = [
                            result[1] for result in standardized_results
                        ]
                else:
                    data_input["standardized_" + key] = data_input[key].apply(
                        lambda x: self.standardize_smiles(x, **kwargs)[0]
                    )
                    if keep_mol:
                        data_input["standardized_mol_" + key] = data_input[key].apply(
                            lambda x: self.standardize_smiles(x, **kwargs)[1]
                        )
        elif isinstance(data_input, list) and all(
            isinstance(item, dict) for item in data_input
        ):
            for key in keys:
                if parallel:
                    standardized_results = Parallel(n_jobs=n_jobs, verbose=1)(
                        delayed(self.standardize_smiles)(
                            reaction_data.get(key, ""), **kwargs
                        )
                        for reaction_data in data_input
                    )
                    for i, reaction_data in enumerate(data_input):
                        reaction_data["standardized_" + key] = standardized_results[i][
                            0
                        ]
                        if keep_mol:
                            reaction_data["standardized_mol_" + key] = (
                                standardized_results[i][1]
                            )
                else:
                    for reaction_data in data_input:
                        smiles = reaction_data.get(key, "")
                        standardized_smiles, standardized_mol = self.standardize_smiles(
                            smiles, **kwargs
                        )
                        reaction_data["standardized_" + key] = standardized_smiles
                        if keep_mol:
                            reaction_data["standardized_mol_" + key] = standardized_mol
        else:
            raise TypeError(
                "Input must be either a pandas DataFrame or a list of dictionaries."
            )

        return data_input
