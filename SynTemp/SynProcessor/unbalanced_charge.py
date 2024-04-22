from typing import Dict, Any, List
from joblib import Parallel, delayed


class UnbalancedCharge:
    @staticmethod
    def fix_negative_charge(
        reaction_dict: Dict[str, any],
        charges_column: str = "total_charge_in_products",
        id_column: str = "R-id",
        reaction_column: str = "reactions",
    ) -> Dict[str, any]:
        """
        Adjusts a reaction dictionary to compensate for a negative charge in the products by adding [Na+] ions.

        This function calculates the number of sodium ions ([Na+]) needed to neutralize negative charges in the reaction products.
        It then adds the appropriate number of sodium ions to both the reactants and products.

        Args:
            reaction_dict (Dict[str, any]): A dictionary representing a chemical reaction. Must include keys for
                                            'total_charge_in_products', 'reactants', 'products', 'R-id', and 'label'.

        Returns:
            Dict[str, any]: A new reaction dictionary with adjusted reactants and products to neutralize the negative charge.
                            The 'total_charge_in_products' is set to 0, assuming the charge has been neutralized.
        """
        # Calculate the number of sodium ions to add based on the absolute value of total_charge_in_products
        num_na_to_add = abs(reaction_dict[charges_column])
        sodium_ion = "[Na+]"

        # Generate the string to add, with the correct number of sodium ions
        sodium_addition = (
            "." + ".".join([sodium_ion] * num_na_to_add) if num_na_to_add > 0 else ""
        )

        # Add the sodium ions to reactants and products
        new_reactants = reaction_dict["reactants"] + sodium_addition
        new_products = reaction_dict["products"] + sodium_addition

        # Generate the new reaction string
        new_reactions = new_reactants + ">>" + new_products

        # Create the new reaction dictionary
        new_reaction_dict = {
            id_column: reaction_dict["R-id"],
            reaction_column: new_reactions,
            "label": reaction_dict["label"],
            "reactants": new_reactants,
            "products": new_products,
            charges_column: 0,  # Assuming the charge is neutralized
        }

        return new_reaction_dict

    @staticmethod
    def fix_positive_charge(
        reaction_dict: Dict[str, any],
        charges_column: str = "total_charge_in_products",
        id_column: str = "R-id",
        reaction_column: str = "reactions",
    ) -> Dict[str, any]:
        """
        Adjusts a reaction dictionary to compensate for a positive charge in the products by adding [Cl-] ions. The function
        takes into account the total positive charge indicated in the reaction dictionary and adds an equivalent number of
        chloride ions ([Cl-]) to both reactants and products to neutralize the charge.

        Args:
            reaction_dict (Dict[str, any]): A dictionary representing a chemical reaction. This dictionary must include
                                            keys for reactants, products, and a specified charge column (default is
                                            'total_charge_in_products') which contains the total charge of the products.
            charges_column (str, optional): The key in `reaction_dict` that contains the total charge of the products.
                                            Defaults to 'total_charge_in_products'.

        Returns:
            Dict[str, any]: A modified reaction dictionary with added [Cl-] ions to neutralize the positive charge. The
                            'total_charge_in_products' is updated to 0, indicating that the reaction's charge has been
                            neutralized. The dictionary includes updated 'reactants', 'products', and a new reaction string.
        """
        # Calculate the number of chloride ions to add based on the positive total charge in the specified column
        num_cl_to_add = abs(reaction_dict[charges_column])
        chloride_ion = "[Cl-]"

        # Generate the string to add, with the correct number of chloride ions
        chloride_addition = (
            "." + ".".join([chloride_ion] * num_cl_to_add) if num_cl_to_add > 0 else ""
        )

        # Add the chloride ions to reactants and products
        new_reactants = reaction_dict["reactants"] + chloride_addition
        new_products = reaction_dict["products"] + chloride_addition

        # Generate the new reaction string
        new_reactions = new_reactants + ">>" + new_products

        # Create and return the new reaction dictionary with the neutralized charge
        new_reaction_dict = {
            "R-id": reaction_dict[id_column],
            "new_reaction": new_reactions,
            "label": reaction_dict["label"],
            "reactants": new_reactants,
            "products": new_products,
            charges_column: 0,  # Update the charge column to reflect the neutralized charge
        }

        return new_reaction_dict

    @staticmethod
    def fix_unbalanced_charged(
        reaction_dict: Dict[str, any], charges_column: str = "total_charge_in_products"
    ) -> Dict[str, any]:
        """
        Adjusts a reaction dictionary to compensate for an unbalanced charge in the products by adding either [Cl-] ions
        for a positive charge or [Na+] ions for a negative charge. The function determines the direction of the charge
        imbalance using the specified charges column and applies the appropriate correction.

        Args:
            reaction_dict (Dict[str, any]): A dictionary representing a chemical reaction. This dictionary must include
                                            keys for reactants, products, and a specified charge column which contains
                                            the total charge of the products.
            charges_column (str, optional): The key in `reaction_dict` that contains the total charge of the products.
                                            Defaults to 'total_charge_in_products'.

        Returns:
            Dict[str, any]: A modified reaction dictionary with added ions to neutralize the charge imbalance. The
                            returned dictionary will have its charge neutralized and include updated 'reactants',
                            'products', and a new reaction string. The specific ions added ([Cl-] for positive charges
                            or [Na+] for negative charges) depend on the initial charge imbalance.
        """
        if reaction_dict[charges_column] > 0:
            return UnbalancedCharge.fix_positive_charge(reaction_dict, charges_column)
        elif reaction_dict[charges_column] < 0:
            return UnbalancedCharge.fix_negative_charge(reaction_dict, charges_column)
        else:
            # If the charge is already balanced, return the original dictionary without modification.
            return reaction_dict

    @classmethod
    def parallel_fix_unbalanced_charge(
        cls,
        reaction_dicts: List[Dict[str, Any]],
        charges_column: str = "total_charge_in_products",
        n_jobs: int = -1,
    ) -> List[Dict[str, Any]]:
        """
        Processes a list of reaction dictionaries in parallel to compensate for unbalanced charges in the products,
        adding either [Cl-] ions for positive charges or [Na+] ions for negative charges.

        Args:
            reaction_dicts (List[Dict[str, Any]]): A list of dictionaries, each representing a chemical reaction that may have an unbalanced charge.
            charges_column (str): The key in each reaction dictionary that contains the total charge of the products. Defaults to 'total_charge_in_products'.
            n_jobs (int): The number of CPU cores to use for parallel processing. -1 means using all available cores.

        Returns:
            List[Dict[str, Any]]: A list of modified reaction dictionaries with charges neutralized, reflecting the addition of necessary ions.

        Note:
            This function requires the joblib library for parallel execution. Ensure joblib is installed and available for import.
        """
        # Use joblib.Parallel and joblib.delayed to parallelize the charge fixing
        fixed_reactions = Parallel(n_jobs=n_jobs)(
            delayed(cls.fix_unbalanced_charged)(reaction_dict, charges_column)
            for reaction_dict in reaction_dicts
        )
        return fixed_reactions
