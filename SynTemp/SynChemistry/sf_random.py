from typing import List, Optional
import random
import copy

class SFRandom:
    def sort_reactions(self, reactions: List[str], seed: Optional[int] = 42) -> List[str]:
        """
        Shuffle a list of reaction strings randomly using a seed for reproducibility and return the new order.

        Parameters:
        reactions (List[str]): A list of reactions represented as strings.
        seed (Optional[int]): An optional integer seed for random number generator to ensure reproducibility.

        Returns:
        List[str]: A list containing the reactions in a randomly shuffled order.
        """
        if seed is not None:
            random.seed(seed)  # Set the random seed for reproducibility

        shuffled_reactions = copy.deepcopy(reactions)  # Create a copy to shuffle
        random.shuffle(shuffled_reactions)
        return shuffled_reactions