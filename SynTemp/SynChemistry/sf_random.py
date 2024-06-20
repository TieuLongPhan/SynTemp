from typing import List, Optional
import random
import copy


class SFRandom:
    def __init__(self, seed: Optional[int] = 42) -> None:
        """
        Initialize the randomizer with an optional seed to ensure reproducibility.

        Parameters:
        seed (Optional[int]): The seed for the random number generator.
        """
        self.seed = seed
        if self.seed is not None:
            random.seed(
                self.seed
            )  # Set the random seed for reproducibility once upon initialization

    def sort_reactions(self, reactions: List[str]) -> List[str]:
        """
        Shuffle a list of reaction strings randomly using the seed provided at initialization for reproducibility,
        and return the new order.

        Parameters:
        reactions (List[str]): A list of reactions represented as strings.

        Returns:
        List[str]: A list containing the reactions in a randomly shuffled order.
        """
        shuffled_reactions = copy.deepcopy(reactions)  # Create a deep copy to shuffle
        random.shuffle(shuffled_reactions)
        return shuffled_reactions
