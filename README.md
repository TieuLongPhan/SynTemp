# SynTemp
Graph based Reaction Template Extraction 

## Overview
This repository is dedicated to the systematic extraction of reaction rules from chemical processes. The primary focus lies in the computational analysis and transformation of molecular reactions into a structured set of rules, facilitating a deeper understanding of reaction mechanisms and pathways.

![screenshot](./Docs/Image/template_extraction_example.png)

### Step 1: Reaction Representation
The initial step involves the representation of chemical reactions, as illustrated below:
![Chemical Reactions](./Docs/Image/reactions.png)

In this stage, we detail the reactants, products, and the overall reaction scheme, laying the foundation for subsequent computational analysis.

### Step 2: Atom-Atom Mapping (AAM) and ITS Graph Generation
Utilizing the Atom-Atom Mapping (AAM) technique, we identify and map corresponding atoms across reactants and products. This mapping is crucial for accurately delineating the Imaginary Transitional State (ITS) of the reaction. The ITS graph, derived from AAM, represents the transitional phase of the reaction, highlighting key changes and interactions at the atomic level.
![Imaginary Transitional State Graph](./Docs/Image/graph_its.png)

### Step 3: Rules Extraction
In the final step, we focus on extracting the underlying rules from the ITS graph. This involves identifying significant nodes and edges that represent the core transformation mechanisms within the reaction. The extracted rules provide a concise and quantifiable description of the reaction process, instrumental for further analysis and application in computational chemistry.
![Extracted Rules](./Docs/Image/rules.png)

### Step 4: Rules Application
![DPO](./Docs/Image/dpo_example.png)


## Table of Contents
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)


## Repository Structure

SynTemp is organized into several key components, each dedicated to a specific aspect of chemical data processing:


## Installation

To install and set up the SynTemp framework, follow these steps. Please ensure you have Python 3.9 or later installed on your system.

### Prerequisites

- Python 3.11
- RDKit 2023.9.5
- networkx 3.2.1


### Step-by-Step Installation Guide

1. **Python Installation:**
  Ensure that Python 3.11 or later is installed on your system. You can download it from [python.org](https://www.python.org/downloads/).

2. **Creating a Virtual Environment (Optional but Recommended):**
  It's recommended to use a virtual environment to avoid conflicts with other projects or system-wide packages. Use the following commands to create and activate a virtual environment:

  ```bash
  python -m venv syntemp-env
  source syntemp-env/bin/activate  # On Windows use `syntemp-env\Scripts\activate`
  ```
  Or Conda

  ```bash
  conda create --name syntemp-env python=3.11
  conda activate syntemp-env
  ```

3. **Cloning and Installing SynTemp:**
  Clone the SynTemp repository from GitHub and install it:

  ```bash
  git clone https://github.com/TieuLongPhan/SynTemp.git
  cd SynTemp
  pip install -r requirements.txt
  ```

4. **Verify Installation:**
  After installation, you can verify that SynTemp is correctly installed by running a simple test or checking the package version.

  ```python
  python -c "import SynTemp; print(SynTemp.__version__)"
  ```

## Usage

### Use in script
  ```python
  from SynTemp.auto_template import AutoTemp

  smiles = (
      "COC(=O)[C@H](CCCCNC(=O)OCc1ccccc1)NC(=O)Nc1cc(OC)cc(C(C)(C)C)c1O>>"
      + "COC(=O)[C@H](CCCCN)NC(=O)Nc1cc(OC)cc(C(C)(C)C)c1O"
  )

  data = [{'R-id': '1', 'reactions': smiles}]

  auto = AutoTemp(
      rebalancing=True,
      mapper_types=["rxn_mapper", "graphormer", "local_mapper"],
      id="R-id",
      rsmi="reactions",
      n_jobs=1,
      verbose=2,
      batch_size=1,
      job_timeout=None,
      safe_mode=False,
      save_dir=None,
      fix_hydrogen=True,
      refinement_its=False,
  )

  (gml_rules, reaction_dicts, templates, hier_templates,
  its_correct, uncertain_hydrogen,) = auto.temp_extract(data, lib_path=None)

  print(gml_rules[0][0])
  >> '''rule [
   ruleID "0"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "N" ]
      node [ id 2 label "C" ]
      node [ id 3 label "O" ]
      node [ id 4 label "H" ]
   ]
   right [
      edge [ source 1 target 4 label "-" ]
      edge [ source 2 target 3 label "-" ]
   ]
]'''
  ```
  

### Use in command line
  ```bash
  echo -e "R-id,reaction\n0,COC(=O)[C@H](CCCCNC(=O)OCc1ccccc1)NC(=O)Nc1cc(OC)cc(C(C)(C)C)c1O>>COC(=O)[C@H](CCCCN)NC(=O)Nc1cc(OC)cc(C(C)(C)C)c1O" > test.csv
  python -m SynTemp --data_path test.csv --rebalancing --id 'R-id' --rsmi 'reaction' --rerun_aam --fix_hydrogen --log log.txt --save_dir ./
  ```

### Reproduce templates extraction
  Run these commands from the root of the cloned repository.
  ```bash
  python -m syntemp --data_path Data/USPTO_50K_original.csv --log_file Data/Test/log.txt --save_dir Data/Test/ --rebalancing --fix_hydrogen --rerun_aam --n_jobs 3 --batch_size 1000 --rsmi reactions --id ID
  ```
    
## Publication

[SynTemp: Efficient Extraction of Graph-Based Reaction Rules from Large-Scale Reaction Databases](https://chemrxiv.org/engage/chemrxiv/article-details/66f677b751558a15ef4cf5f7)

1. . . . 2024; doi:10.26434/chemrxiv-2024-tkm36  This content is a preprint and has not been peer-reviewed.

### Citation
```
@Article{Phan2024,
  author={Phan T-L, Weinbauer K, Gonzalez Laffitte ME, Pan Y, Merkle D, Andersen JL, et al},
  title={SynTemp: Efficient Extraction of Graph-Based Reaction Rules from Large-Scale Reaction Databases},
  journal={ChemRxiv},
  year={2024},
  doi={10.26434/chemrxiv-2024-tkm36},
  url={https://chemrxiv.org/engage/chemrxiv/article-details/66f677b751558a15ef4cf5f7}
}
```


## Contributing
- [Tieu-Long Phan](https://tieulongphan.github.io/)


## License

This project is licensed under MIT License - see the [License](LICENSE) file for details.

## Acknowledgments

This project has received funding from the European Unions Horizon Europe Doctoral Network programme under the Marie-Skłodowska-Curie grant agreement No 101072930 (TACsy -- Training Alliance for Computational)