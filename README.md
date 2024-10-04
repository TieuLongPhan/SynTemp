# SynTemp
Graph-based Reaction Templates/Rules Extraction 

## Overview

This repository is dedicated to the systematic extraction of reaction rules from reaction databases. Our primary focus is the computational analysis and transformation of molecular reactions into a structured set of rules. This work facilitates a deeper understanding of reaction mechanisms and pathways. The `SynTemp` framework is organized into four main phases:

1. **AAMs Inference**: Based on ensemble AAMs for accurate atom mapping.
2. **Imaginary Transition State (ITS) Completion**: Enhances ITS by incorporating hydrogen inference to fully capture the reaction mechanism.
3. **Reaction Center Detection and Extension**: Focuses on identifying and extending the core active sites of reactions.
4. **Hierarchical Clustering**: Groups extended reaction centers or partial ITS to analyze reaction patterns.

The general framework and its components are depicted in Figures A, B, and C below.

![screenshot](https://github.com/TieuLongPhan/SynTemp/raw/main/Docs/Image/TOC.png)

### Downstream Applications

- **Templates Analysis**: We have developed topological descriptors for ITS graphs to encapsulate the essential information of templates.
- **Rules Application**: Observes the trade-off between radii and coverage/novelty metrics. Increased coverage tends to reduce the number of output solutions due to the complexities of subgraph matching within the DPO framework. However, it also decreases novelty. This trade-off serves as a precursor to our forthcoming research, which will focus on developing a constrained framework for synthesis planning.


## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [Publication](#publication)
- [License](#license)
- [Acknowledgments](#acknowledgments)


## Installation

To install and set up the SynTemp framework, follow these steps. Please ensure you have Python 3.9 or later installed on your system.

### Prerequisites

- Python 3.11
- rdkit>=2023.9.5
- networkx>=3.3
- fgutils==0.0.17
- seaborn==0.13.2
- joblib==1.3.2

If you want to run ensemble AAMs

- dgl==2.1.0
- dgllife==0.3.2
- localmapper==0.1.3
- rxn-chem-utils==1.5.0
- rxn-utils==2.0.0
- rxnmapper==0.3.0
- chython==1.75
- chytorch==1.60
- chytorch-rxnmap==1.4
- torchdata==0.7.1


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

3. **Install from PyPi:**
  The easiest way to use SynTemp is by installing the PyPI package 
  [syntemp](https://pypi.org/project/syntemp/).

  ```
  pip install syntemp
  ```

4. **Verify Installation:**
  After installation, you can verify that Syn Temp is correctly installed by running a simple test

  ```bash
  echo -e "R-id,reaction\n0,COC(=O)[C@H](CCCCNC(=O)OCc1ccccc1)NC(=O)Nc1cc(OC)cc(C(C)(C)C)c1O>>COC(=O)[C@H](CCCCN)NC(=O)Nc1cc(OC)cc(C(C)(C)C)c1O" > test.csv
  python -m syntemp --data_path test.csv --rebalancing --id 'R-id' --rsmi 'reaction' --rerun_aam --fix_hydrogen --log_file ./log.txt --save_dir ./
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
  python -m syntemp --data_path test.csv --rebalancing --id 'R-id' --rsmi 'reaction' --rerun_aam --fix_hydrogen --log_file ./log.txt --save_dir ./
  ```

### Reproduce templates extraction
  Run these commands from the root of the cloned repository.
  ```bash
  python -m syntemp --data_path Data/USPTO_50K_original.csv --log_file Data/Test/log.txt --save_dir Data/Test/ --rebalancing --fix_hydrogen --rerun_aam --n_jobs 3 --batch_size 1000 --rsmi reactions --id ID
  ```
    
## Publication

[SynTemp: Efficient Extraction of Graph-Based Reaction Rules from Large-Scale Reaction Databases](https://chemrxiv.org/engage/chemrxiv/article-details/66f677b751558a15ef4cf5f7)


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

This project has received funding from the European Unions Horizon Europe Doctoral Network programme under the Marie-Skłodowska-Curie grant agreement No 101072930 ([TACsy](https://tacsy.eu/) -- Training Alliance for Computational)