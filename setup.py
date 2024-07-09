from setuptools import setup, find_packages

setup(
    name="SynTemp",
    version="0.0.1",
    description="Graph Rules Extraction from Reaction Graphs",
    author="TieuLongPhan",
    author_email="ptlong8995@gmail.com",
    url="https://github.com/TieuLongPhan/SynTemp",
    packages=find_packages(),
    install_requires=["synrbl==0.0.22", "rdkit==2023.9.5", "networkx==3.2.1", "chython==1.75", "chytorch==1.60",
                      "chytorch-rxnmap==1.4", "dgl==2.0.0", "dgllife==0.3.2", "localmapper==0.1.3",
                      "rxn-chem-utils==1.5.0", "rxn-utils==2.0.0", "rxnmapper==0.3.0"],
    python_requires=">=3.11",
)
