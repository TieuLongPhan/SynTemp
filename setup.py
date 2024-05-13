from setuptools import setup, find_packages

setup(
    name='SynTemp',
    version='0.0.1',
    description='Graph Rules Extraction from Reaction Graphs',
    author='TieuLongPhan',
    author_email='ptlong8995@gmail.com',
    url='https://github.com/TieuLongPhan/SynTemp',
    packages=find_packages(),
    install_requires=[
        'rdkit==2023.9.5',
        'networkx==3.2.1'
   
    ],
    python_requires='>=3.11',
   
)
