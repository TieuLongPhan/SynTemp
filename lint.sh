#!/bin/bash

# Run flake8 with specified rules
flake8 . --count --max-complexity=13 --max-line-length=100 \
    --exclude='./Docs' \
    --per-file-ignores="__init__.py:F401,rule_engine.py:F401,F403,F405,hier_engine.py:F401,F403,F405, rule_cluster.py:E203, graphormer_wrapper.py:E203, local_mapper_wrapper.py:E203, rdt_wrapper.py:E203, rxn_mapper_wrapper.py:E203, test_its_extraction.py:E501, test_rxn_mapper_wrapper.py:E501, test_rdt_wrapper.py:E501, test_local_mapper_wrapper.py:E501, test_graphormer_wrapper.py:E501, test_aam_validator.py:E501, test_aam_postprocess.py:E501, chemical_reaction_visualizer.py:E501" \
    --statistics
