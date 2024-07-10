#!/bin/bash

# Run flake8 with specified rules
flake8 . --count --max-complexity=13 --max-line-length=100 \
    --per-file-ignores="__init__.py:F401,rule_engine.py:F401,F403,F405,hier_engine.py:F401,F403,F405" \
    --statistics
