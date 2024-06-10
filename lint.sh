#!/bin/bash

flake8 . --count --max-complexity=13 --max-line-length=120 \
	--per-file-ignores="__init__.py:F401" \
	--per-file-ignores="rule_engine.py:F403" \
	--statistics