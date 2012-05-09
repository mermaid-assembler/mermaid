#!/bin/sh

export PYTHONPATH=$PYTHONPATH:$CWD
python tests/bwt_test.py
python tests/mtf_test.py
