#!/bin/bash

python -W ignore ../../plot_SBOL_designs.py -params plot_parameters.csv -parts part_information.csv -designs dna_designs_riboj.csv -output cello_gates_riboj.pdf -regulation reg_information.csv
