#!/bin/bash

python -W ignore ../../plot_SBOL_designs.py -params plot_parameters.csv -parts part_information.csv -designs dna_designs.csv -output cello_eugene.pdf -regulation reg_information.csv
