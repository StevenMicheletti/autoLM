#!/bin/bash

OUTPUT=autolm.out
FREQ_FILE=all_freq.txt
ENV_FILE=env.txt
COORDS_FILE=locations.txt #remove if 
STANDARDIZE=TRUE

Rscript autoLM.R $OUTPUT $ENV_FILE $FREQ_FILE $STANDARDIZE $COORDS_FILE
