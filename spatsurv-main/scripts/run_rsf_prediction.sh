#!/bin/bash

source activate survival

python rsf_prediction.py --cv-folds patient --output-dir ./output/patient_exnormal_labeled_samples --features biomarker_region biomarker_cell celltype_prop neighbor_mat ripley --exclude-normal
