# Author: Alejandro Santos Mejías
# Date last update: 2026-06-11
# Script name: main.R
# Aims: Main script to orchestrate all steps


rm(list = ls())
gc()

source("98_packages.R")

step_0 <- F
step_1 <- T
step_1b <- F

if(step_0){
  source("p_step/step_0_null_model.R", local = new.env())
}

if(step_1){
  source("p_step/step_1_run_communities.R", local = new.env())
}

if(step_1b){
  source("p_step/step_1b_consensus_communities.R", local = new.env())
}
