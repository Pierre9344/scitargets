library(targets)
library(tarchetypes)
library(scitargets)
# library(crew)
# library(autometric)

# tar_source() # Source all the R scripts in the ./R/ directory in the environment. Can be used to load libraries and functions when expending the pipeline.

###################################################
###  PIPELINE CONTROLLER FOR PARALLEL COMPUTING ###
###################################################
# my_controller <- crew_controller_local(
#   name = "my_pc_local",
#   workers = 2,
#   options_local = crew_options_local(log_directory = "./logs/crew/")
# )
#
# if (tar_active()) {
#   my_controller$start()
#   log_start(
#     path = "./logs/autometric.txt",
#     seconds = 1
#   )
# }

##########################
###  PIPELINE OPTIONS  ###
##########################

tar_option_set(
  tidy_eval = T,
  # controller = my_controller, # parallel processing if the targets need it
  seed = 5114, # main seed (step specific seed are derived from it and step name)
  format = "qs", # more efficiant than default rds
  memory = "transient", # unload the targets after the steps are completed
  # the next argument indicate the packages that need to be loaded by the pipeline.This override the packages loaded in the main environment
  packages = c(
    "scitargets"
  )
)

##################
###  PIPELINE  ###
##################
list(
  scitargets::tar_demultiplex_hto(
    run_id = "",
    project_id = "",
    min_genes_detected = 200,
    max_genes_detected = 4000,
    mt_percent_cutoff = 5,
    singlets_dim_to_use = 1:15
  ),

  tarchetypes::tar_quarto(
    name = report,
    path = "."
    )
)
