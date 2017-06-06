PLOT_TO_FILE <- TRUE
FOLDER <- "results"

if(!file.exists(FOLDER))
    dir.create(FOLDER)

ABOVE <- FALSE

SPECIES <- "Pleuronectes platessa"
MINDSTEMAAL <- 26
source("model.R")

SPECIES <- "Gadus morhua"
MINDSTEMAAL <- 30
source("model.R")

ABOVE <- TRUE

SPECIES <- "Pleuronectes platessa"
MINDSTEMAAL <- 26
source("model.R")

SPECIES <- "Gadus morhua"
MINDSTEMAAL <- 30
source("model.R")
