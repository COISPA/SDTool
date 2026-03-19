

# 1) INSTALL the ZIP PACKAGE


# 2) Load package
library(SDTool)
# 3) set your working directory

setwd("C:/_LAVORO CORRENTE")

#4) define the folder where the results will be saved

out_dir="C:\\_LAVORO CORRENTE"


# 5) Load dummy datasets shipped with the package
data("cs_gsa99_dummy", package = "SDTool")
data("cl_gsa99_dummy", package = "SDTool")

head(cs_gsa99_dummy)

head(cl_gsa99_dummy)

# 6) Run functions to convert CS and CL in COST objects

costCS_from_CS <- convert_cs_to_costCS(cs_gsa99_dummy,write_csv = TRUE,out_dir=out_dir)

costCL_from_CL <- convert_cl_to_costCL(cl_gsa99_dummy,write_csv = TRUE,out_dir=out_dir)


# =========================================================
# 7) run data_investigation()
# =========================================================

# IMPORTANT:
# the function needs a metier shortcode table. An example is provided with the package
# the functions work on one species at a time

# Read from csv:
metier_shortcode <- read.csv(
    file.path( "C:\\_LAVORO CORRENTE\\input files_OPT\\metier_shortcode.csv"),
    sep = ";",
    stringsAsFactors = FALSE
)

available_species <- unique(costCS_from_CS@hl$spp)

print(available_species)

# Set target species manually after checking the printed species
target_species <- available_species[1]

# Run the function data_investigation and the results will be saved in the folder DataInvestigation, autamitically created in your out_dir directory.

res_investigation <- data_investigation(
    CS = costCS_from_CS,
    CL = costCL_from_CL,
    metier_shortcode = metier_shortcode,
    space = "GSA",
    time = "Q",
    technical = "lev4",
    target_SPECIES = target_species,
    output_dir = out_dir
)

?data_investigation

# =========================================================
# 8) run sampling_optimization()
# =========================================================
# The results will be saved in the folder DataInvestigation, automatically created in your out_dir directory.


# The path of the past situation csv file, produced by the data_investigation function need to be reported below.

past_situations <- read.csv("DataInvestigation/Past_situation.csv", sep = ";")
head(past_situations)

# the minimum number of trips considered acceptable for each stratum (e.g. quarter-gear) needs to be indicated (e.g. 3, corresponding to one trip per month)

table_min_accepted_sample_size <- read.csv("C:\\_LAVORO CORRENTE\\05_1_min_accepted_sample_size.csv", sep = ";")
head(table_min_accepted_sample_size)

res <- sampling_optimization(
    CS = costCS_from_CS,
    CL = costCL_from_CL,
    metier_shortcode = metier_shortcode,
    table_min_accepted_sample_size=table_min_accepted_sample_size,
    past_situations = past_situations,
    space = "GSA",
    time = "Q",
    technical = "lev4",
    target_SPECIES = target_species,
    nIter = 20,
    output_dir = out_dir
)

?sampling_optimization
