#install.packages("remotes")

# remotes::install_github("https://github.com/COISPA/SDTool/tree/main", dependencies = TRUE, upgrade = "never")

# INSTALLA PACCHETTO ZIP


# 1) Load package
library(SDTool)
setwd("C:/_LAVORO CORRENTE")
# 2) Load dummy datasets shipped with the package
data("cs_gsa99_dummy", package = "SDTool")
data("cl_gsa99_dummy", package = "SDTool")

out_dir="C:\\_LAVORO CORRENTE"

# 3) Run functions
# --- A) CS -> CL (replace with your real function name)
costCS_from_CS <- convert_cs_to_costCS(cs_gsa99_dummy,write_csv = TRUE,out_dir=out_dir)

# --- B) CL -> COST clData (replace with your real function name)
costCL_from_CL <- convert_cl_to_costCL(cl_gsa99_dummy,write_csv = TRUE,out_dir=out_dir)




# =========================================================
# TEST data_investigation()
# =========================================================

# IMPORTANT:
# the function needs a metier shortcode table.

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

# Run the function
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


# =========================================================
# TEST sampling_optimization()
# =========================================================

past_situations <- read.csv("DataInvestigation/Past_situation.csv", sep = ";")
table_min_accepted_sample_size <- read.csv("C:\\_LAVORO CORRENTE\\05_1_min_accepted_sample_size.csv", sep = ";")


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
    nIter = 5,      # basso per test veloce
    output_dir = out_dir
)
