library(data.table)
library(ggplot2)
library(quantreg)

df <- fread(here::here("data-raw/sim_data.csv"))
qlist <- seq(0.1, 0.5, by = 0.1)

# Average Mean Square Error optimal bandwidth
AMSE_Optimal_BW(df, qlist, yname = "Y", tname = "T", rname = "R")

# Q-LATE with Bias-Correction
QLATE_bc_se(df, qlist, yname = "Y", tname = "T", rname = "R")


