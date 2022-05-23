library(data.table)
library(ggplot2)
library(quantreg)

df <- fread(here::here("data/sim_data.csv"))
qlist <- seq(0.1, 0.5, by = 0.1)

# Average Mean Square Error optimal bandwidth
AMSE_Optimal_BW(df, qlist)

#
QLATE_bc_se(df, qlist)


