## sim-data.R ------------------------------------------------------------------
## Kyle Butts, CU Boulder Economics
##
## Simulate dataset for testing.

library(data.table)

# ---- Generate Dataset --------------------------------------------------------

set.seed(5)
sim_data <- data.table::CJ(id = 1:10000)
sim_data[, let(u = runif(.N, 0, 1), R = runif(.N, -1, 1))]
sim_data[, Z := as.numeric(R > 0)]
sim_data[, let(
	# Strictly increasing in u
	# Treatment effect up to about 30th percentile
	T = (10/3*u) * (u <= 0.3) * (Z == 0)
	+ (0.75*u + 0.8) * (u > 0.3) * (Z == 0)
	+ (0.75*u + 0.8) * (Z == 1)
)]
sim_data[, let(
	# Marginal effect of 1
	Y = T + rnorm(.N, 0, 0.005)
)]

# Looks similar to Figure 2
# ggplot(sim_data) +
# 	geom_point(
# 		aes(R = u, y = T, color = u),
# 		alpha = 0.4
# 	)

# "treatment Effects" test
# plot(sim_data[u >  0.3 & abs(R) < 0.1]$R, sim_data[u >  0.3 & abs(R) < 0.1]$Y)
# plot(sim_data[u <= 0.3 & abs(R) < 0.1]$R, sim_data[u <= 0.3 & abs(R) < 0.1]$Y)


# True treatment effect = 1 if u < 0.3
# summ <- sim_data[
# 	u <= 0.3 & abs(R) < 0.1,
# 	.(T = mean(T), Y = mean(Y)),
# 	by = Z
# ]
# (summ[Z == 1, Y] - summ[Z == 0, Y]) / (summ[Z == 1, T] - summ[Z == 0, T])

# ---- ERport dataset ----------------------------------------------------------

fwrite(sim_data, here::here("data-raw/sim_data.csv"))
usethis::use_data(sim_data, overwrite = TRUE)
