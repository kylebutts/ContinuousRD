## sim-data.R ------------------------------------------------------------------
## Kyle Butts, CU Boulder Economics
##
## Simulate dataset for testing.

library(data.table)

# ---- Generate Dataset --------------------------------------------------------

set.seed(5)
df <- data.table::CJ(id = 1:5000)
df[, let(u = runif(.N, 0, 1), x = runif(.N, -1, 1))]
df[, Z := as.numeric(x > 0)]
df[, let(
	# Strictly increasing in u
	# Treatment effect up to about 30th percentile
	T = (10/3*u) * (u <= 0.3) * (Z == 0)
	+ (0.75*u + 0.8) * (u > 0.3) * (Z == 0)
	+ (0.75*u + 0.8) * (Z == 1)
)]
df[, let(
	R = x,
	R2 = x^2,
	ZR = Z * x,
	ZR2 = Z * x^2,
	# Marginal effect of 1
	Y = T + rnorm(.N, 0, 0.005)
)]

# Looks similar to Figure 2
# ggplot(df) +
# 	geom_point(
# 		aes(x = u, y = T, color = u),
# 		alpha = 0.4
# 	)

# "treatment Effects" test
# plot(df[u >  0.3 & abs(x) < 0.1]$x, df[u >  0.3 & abs(x) < 0.1]$Y)
# plot(df[u <= 0.3 & abs(x) < 0.1]$x, df[u <= 0.3 & abs(x) < 0.1]$Y)


# True treatment effect = 1 if u < 0.3
# summ <- df[
# 	u <= 0.3 & abs(x) < 0.1,
# 	.(T = mean(T), Y = mean(Y)),
# 	by = Z
# ]
# (summ[Z == 1, Y] - summ[Z == 0, Y]) / (summ[Z == 1, T] - summ[Z == 0, T])

# ---- Export dataset ----------------------------------------------------------
fwrite(df, here::here("data/sim_data.csv"))
usethis::use_data(df, overwrite = TRUE)
