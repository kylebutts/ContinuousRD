## -----------------------------------------------------------------------------
## Utility functions for `QLATE_bc_se` and `AMSE_Optimal_BW` to hopefully
## increase readability of main functions and to not repeat myself
## -----------------------------------------------------------------------------


create_variables = function(df, yname, tname, rname, c) {

	# Use data.table
	data.table::setDT(df) # data.table


	# Center running variable
	df$Y <- df[[yname]]
	df$T <- df[[tname]]
	df$R <- df[[rname]] - c
	df$Z <- as.numeric(df$R > 0)
	df$ZR <- df$Z * df$R
	df$R2 <- df$R^2
	df$ZR2 <- df$Z * df$R^2

	return(df)
}


calculate_summary_values = function(df, qlist) {
	# Fix `no visible binding for global variables` problem
	R <- NULL

	# Calculate summary values
	qnum <- length(qlist)
	qgrid <- qlist[2] - qlist[1]
	n <- nrow(df)
	n_plus <- nrow(df[R >= 0, ])
	n_minus <- nrow(df[R < 0, ])
	sd_R <- df[, stats::sd(R)]
	sd_T <- df[, stats::sd(T)]

	# For the adjustment term in VpimA, the upper and lower bounds of the support of T1 and T0
	Tubp <- df[R >= 0, max(T)]
	Tubm <- df[R < 0, max(T)]
	Tlbp <- df[R >= 0, min(T)]
	Tlbm <- df[R < 0, min(T)]

	return(list(
		qnum = qnum, qgrid = qgrid, n = n, n_plus = n_plus, n_minus = n_minus,
		sd_R = sd_R, sd_T = sd_T, Tubp = Tubp, Tubm = Tubm, Tlbp = Tlbp,
		Tlbm = Tlbm
	))
}


# Generate trimming threshold with uniform kernel
estimate_trim <- function(df, qlist, h_R, h_R0, h_R1) {
	# Fix `no visible binding for global variables` problem
	R <- kern_R <- kern_R0 <- kern_R1 <- rq <- NULL

	df[, kern_R := ifelse(abs(R) > h_R, 0, 0.5)]
	df[, kern_R0 := ifelse(abs(R) > h_R0, 0, 0.5)]
	df[, kern_R1 := ifelse(abs(R) > h_R1, 0, 0.5)]

	trim <- quantreg::rq(
		T ~ Z + R + ZR,
		tau = qlist,
		data = df[kern_R0 > 0, ],
		weights = df[kern_R0 > 0, ]$kern_R0
	) |>
		summary() |>
		lapply(\(x) {
			return(1.96 * x$coefficients["Z", "Std. Error"])
		}) |>
		unlist()

	return(list(
		df = df, trim = trim
	))
}

# Estimate $tau(u)$
estimate_tau_u = function(df, i, q, qlist, qte, qminus, qplus, trim, tau_u,
													h_T, h_T1) {

	# Fix `no visible binding for global variables` problem
	R <- Z <- T_c <- kern_T <- kern_prod <- kern_R <- kern_T1 <- kern_prod1 <- NULL

	model <- quantreg::rq(T ~ Z + R + ZR, q, df, weights = df$kern_R)

	# Get estimates of q^{+}(u) and q^{-}(u)
	qte[i] <- stats::coef(model)["Z"]
	qminus[i] <- stats::coef(model)["(Intercept)"]
	qplus[i] <- qminus[i] + qte[i]

	# Purge T of q^{+} and q^{-}
	df[
		, T_c := ifelse(R >= 0, T - qplus[i], T - qminus[i])
	]

	df[, `:=`(
		T_c2 = T_c^2,
		ZT = T_c * Z,
		ZT2 = T_c^2 * Z,
		RT_c = T_c * R,
		ZRT_c = T_c * R * Z
	)]

	# In the following, `qte` is set to be missing at the trimmed quantiles so they will drop in the integration of tau(u) `tau_u` to obtain Pi. In addition, quantile index `qlist` is set to be zero, so in the double integration for variance calculation, they will drop
	if (abs(qte[i]) < trim[i]) {
		qte[i] <- NA
		qlist[i] <- 0
		qplus[i] <- NA
		qminus[i] <- NA
	}

	# use uniform kernel
	df[, kern_T := ifelse(abs(T_c) > h_T, 0, 0.5)]
	df[, kern_prod := kern_R * kern_T]
	df[, kern_T1 := ifelse(abs(T_c) > h_T1, 0, 0.5)]
	# TODO: Check this. In AMSE_Optimal_BW, uses kern_R
	# In QLATE_bc_se, uses kern_R1
	df[, kern_prod1 := kern_R * kern_T1]

	# Estimate non bias-corrected tau(u)
	reg <- fixest::feols(Y ~ R + T_c + Z + ZR + ZT, df[kern_prod > 0, ], weights = ~kern_prod)
	tau_u[i] <- stats::coef(reg)["Z"] / qte[i]

	return(list(
		df = df, qte = qte, qminus = qminus, qplus = qplus, trim = trim, tau_u = tau_u
	))
}

# Estimate $sigma^{2+}$, $sigma^{2-}$ in V_{tau}
estimate_Vtau_u = function(df, i, qte, qgrid, Vtau,
													 ifT_Rplus, ifT_Rminus, fRplus, fRminus, fR,
													 h_fplus, h_fminus, sd_R, sd_T,
													 sigma2p, sigma2m) {

	# Fix `no visible binding for global variables` problem
	Y <- Y2 <- R <- Z <- T_c <- kern_Rf <- kern_prod1 <- NULL

	reg <- fixest::feols(
		Y ~ R + T_c + Z + ZR + ZT,
		df[kern_prod1 > 0, ],
		weights = ~kern_prod1,
		warn = F
	)
	df[, Y2 := (Y - stats::coef(reg)["(Intercept)"] - stats::coef(reg)["Z"])^2]
	df[R < 0, Y2 := (Y - stats::coef(reg)["(Intercept)"])^2]

	reg <- fixest::feols(
		Y2 ~ R + T_c + Z + ZR + ZT,
		df[kern_prod1 > 0, ],
		weights = ~kern_prod1,
		warn = F
	)
	sigma2plus <- stats::coef(reg)["(Intercept)"] + stats::coef(reg)["Z"]
	sigma2minus <- stats::coef(reg)["(Intercept)"]

	# Integration in VpimA+ and VpimA-
	if (!is.na(qte[i])) {
		sigma2p <- sigma2p + sigma2plus * qgrid
		sigma2m <- sigma2m + sigma2minus * qgrid
	}

	# Parameters in V_tau
	# Estimate the joint density of R and T right above r0
	df[, kern_Rf := ifelse(
		abs(R) > (h_fplus * sd_R) | abs(T_c) > (h_fplus * sd_T),
		0, 0.25
	)]
	fTRplus <- df[Z == 1, mean(kern_Rf)] / (h_fplus^2 * sd_R * sd_T)

	# Estimate the reciprocal of the conditional density of T given R right above r0.
	# That is, ifT_Rplus = 1/f_{TR}^+
	ifT_Rplus[i] <- fRplus / fTRplus

	# Estimate the joint density of R and T right below r0
	df[, kern_Rf := ifelse(
		abs(R) > (h_fminus * sd_R) | abs(T_c) > (h_fminus * sd_T),
		0, 0.25
	)]
	fTRminus <- df[Z == 0, mean(kern_Rf)] / (h_fminus^2 * sd_R * sd_T)

	# Estimate the reciprocal of the conditional density of T given R right below r0.
	# That is, ifT_R_minus = 1/f_{TR}^-
	ifT_Rminus[i] <- fRminus / fTRminus

	Vtau[i] <- 2 * (sigma2plus * ifT_Rplus[i] + sigma2minus * ifT_Rminus[i]) /
		(fR * qte[i]^2)

	return(list(
		df = df, ifT_Rplus = ifT_Rplus, ifT_Rminus = ifT_Rminus,
		sigma2p = sigma2p, sigma2m = sigma2m, Vtau = Vtau
	))
}



# Estimate bias-corrected $tau(u)$ and optimal bandwidth (for AMSE_Optimal_BW)
estimate_tau_u_bc = function(df, i, q, qte, tau_u, tau_u_bc,
														 n, h_R, h_T, Vtau,
														 ifT_Rplus, ifT_Rminus, L1_plus, L1_minus,
														 biasPi_u1, biasRPi_u2, biasTPi_u,
														 opt_h_Rtau, opt_h_Ttau) {

	# Fix `no visible binding for global variables` problem
	kern_R1 <- kern_prod1 <- NULL

	model <- quantreg::rq(T ~ Z + R + R2 + ZR + ZR2, q,
							df[kern_R1 > 0, ], weights = df[kern_R1 > 0, ]$kern_R1
					 )
	biasT_minus <- 2 * stats::coef(model)["R2"] * (-1 / 12)
	biasT_plus <- 2 * (stats::coef(model)["R2"] + stats::coef(model)["ZR2"]) * (-1 / 12)

	lm <- fixest::feols(
		Y ~ R + R2 + T_c + T_c2 + RT_c + Z + ZR + ZT + ZR2 + ZT2 + ZRT_c,
		df[kern_prod1 > 0,  ],
		weights = ~kern_prod1,
		warn = F
	)

	# Auxiliary terms for bias
	BR2 <- 2 * stats::coef(lm)["ZR2"] * (-1 / 12)
	BT2 <- 2 * stats::coef(lm)["ZT2"] * (1 / 6)
	dmdt_minus <- stats::coef(lm)["T_c"]
	dmdt_plus <- stats::coef(lm)["T_c"] + stats::coef(lm)["ZT"]

	L1_minus[i] <- dmdt_minus * ifT_Rminus[i]
	L1_plus[i] <- dmdt_plus * ifT_Rplus[i]

	biasRTau_u <- (BR2 + biasT_plus * (dmdt_plus - tau_u[i]) - biasT_minus * (dmdt_minus - tau_u[i])) / qte[i]
	biasTTau_u <- BT2 / qte[i]

	# Bias-corrected $tau(u)$
	tau_u_bc[i] = tau_u[i] - h_R^2 * biasRTau_u - h_T^2 * biasTTau_u

	# Estimated bias terms in biasRPi
	biasPi_u1[i] <- (biasT_plus - biasT_minus) / qte[i]
	biasRPi_u2[i] <- biasRTau_u + (biasT_plus - biasT_minus) * tau_u[i] / qte[i]
	biasTPi_u[i] <- biasTTau_u

	# NB: These terms are only used by AMSE_Optimal_BW, but easier than writing seperate functions at minimal cost
	# Optimal bandwidth terms
	opt_h_Rtau[i] <- (Vtau[i] / (8 * n))^(1 / 6) * (abs(biasTTau_u) / abs(biasRTau_u)^5)^(1 / 12)
	opt_h_Ttau[i] <- (Vtau[i] / (8 * n))^(1 / 6) * (abs(biasRTau_u) / abs(biasTTau_u)^5)^(1 / 12)

	return(list(
		df = df, L1_plus = L1_plus, L1_minus = L1_minus,
		tau_u = tau_u, tau_u_bc = tau_u_bc,
		biasPi_u1 = biasPi_u1, biasRPi_u2 = biasRPi_u2, biasTPi_u = biasTPi_u,
		opt_h_Rtau = opt_h_Rtau, opt_h_Ttau = opt_h_Ttau
	))
}


estimate_tau_u_bc_se = function(Vtau, n, h_R, h_T, rho_1) {
	SEtau = sqrt(Vtau * 1 / (n * h_R * h_T) )

	if(rho_1 >= 1) {
		contCtau = 37.5 * (rho_1/3 - 0.25)
	} else {
		contCtau = 3.125 * rho_1^3
	}

	# Bias-corrected Standard error for Tau
	SEtau_bc = sqrt(1 + rho_1^6 * 9.765625 + rho_1 * contCtau) * SEtau

	return(list(
		SEtau = SEtau, SEtau_bc = SEtau_bc
	))
}


# Estimate overall treatment effect Pi, bias-corrected Pi, and V_{Pi}
estimate_pi <- function(qlist, qte, tau_u, biasPi_u1, biasRPi_u2, biasTPi_u,
												n, qgrid, qnum, qplus, qminus, h_R, h_T, fR,
												Tubp, Tlbp, Tubm, Tlbm, sigma2p, sigma2m,
												L1_plus, L1_minus, ifT_Rplus, ifT_Rminus) {
	# Estimate overall treatment effect pi
	w_sum <- sum(abs(qte), na.rm = T)
	pi <- sum(tau_u * abs(qte) / w_sum, na.rm = T)

	# Estimate bias of pi
	biasRPi_u <- biasRPi_u2 - pi * biasPi_u1
	biasRPi <- sum(biasRPi_u * abs(qte) / w_sum, na.rm = T)
	biasTPi <- sum(biasTPi_u * abs(qte) / w_sum, na.rm = T)

	pi_bc <- pi - h_R^2 * biasRPi - h_T^2 * biasTPi

	# Estimate Variance of pi: Vpi = VpimA + VpiqA
	# Compute the adjustment term for VpimA
	qub <- max(qplus, na.rm = T)
	qlb <- min(qplus, na.rm = T)
	tsupp <- abs(qub - qlb) / h_T

	u1 <- (Tubp - qub) / h_T
	G1 <- max(min((u1 + 1) / 2, 1), 0)
	u2 <- (Tlbp - qub) / h_T
	G2 <- max(min((u2 + 1) / 2, 1), 0)

	if (tsupp >= 2) {
		Ap <- G1 - G2
	} else {
		u2 <- (qub - Tlbp) / h_T
		U <- min(u2, 1)
		u1 <- (qub - Tubp) / h_T
		L <- u1
		if (u1 < tsupp - 1) L <- tsupp - 1
		A1 <- U^2 / 8 + (1 - tsupp) * U / 4 - L^2 / 8 - (1 - tsupp) * L / 4

		u2 <- (Tubp - qlb) / h_T
		U <- min(u2, 1)
		u1 <- (Tlbp - qlb) / h_T
		L <- u1
		if (u1 < tsupp - 1) L <- tsupp - 1
		A2 <- U^2 / 8 + (1 - tsupp) * U / 4 - L^2 / 8 - (1 - tsupp) * L / 4

		Ap <- G1 - G2 - A1 - A2
	}

	# Compute the adjustment term for VpiqA
	qub <- max(qminus, na.rm = T)
	qlb <- min(qminus, na.rm = T)
	tsupp <- abs(qub - qlb) / h_T

	u1 <- (Tubm - qub) / h_T
	G1 <- max(min((u1 + 1) / 2, 1), 0)
	u2 <- (Tlbm - qub) / h_T
	G2 <- max(min((u2 + 1) / 2, 1), 0)

	if (tsupp >= 2) {
		Am <- G1 - G2
	} else {
		u2 <- (qub - Tlbm) / h_T
		U <- min(u2, 1)
		u1 <- (qub - Tubm) / h_T
		L <- u1
		if (u1 < tsupp - 1) L <- tsupp - 1
		A1 <- U^2 / 8 + (1 - tsupp) * U / 4 - L^2 / 8 - (1 - tsupp) * L / 4

		u2 <- (Tubm - qlb) / h_T
		U <- min(u2, 1)
		u1 <- (Tlbm - qlb) / h_T
		L <- u1
		if (u1 < tsupp - 1) L <- tsupp - 1
		A2 <- U^2 / 8 + (1 - tsupp) * U / 4 - L^2 / 8 - (1 - tsupp) * L / 4

		Am <- G1 - G2 - A1 - A2
	}

	sigma2A <- sigma2p * Ap + sigma2m * Am
	VpimA <- 4 * sigma2A / (fR * w_sum^2 * qgrid^2)

	# Estimate parameters in VpiqA
	Lplus <- (L1_plus - pi * ifT_Rplus) / (w_sum * qgrid)
	Lminus <- (L1_minus - pi * ifT_Rminus) / (w_sum * qgrid)

	Ing <- 0
	for (j in 1:qnum) {
		for (k in 1:qnum) {
			uj <- qlist[j]
			vk <- qlist[k]
			Lkj <- Lplus[j] * Lplus[k] + Lminus[j] * Lminus[k]

			if (k <= j) {
				Ing <- Ing + vk * (1 - uj) * Lkj * qgrid^2
			}
			if (k > j) {
				Ing <- Ing + uj * (1 - vk) * Lkj * qgrid^2
			}
		}
	}

	VpiqA = 4 * Ing/fR
	Vpi = VpimA + VpiqA
	SEpi = sqrt(Vpi / (n * h_R))

	return(list(
		pi = pi, pi_bc = pi_bc,
		VpimA = VpimA, VpiqA = VpiqA, Vpi = Vpi, SEpi = SEpi,
		biasRPi = biasRPi, biasTPi = biasTPi
	))
}


# Estimate optimal bandwidth for Pi
estimate_optimal_bw_pi <- function(opt_h_Rtau, opt_h_Ttau, Vpi, biasRPi,
																	 n, sd_T, sd_R, h_R1) {
	h_R_Rtau <- opt_h_Rtau
	h_R_Ttau <- opt_h_Ttau
	opt_h_RPiA <- (Vpi / (4 * n * biasRPi^2))^(1 / 5)
	opt_h_TPiA <- opt_h_RPiA * n^(-1 / 30) * sd_T / sd_R
	rho1_inf <- opt_h_RPiA / h_R1

	return(list(
		opt_h_RPiA = opt_h_RPiA, opt_h_TPiA = opt_h_TPiA, rho1_inf = rho1_inf
	))
}

# Estimate Bias-corrected Pi and V_Pi
estimate_Vpi_bc <- function(pi, n, h_R, h_T, biasRPi, biasTPi,
													 rho_1, VpimA, VpiqA) {
	if(rho_1 >= 1) {
		contCpi = 2.5 - 1.875 / rho_1
	} else {
		contCpi = 3.125 * rho_1 - 2.5 * rho_1^3
	}

	Vpi_bc = (1 + 1.640625 * rho_1^5 + contCpi * rho_1^2) * VpimA + VpiqA
	SEpi_bc = sqrt(Vpi_bc / (n * h_R))

	return(list(
		Vpi_bc = Vpi_bc, SEpi_bc = SEpi_bc
	))
}






