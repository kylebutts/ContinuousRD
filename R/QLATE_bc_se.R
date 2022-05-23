
#' Implements Quantile-LATE from Dong, Lee, and Gou (2021)
#'
#' @description
#' Implements the QLATE estimate proposed in Dong, Lee, and Gou and uses a
#' bias-corrected standard error
#'
#'
#' @param qlist Vector of quantiles to estimate QLATE at
#' @param df Data.frame. For now, need the following columns:
#'   Y (outcome variable), T (continuous treatment),
#'   R (running variable centered at 0), Z (cutoff indicator),
#'   ZR (cutoff indicator x running variable),
#'   R2 (running variable^2),
#'   and ZR2 (cutoff indicator x running variable^2)
#'
QLATE_bc_se <- function(df, qlist) {

	data.table::setDT(df) # data.table
	df[, Z := as.numeric(df$Z)] # if F/T make 0/1

  # Calculate summary values
  qnum <- length(qlist)
	grid <- qlist[2] - qlist[1]
  n <- nrow(df)
	n_plus <- nrow(df[R >= 0])
	n_minus <- nrow(df[R < 0])
	sd_R <- df[, stats::sd(R)]
	sd_T <- df[, stats::sd(T)]

	# TODO: Grab optimal bandwidths from AMSE_Optimal_BW programatically
	# AMSE optimal bandwidth for Pi
	h_RPi = 2.2288
	h_TPi = 1.0049

	# constant for the bandwidth for bias and variance estimation
	c = 4.5
	# rho_0 is the ratio of the main bandwidth to the preliminary bandwidth used to estimate the trimming thresholds
	rho_0 = 4 / 3
	# rho_1 is the ratio of the main bandwidth to the bandwidth used to estimate biases and variances
	rho_1 = h_RPi / (c * n^(-1/8) * sd_R)

	# h_R and h_T are the bandwidths for the main estimation
	h_R = h_RPi
	h_T = h_TPi

	# h_R0 is the preliminary bandwidth for generating the triming parameter
	h_R0 = h_R / rho_0

	# h_R1 and h_T1 the bandwidths for estimating the biases and variances (h/b)
	h_R1 = h_R / rho_1
	h_T1 = h_R1 * sd_T / sd_R

  df[, kern_R := ifelse(abs(R) > h_R, 0, 0.5)]
  df[, kern_R0 := ifelse(abs(R) > h_R0, 0, 0.5)]
  df[, kern_R1 := ifelse(abs(R) > h_R1, 0, 0.5)]

  trim <- rq(
  	T ~ Z + R + ZR,
  	tau = qlist,
  	data = df,
  	weights = df$kern_R0
  ) |>
  	summary() |>
  	lapply(\(x) {
  		return(1.96 * x$coefficients["Z", "Std. Error"])
  	}) |>
  	unlist()

  # use uniform kernel
  # Silverman rule-of-thumb bandwidth for a unifrom kernel for density estimation
  h_Rf <- 1.843 * n^(-0.2) * sd_R
  df[, kern_Rf := ifelse(abs(R) > h_Rf, 0, 0.5)]
  fR <- df[, mean(kern_Rf)] / h_Rf

  # Estimate the density of R right above r0
  # Silverman rule-of-thumb bandwidth for a unifrom kernel
  h_fplus <- 0.7344 * n_plus^(-1 / 6)
  df[, kern_Rf := ifelse(abs(R) > h_fplus * sd_R, 0, 0.5)]
  fRplus <- df[Z == 1, mean(kern_Rf)] / (h_fplus * sd_R)

  # Estimate the density of R right below r0
  # Silverman rule-of-thumb bandwidth for a unifrom kernel
  h_fminus <- 0.7344 * n_minus^(-1 / 6)
  df[, kern_Rf := ifelse(abs(R) > h_fminus * sd_R, 0, 0.5)]
  fRminus <- df[Z == 0, mean(kern_Rf)] / (h_fminus * sd_R)




  # Create empty values
  r <- 0
  sigma2p <- 0
  sigma2m <- 0

  # Create empty vectors
  qte <- c()        # q^{+}(u) - q^{-}(u)
  qplus <- c()      # q^{+}(u)
  qminus <- c()     # q^{-}(u)
  tau_u <- c()      # \tau(u) = [m^{+}(u) - m^{-}(u)] / [q^{+}(u) - q^{-}(u)]
  tau_u_bc <- c()   # bias-corrected version of \tau(u)
  Vtau <- c() # Variance of \tau(u)
  biasPi_u1 <- c() # Estimated bias terms
  biasRPi_u2 <- c() # Estimated bias terms
  biasTPi_u2 <- c() # Estimated bias terms
  dmdt_minus <- c()
  dmdt_plus <- c()
  L1_minus <- c()
  L1_plus <- c()

  for(i in 1:length(qlist)) {

  	# ---- Estimate \tau(u) ----------------------------------------------------
  	q <- qlist[i]
  	model <- rq(T ~ Z + R + ZR, q, df, weights = df$kern_R)

  	# Get estimates of q^{+}(u) and q^{-}(u)
  	qte[i] <- stats::coef(model)["Z"]
  	qminus[i] <- stats::coef(model)["(Intercept)"]
  	qplus[i] <- qminus[i] + qte[i]

		# Purge T of q^{+} and q^{-}
		df[, T_c := ifelse(R >= 0, T - qplus[i], T - qminus[i]) ]
		df[, `:=`(
			T_c2 = T_c^2,
			ZT = T_c * Z,
			ZT2 = T_c^2 * Z,
			RT_c = T_c * R,
			ZRT_c = T_c * R * Z
		)]

		# In the following, `qte` is set to be missing at the trimmed quantiles so they will drop in the integration of \tau(u) `tau_u` to obtain \pi. In addition, quantile index `qlist` is set to be zero, so in the double integration for variance calculation, they will drop
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
		df[, kern_prod1 := kern_R1 * kern_T1]

		# Estimate uncorrected \tau(u)
		reg <- fixest::feols(Y ~ R + T_c + Z + ZR + ZT, df, weights = ~kern_prod)
		tau_u[i] <- stats::coef(reg)["Z"] / qte[i]

		# ---- Estimate Variance of \tau(u) ----------------------------------------
		# Estimate sigma^{2+}, sigma^{2-} in V_\tau
		reg <- fixest::feols(
			Y ~ R + T_c + Z + ZR + ZT,
			df,
			weights = ~kern_prod1,
			warn = F
		)
		df[, Y2 := (Y - stats::coef(reg)["(Intercept)"] - coef(reg)["Z"])^2]
		df[R < 0, Y2 := (Y - stats::coef(reg)["(Intercept)"])^2]
		reg <- fixest::feols(
			Y2 ~ R + T_c + Z + ZR + ZT,
			df,
			weights = ~kern_prod1,
			warn = F
		)
		sigma2plus <- stats::coef(reg)["(Intercept)"] + coef(reg)["Z"]
		sigma2minus <- stats::coef(reg)["(Intercept)"]

		# Integration in Vpim+ and Vpim-
		if (!is.na(qte[i])) {
			sigma2p <- sigma2p + sigma2plus * grid
			sigma2m <- sigma2m + sigma2minus * grid
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


		# ---- Bias-corrected \tau(u) ------------------------------------------------
		# Estimate Bias of the estimated tau_u and optimal bandwidth for estimating tau_u
		model <- rq(T ~ Z + R + R2 + ZR + ZR2, q, df, weights = df$kern_R1)
		biasT_minus <- 2 * stats::coef(model)["R2"] * (-1 / 12)
		biasT_plus <- 2 * (stats::coef(model)["R2"] + stats::coef(model)["ZR2"]) * (-1 / 12)

		lm <- fixest::feols(
			Y ~ R + R2 + T_c + T_c2 + RT_c + Z + ZR + ZT + ZR2 + ZT2 + ZRT_c,
			df,
			weights = ~kern_prod1,
			warn = F
		)

		# Auxiliary terms for bias
		BR2 <- 2 * stats::coef(lm)["ZR2"] * (-1 / 12)
		BT2 <- 2 * stats::coef(lm)["ZT2"] * (1 / 6)
		dmdt_minus[i] <- stats::coef(lm)["T_c"]
		dmdt_plus[i] <- stats::coef(lm)["T_c"] + stats::coef(lm)["ZT"]

		L1_minus[i] = dmdt_minus[i] * ifT_Rminus[i]
		L1_plus[i] = dmdt_plus[i] * ifT_Rplus[i]

		biasRTau_u <- (BR2 + biasT_plus * (dmdt_plus[i] - tau_u[i]) - biasT_minus * (dmdt_minus[i] - tau_u[i])) / qte[i]
		biasTTau_u <- BT2 / qte[i]

		tau_u_bc[i] = tau_u[i] - h_R^2 * biasRTau_u - h_T^2 * biasTTau_u

		# Estimated bias terms in biasRPi
		biasPi_u1[i] <- (biasT_plus - biasT_minus) / qte[i]
		biasRPi_u2[i] <- biasRTau_u + (biasT_plus - biasT_minus) * tau_u[i] / qte[i]
		biasTPi_u2[i] <- biasTTau_u

  } # End for loop

  # Standard Error for Tau
  SEtau = sqrt(Vtau * 1 / (n * h_R * h_T) )

  if(rho_1 >= 1) {
  	contCtau = 37.5 * (rho_1/3 - 0.25)
  } else {
  	contCtau = 3.125 * rho_1^3
  }

  # Bias-corrected Standard error for Tau
  SEtau_bc = sqrt(1 + rho_1^6 * 9.765625 + rho_1 * contCtau) * SEtau


  # Estimate overall treatment effect pi
  w_sum <- sum(abs(qte), na.rm = T)
  pi <- sum(tau_u * abs(qte) / w_sum, na.rm = T)

  # Estimate bias of pi
  biasRPi_u <- biasRPi_u2 - pi * biasPi_u1
  biasRPi <- sum(biasRPi_u * abs(qte) / w_sum, na.rm = T)
  biasTPi_u <- biasTPi_u2
	biasTPi <- sum(biasTPi_u * abs(qte) / w_sum, na.rm = T)

	pi_bc <- pi - h_R^2 * biasRPi - h_T^2 * biasTPi

	# Estimate Variance of pi: Vpi = Vpim + Vpiq
	# Compute the adjustment term for Vpim
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

	# Compute the adjustment term for Vpiq
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
	VpimA <- 4 * sigma2A / (fR * w_sum^2 * grid^2)

	# Estimate parameters in Vpiq
	Lplus <- (L1_plus - pi * ifT_Rplus) / (w_sum * grid)
	Lminus <- (L1_minus - pi * ifT_Rminus) / (w_sum * grid)

	Ing <- 0
	for (j in 1:qnum) {
		for (k in 1:qnum) {
			uj <- qlist[j]
			vk <- qlist[k]
			Lkj <- Lplus[j] * Lplus[k] + Lminus[j] * Lminus[k]

			if (k <= j) {
				Ing <- Ing + vk * (1 - uj) * Lkj * grid^2
			}
			if (k > j) {
				Ing <- Ing + uj * (1 - vk) * Lkj * grid^2
			}
		}
	}

	Vpiq = 4 * Ing/fR
	Vpi = VpimA + Vpiq
	SEpi = sqrt(VpiA / (n * h_R))

	if(rho_1 >= 1) {
		contCpi = 2.5 - 1.875 / rho_1
	} else {
		contCpi = 3.125 * rho_1 - 2.5 * rho_1^3
	}

	Vpi_bc = (1 + 1.640625 * rho_1^5 + contCpi * rho_1^2)*VpimA + Vpiq
	SEpi_bc = sqrt(Vpi_bc / (n * h_R))

  # Print results
	cli::cli({
		cli::cli_text("{.code tau_u}: {round(tau_u, 3)}")
		cli::cli_text("{.code SEtau}: {round(SEtau, 3)}")
		cli::cli_text("{.code tau_u_bc}: {round(tau_u_bc, 3)}")
		cli::cli_text("{.code SEtau_bc}: {round(SEtau_bc, 3)}")
		cli::cli_text("{.code pi}: {round(pi, 3)}")
		cli::cli_text("{.code SEpi}: {round(SEpi, 3)}")
		cli::cli_text("{.code pi_bc}: {round(pi_bc, 3)}")
		cli::cli_text("{.code SEpi_bc}: {round(SEpi_bc, 3)}")
	})
}
