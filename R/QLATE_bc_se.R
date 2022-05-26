
#' Implements Quantile-LATE from Dong, Lee, and Gou (2021)
#'
#' @description
#' Implements the QLATE estimate proposed in Dong, Lee, and Gou and uses a
#' bias-corrected standard error
#'
#' @param df Data.frame.
#' @param qlist Vector of quantiles to estimate QLATE at. Must be evenly spaced
#' @param yname Character. Variable name for outcome variable
#' @param tname Character. Variable name for continuous treatment
#' @param rname Character. Variable name for running variable.
#' @param c Numeric. Cutoff for running variable. Default is 0.
#' @param verbose Logical. If TRUE, then print results to the console
#'
#' @export
QLATE_bc_se <- function(df, qlist, yname, tname, rname, c = 0, verbose = TRUE) {

	df = create_variables(df, yname, tname, rname, c)

	# Fix `no visible binding for global variables` problem
	qnum <-	qgrid <- n <- n_plus <- n_minus <- sd_R <- sd_T <- Z <- R <- NULL
	Tubp <-	Tubm <-	Tlbp <-	Tlbm <- NULL
	kern_Rf <- biasRPi <- biasTPi <- VpimA <- VpiqA <- NULL
	SEtau <- SEtau_bc <- SEpi <- pi_bc <- SEpi_bc <- NULL

	# Create empty vectors
	qte <- c() # q^{+}(u) - q^{-}(u)
	qplus <- c() # q^{+}(u)
	qminus <- c() # q^{-}(u)
	tau_u <- c() # tau(u) = [m^{+}(u) - m^{-}(u)] / [q^{+}(u) - q^{-}(u)]
	tau_u_bc <- c() # bias-corrected version of tau(u)
	Vtau <- c() # Variance of tau(u)
	ifT_Rplus <- c()
	ifT_Rminus <- c()
	Vtau <- c() # Variance of tau(u)
	dmdt_minus <- c()
	dmdt_plus <- c()
	L1_minus <- c()
	L1_plus <- c()
	biasRTau_u <- c() # Estimated bias for Reduced form tau(u)
	biasTTau_u <- c() # Estimated bias for T tau(u)
	opt_h_Rtau <- c() # Optimal bandwidth for Reduced form  tau(u)
	opt_h_Ttau <- c() # Optimal bandwidth for T tau(u)
	biasPi_u1 <- c() # Estimated bias terms in biasRPi (Pi)
	biasRPi_u2 <- c() # Estimated bias terms in biasRPi (Pi)
	biasTPi_u <- c() # Estimated bias terms in biasRPi (Pi)

	# Create empty values
	r <- 0
	sigma2p <- 0
	sigma2m <- 0

	summ <- calculate_summary_values(df, qlist)
	list2env(summ, environment())

	# ---- QLATE_bc_se Constants -------------------------------------------------

	# AMSE optimal bandwidth for Pi
	# Note that c = 0 now that I recentered above!
	optimal_bw_est = AMSE_Optimal_BW(df, qlist, yname, tname, rname, c=0, verbose = F)
	h_RPi = optimal_bw_est$h_Rpi
	h_TPi = optimal_bw_est$h_Tpi

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



	# ---- Calculate trimming parameter ------------------------------------------

	ret <- estimate_trim(df, qlist, h_R, h_R0, h_R1)
	df <- ret$df
	trim <- ret$trim


  # ---- Loop through qlist to estimate $tau(u)$ ------------------------------

  for(i in 1:length(qlist)) {

  	# ---- Estimate tau(u) ----------------------------------------------------

  	q <- qlist[i]

  	tau_est <- estimate_tau_u(df, i, q, qlist, qte, qminus, qplus, trim, tau_u,
  														h_T, h_T1)
  	list2env(tau_est, environment())


  	# ---- Estimate Variance of tau(u) ----------------------------------------
  	# Estimate sigma^{2+}, sigma^{2-} in V_tau

  	Vtau_est <- estimate_Vtau_u(df, i, qte, qgrid, Vtau,
  															ifT_Rplus, ifT_Rminus, fRplus, fRminus, fR,
  															h_fplus, h_fminus, sd_R, sd_T,
  															sigma2p, sigma2m)
  	list2env(Vtau_est, environment())



  	# ---- Bias-corrected tau(u) ----------------------------------------------
		# Estimate Bias of the estimated tau_u and optimal bandwidth for estimating tau_u

  	tau_bc_est <- estimate_tau_u_bc(df, i, q, qte, tau_u, tau_u_bc,
  																	n, h_R, h_T, Vtau,
  																	ifT_Rplus, ifT_Rminus, L1_plus, L1_minus,
  																	biasPi_u1, biasRPi_u2, biasTPi_u,
  																	opt_h_Rtau, opt_h_Ttau)
  	list2env(tau_bc_est, environment())

  } # End for loop


  # ---- Standard Error for $tau(u)$ ------------------------------------------

  tau_u_bc_se_est = estimate_tau_u_bc_se(Vtau, n, h_R, h_T, rho_1)
  list2env(tau_u_bc_se_est, environment())


  # ---- Estimate Pi and VPi ---------------------------------------------------

  pi_est = estimate_pi(qlist, qte, tau_u, biasPi_u1, biasRPi_u2, biasTPi_u,
  										 n, qgrid, qnum, qplus, qminus, h_R, h_T, fR,
  										 Tubp, Tlbp, Tubm, Tlbm, sigma2p, sigma2m,
  										 L1_plus, L1_minus, ifT_Rplus, ifT_Rminus)
  list2env(pi_est, environment())


  # ---- Bias-corrected Pi and VPi ---------------------------------------------

	Vpi_bc_est = estimate_Vpi_bc(pi, n, h_R, h_T, biasRPi, biasTPi,
														 rho_1, VpimA, VpiqA)
  list2env(Vpi_bc_est, environment())


  # Print results
  if(verbose) {
		cli::cli({
			cli::cli_text("{.code qlist}: {qlist}")
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

  return(list(
  	qlist = qlist,
  	tau_u = tau_u,
  	SEtau = SEtau,
  	tau_u_bc = tau_u_bc,
  	SEtau_bc = SEtau_bc,
  	pi = pi,
  	SEpi = SEpi,
  	pi_bc = pi_bc,
  	SEpi_bc = SEpi_bc
  ))
}
