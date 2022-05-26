
#' AMSE-Optimal Bandwidth Selector for Q-LATE Estimator
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
AMSE_Optimal_BW = function(df, qlist, yname, tname, rname, c = 0, verbose = TRUE) {

	df = create_variables(df, yname, tname, rname, c)

	# Fix `no visible binding for global variables` problem
	qnum <-	qgrid <- n <-	n_plus <-	n_minus <-	sd_R <-	sd_T <- Z <- R <- NULL
	Tubp <-	Tubm <-	Tlbp <-	Tlbm <- NULL
	kern_Rf <- Vpi <- biasRPi <- opt_h_RPiA <- opt_h_TPiA <- rho1_inf <- NULL


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


	# ---- Calculate summary values ----------------------------------------------

	summ <- calculate_summary_values(df, qlist)
	list2env(summ, environment())

	# ---- AMSE_Optimal_BW Constants ---------------------------------------------

	# c is constant for the bandwidth for bias and variance estimation
	c = 4.5
	# a is the power in the bandwidth
	a = -1 / 8
	# rho_0 is the ratio of the main bandwidth to the preliminary bandwidth used to estimate the trimming thresholds
	rho_0 = 4 / 3
	# TODO: check this
	# rho_1 is the ratio of the main bandwidth to the bandwidth used to estimate biases and variances
	rho_1 = 1

	# h_std is the standardized bandwidth for the bias and variance estimation (b)
	h_std <- c * n^a
	h_R1 <- h_std * sd_R
	h_T1 <- h_std * sd_T

	# h_R and h_T the bandwidths for estimating tau_u and pi (h)
	h_R <- h_R1 * rho_1
	h_T <- h_T1 * rho_1

	# h_R0 is the bandwidth for generating the trimming parameter
	h_R0 <- h_R / rho_0

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
  for (i in 1:length(qlist)) {

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



    # ---- Bias-corrected tau(u) ------------------------------------------------
    # Estimate Bias of the estimated tau_u and optimal bandwidth for estimating tau_u

    tau_bc_est <- estimate_tau_u_bc(df, i, q, qte, tau_u, tau_u_bc,
    																n, h_R, h_T, Vtau,
    																ifT_Rplus, ifT_Rminus, L1_plus, L1_minus,
    																biasPi_u1, biasRPi_u2, biasTPi_u,
    																opt_h_Rtau, opt_h_Ttau)
    list2env(tau_bc_est, environment())

  } # End for loop


  # ---- Estimate Pi and VPi ---------------------------------------------------

  pi_est = estimate_pi(qlist, qte, tau_u, biasPi_u1, biasRPi_u2, biasTPi_u,
  										 n, qgrid, qnum, qplus, qminus, h_R, h_T, fR,
  										 Tubp, Tlbp, Tubm, Tlbm, sigma2p, sigma2m,
  										 L1_plus, L1_minus, ifT_Rplus, ifT_Rminus)
  list2env(pi_est, environment())

  # ---- Estimate optimal bandwidth for Pi -------------------------------------

  opt_bw_pi_est = estimate_optimal_bw_pi(opt_h_Rtau, opt_h_Ttau, Vpi, biasRPi,
  																			 n, sd_T, sd_R, h_R1)
  list2env(opt_bw_pi_est, environment())

  # ---- Display results -------------------------------------------------------

  if(verbose) {
  	# Display the initial standardized bandwidth used
  	cli::cli({
  		cli::cli_alert_info("Initial standardized bandwidth used")
  		cli::cli_text("{.code c}: {c}")
  		cli::cli_text("{.code h_std}: {h_std}")
  		cli::cli_text("{.code h_R}: {h_R}")
  		cli::cli_text("{.code h_T}: {h_T}")
  		cli::cli_text("{.code opt_h_Rtau}: {opt_h_Rtau}")
  		cli::cli_text("{.code opt_h_Ttau}: {opt_h_Ttau}")
  	})

  	# Display the adjusted optimal bandwidth for pi
  	cli::cli({
  		cli::cli_alert_info("Adjusted optimal bandwidth for pi")
  		cli::cli_text("{.code opt_h_RPiA}: {opt_h_RPiA}")
  		cli::cli_text("{.code opt_h_TPiA}: {opt_h_TPiA}")
  		cli::cli_text("{.code rho1_inf}: {rho1_inf}")
  	})
  }

  return(
  	list(
  		h_Rpi = opt_h_RPiA,
			h_Tpi = opt_h_TPiA,
  		rho1_inf = rho1_inf
  	)
  )

}
