//This program generates Figure 5: bias-corrected Q-LATEs of covariates against quantiles//

cd "C:\Users\yyd\Dropbox\RD_DistributionTreatment\EmpiricalApplication\" 
clear all
set more off
set matsize 8100
set maxvar 30000

capture program drop QLATE_bc_se
program define QLATE_bc_se, eclass 
	
	tempname qte qlist tau_u tau_u_bc biasT_plus biasT_minus dmdt_plus dmdt_minus biasRTau_u biasTTau_u BR2 BT2 c 
	tempname mi fTRplus fTRminus ifT_Rplus ifT_Rminus sigma2plus sigma2minus SEtau SEtau_bc Vtau qplus qminus
	tempvar Y2i	

	mat `qte'=J($qnum,1,.)	
	mat `qplus'=J($qnum,1,.)	
	mat `qminus'=J($qnum,1,.)	
	mat `qlist'=J($qnum,1,.)
	mat `tau_u'=J($qnum,1,.)
	mat `tau_u_bc'=J($qnum,1,.)
	mat `dmdt_plus' = J($qnum,1,.)
	mat `dmdt_minus' = J($qnum,1,.)
	mat `ifT_Rplus' = J($qnum,1,.)
	mat `ifT_Rminus' = J($qnum,1,.)
	mat `mi'=J($n,1,. )
	mat `Vtau' = J($qnum,1,.)
	
	local r=0
	forvalues q=$minq($grid)$maxq {
		tempvar T_c T_c2 ZT ZT2 RT_c ZRT_c Y2 kern_T kern_T1 kern_prod kern_prod1 kern_fTRplus kern_fTRminus w_sum 
		local r=`r'+1
		
		quiet qreg T R Z ZR  [iweight=kern_R], quantile(`q') 
		mat `qlist'[`r',1]=`q'
		mat `qte'[`r',1]=_b[Z]
		
		/* 
			Record q^+ and q^- for choosing T1 and T0 in V_pi^m
		*/
		mat `qplus'[`r',1] = _b[_cons] + _b[Z]
		mat `qminus'[`r',1] = _b[_cons]
		
		quiet gen `T_c'=T-`qminus'[`r',1]
		quiet replace `T_c'=T-`qplus'[`r',1] if R>=0
		quiet gen `T_c2'=`T_c'^2
		quiet gen `ZT'=`T_c'*Z
		quiet gen `ZT2'=`T_c2'*Z
		quiet gen `RT_c'=`T_c'*R
		quiet gen `ZRT_c'=`RT_c'*Z
		
		/*
			In the following, qte is set to be missing at the trimmed quantiles so they 
			will drop in the integration of tau_u to obtain pi. In addition, quantile 
			index is set to be zero, so in the double integration for variance 
			calculation, they will drop 
		*/
		if abs(`qte'[`r',1])<trim[`r',1]{
			mat `qte'[`r',1]=.
			mat `qlist'[`r',1]=0
			mat `qplus'[`r',1]=.
			mat `qminus'[`r',1]=.
		}
		
		/*
			Use uniform kernel
		*/
		quiet gen `kern_T'=0.5
		quiet replace `kern_T'=0 if abs(`T_c')>$h_T	
		quiet gen `kern_prod'=kern_R*`kern_T'
		
		quiet reg Y R `T_c' Z ZR `ZT' [iweight=`kern_prod']
		mat `tau_u'[`r',1]=_b[Z]/`qte'[`r',1]		

		/* 
			Estimate Variance of tau_u 
			Estimate sigma^{2,+/-} in V_tau 
		*/
		quiet gen `kern_T1'=0.5
		quiet replace `kern_T1'=. if abs(`T_c')>$h_T1	
		quiet gen `kern_prod1'=kern_R1*`kern_T1'
		
		quiet reg Y R `T_c' Z ZR `ZT' [iweight=`kern_prod1']
		quiet gen `Y2' = (Y- _b[_cons] - _b[Z] )^2
		quiet replace `Y2' = (Y- _b[_cons])^2 if R < 0
		quiet reg `Y2' R `T_c' Z ZR `ZT' [iweight = `kern_prod1']
		scalar `sigma2plus' = _b[_cons] + _b[Z]
		scalar `sigma2minus' = _b[_cons]
		
		/* 
			Parameters in V_tau. Estimate the joint density of R and T right above r0 
		*/
		quiet gen `kern_fTRplus'=0.25
		quiet replace `kern_fTRplus'=0 if abs(`T_c')>($h_fplus*sd_T)|abs(R)>($h_fplus*sd_R)
		quiet sum `kern_fTRplus' if Z == 1
		scalar `fTRplus' = r(mean)/($h_fplus^2*sd_R*sd_T)
		
		/* 
			Estimate the reciprocal of the conditional density of T given R right 
			above r0. That is, ifT_Rplus = 1/f_{TR}^+ 
		*/
		mat `ifT_Rplus'[`r',1] = fRplus/`fTRplus'		
		
		/* 
			Estimate the joint density of R and T right below r0
		*/
		quiet gen `kern_fTRminus'=0.25
		quiet replace `kern_fTRminus'=0 if abs(`T_c')>($h_fminus*sd_T)|abs(R)>($h_fminus*sd_R)
		quiet sum `kern_fTRminus' if Z == 0
		scalar `fTRminus' = r(mean)/($h_fminus^2*sd_R*sd_T)
		
		/* 
			Estimate the reciprocal of the conditional density of T given R right 
			below r0. That is, ifT_R_minus = 1/f_{TR}^-
		*/
		mat `ifT_Rminus'[`r',1] = fRminus/`fTRminus'		

		mat `Vtau'[`r',1] = 2*(`sigma2plus'*`ifT_Rplus'[`r',1] + `sigma2minus'*`ifT_Rminus'[`r',1])/(fR*(`qte'[`r',1])^2)

		/* 
			Estimate the bias of the estimated tau_u
		*/
		quiet qreg T R R2 Z ZR ZR2 [iweight=kern_R1], quantile(`q') 
		scalar `biasT_minus'=2*_b[R2]*(-1/12)
		scalar `biasT_plus'=2*(_b[R2]+_b[ZR2])*(-1/12)	

		quiet reg Y R R2 `T_c' `T_c2' `RT_c' Z ZR `ZT' ZR2 `ZT2' `ZRT_c' [iweight=`kern_prod1']
		scalar `BR2'=2*(_b[ZR2])*(-1/12) 
		scalar `BT2'=2*(_b[`ZT2'])*(1/6)
		mat `dmdt_minus'[`r',1] =_b[`T_c']
		mat `dmdt_plus'[`r',1] =_b[`T_c']+_b[`ZT']	
				
		scalar `biasRTau_u'=(`BR2'+`biasT_plus'*(`dmdt_plus'[`r',1]-`tau_u'[`r',1])-`biasT_minus'*(`dmdt_minus'[`r',1]-`tau_u'[`r',1]))/`qte'[`r',1]
		scalar `biasTTau_u'= `BT2'/`qte'[`r',1]
		
		mat `tau_u_bc'[`r',1]=`tau_u'[`r',1] - $h_R^2*`biasRTau_u' - $h_T^2*`biasTTau_u'
	}
	 
	svmat `Vtau'
	
	quiet gen `SEtau'= sqrt(`Vtau'1/($n*$h_R*$h_T))	
	
	if $rho_1 >= 1 {
		local contCtau = 37.5*($rho_1/3 - 0.25)
	}
	if $rho_1 < 1 {
		local contCtau = 3.125*$rho_1^3
	}
	
	quiet gen `SEtau_bc' = sqrt( 1 + $rho_1^6*9.765625  + $rho_1*`contCtau'   )*`SEtau'
	mkmat `SEtau_bc' 
	
	ereturn mat qlist=`qlist'	
	ereturn mat qte=`qte'
	ereturn mat tau_u_bc=`tau_u_bc'
	ereturn mat SEtau_bc=`SEtau_bc'
	
end 
	
***************************
***************************
use ***.dta

//Treatment//
	gen T=log_CapitalSurplus

//Covariates tested are bank age (1905-YearofEstablishment) and county characteristics including pct_blackpop, pct_farmland per square miles, and log(mfgout_percapita) for the year 1900 // 
	gen Y=1905-YearofEstablishment


global minq=0.1
global maxq=0.5
global grid=0.01
global qnum=($maxq-$minq)/$grid

	quiet sum R
	global n=r(N)	
	scalar sd_R=r(sd)
	quiet sum R if R>=0
	scalar n_plus=r(N)
	quiet sum R if R<0
	scalar n_minus=r(N)	
	quiet sum T
	scalar sd_T=r(sd)	

// Use the optimal bandwidth generated by RD_conT_OptimalBW_6.do (c=4.5) //
// For tau_u, the mean AMSE optimal bandwidth across quantiles for R is 1220.62 and for T is 0.5664//
// For pi, the AMSE optimal bandwidth is the following//
	global h_RPi=1462.7571
	global h_TPi=0.44081522
	
//rho_0 is the ratio of the optimal main bandwidth to the preliminary bandwidth used to estimate the trimming thresholds//
	global rho_0=4/3
	
//rho_1 is the ratio of the optimal main bandwidth to the bandwidth used to estimate biases and variances// 
	global c = 4.5
	global rho_1= $h_RPi/($c*$n^(-1/8)*sd_R)		

//h_R and h_T are the bandwidths for the main estimation//
	global h_R=$h_RPi
	global h_T=$h_TPi
	
//h_R0 is the preliminary bandwidth for generating the triming parameter//
	global h_R0=$h_R/$rho_0
	
//h_R1 and h_T1 the bandwidths for estimating the biases and variances (h/b) //
	global h_R1=$h_R/$rho_1
	global h_T1=$h_R1*sd_T/sd_R
	
//use uniform kernel//
	gen kern_R0=0.5
	replace kern_R0=0 if abs(R)>$h_R0	
	gen kern_R=0.5
	replace kern_R=0 if abs(R)>$h_R	
	quiet gen kern_R1=0.5
	replace kern_R1=0 if abs(R)>$h_R1	
	
//generate trimming threshold//
			mat trim=J($qnum,1,.)
			local r=0
			forvalues q=$minq($grid)$maxq {
			local r=`r'+1
			quiet qreg T R Z ZR  [iweight=kern_R0], quantile(`q') 
			mat trim[`r',1]=_se[Z]*1.96
			}
	
//Silverman rule-of-thumb bandwidth for a unifrom kernel for density estimation//
	global h_fR = 1.843*$n^(-0.2)*sd_R   
	gen kern_fR=0.5
	replace kern_fR=0 if abs(R)>$h_fR	
	quiet sum kern_fR
	scalar fR = r(mean)/$h_fR
	
//Estimate the density of R right above r0//
//Silverman rule-of-thumb bandwidth for a unifrom kernel//
	global h_fplus = 0.7344*n_plus^(-1/6)   
	quiet replace kern_fR=0.5
	quiet replace kern_fR=0 if abs(R)> ($h_fplus*sd_R)	
	quiet sum kern_fR if Z == 1
	scalar fRplus = r(mean)/($h_fplus*sd_R)
	
//Estimate the density of R right below r0//
//Silverman rule-of-thumb bandwidth for a unifrom kernel//		
	global h_fminus = 0.7344*n_minus^(-1/6)   
	quiet replace kern_fR=0.5
	quiet replace kern_fR=0 if abs(R)>$h_fminus*sd_R
	quiet sum kern_fR if Z == 0
	scalar fRminus = r(mean)/($h_fminus*sd_R)

QLATE_bc_se
		
mat qlist=e(qlist)
mat qte=e(qte)
mat tau_bc=e(tau_u_bc)
mat SEtau_bc=e(SEtau_bc)

svmat qlist, name(qlist_a)
svmat qte, name(qte_a)
svmat tau_bc, name(tau_bc_a)
svmat SEtau_bc, name (SEtau_bc_a)

quiet gen cil_a=tau_bc_a1-1.96*SEtau_bc_a1
quiet gen ciu_a=tau_bc_a1+1.96*SEtau_bc_a1

*note ylabel needs to set accordingly for different outcome variables
twoway (rarea cil_a ciu_a qlist_a1 if qte_a1!=., color(gs12)) (line tau_bc_a1 qlist_a1 if qte_a1!=., lcolor(blue) lwidth(medium)) (function y=0 if qte_a1!=., range(qlist_a1) clcolor(red) clpattern(dash) clwidth(medium)), ylabel(-3(1)3) title(Bank age) xtitle("Quantiles") ytitle("Bias corrected Q-LATE") legend(off) graphregion(color(white) ) plotregion(color(white)  style(none)) 
graph save "JASA_Tables&Figures\QLATE_bc_BankAge.gph", replace 


replace Y=pct_blackpop

QLATE_bc_se

mat qlist=e(qlist)
mat qte=e(qte)
mat tau_bc=e(tau_u_bc)
mat SEtau_bc=e(SEtau_bc)

svmat qlist, name(qlist_p)
svmat qte, name(qte_p)
svmat tau_bc, name(tau_bc_p)
svmat SEtau_bc, name (SEtau_bc_p)

quiet gen cil_p=tau_bc_p1-1.96*SEtau_bc_p1
quiet gen ciu_p=tau_bc_p1+1.96*SEtau_bc_p1

twoway (rarea cil_p ciu_p qlist_p1 if qte_p1!=., color(gs12)) (line tau_bc_p1 qlist_p1 if qte_p1!=., lcolor(blue) lwidth(medium)) (function y=0 if qte_p1!=., range(qlist_p1) clcolor(red) clpattern(dash) clwidth(medium)),  ylabel(-0.6(0.2)0.6) title(Black population (%)) xtitle("Quantiles") ytitle("Bias corrected Q-LATE") legend(off) graphregion(color(white) ) plotregion(color(white)  style(none)) 
graph save "JASA_Tables&Figures\QLATE_bc_BlackPop.gph", replace

replace Y=pct_farmland

QLATE_bc_se

mat qlist=e(qlist)
mat qte=e(qte)
mat tau_bc=e(tau_u_bc)
mat SEtau_bc=e(SEtau_bc)

svmat qlist, name(qlist_f)
svmat qte, name(qte_f)
svmat tau_bc, name(tau_bc_f)
svmat SEtau_bc, name (SEtau_bc_f)

quiet gen cil_f=tau_bc_f1-1.96*SEtau_bc_f1
quiet gen ciu_f=tau_bc_f1+1.96*SEtau_bc_f1

*note ylabel needs to set accordingly for different outcome variables
twoway (rarea cil_f ciu_f qlist_f1 if qte_f1!=., color(gs12)) (line tau_bc_f1 qlist_f1 if qte_f1!=., lcolor(blue) lwidth(medium)) (function y=0 if qte_f1!=., range(qlist_f1) clcolor(red) clpattern(dash) clwidth(medium)),  ylabel(-0.6(0.2)0.6)  title(Farmland (%)) xtitle("Quantiles") ytitle("Bias corrected Q-LATE") legend(off) graphregion(color(white) ) plotregion(color(white)  style(none)) 
graph save "JASA_Tables&Figures\QLATE_bc_Farmland.gph", replace 

replace Y=log(mfgout_percapita)

QLATE_bc_se

mat qlist=e(qlist)
mat qte=e(qte)
mat tau_bc=e(tau_u_bc)
mat SEtau_bc=e(SEtau_bc)

svmat qlist, name(qlist_m)
svmat qte, name(qte_m)
svmat tau_bc, name(tau_bc_m)
svmat SEtau_bc, name (SEtau_bc_m)

quiet gen cil_m=tau_bc_m1-1.96*SEtau_bc_m1
quiet gen ciu_m=tau_bc_m1+1.96*SEtau_bc_m1

twoway (rarea cil_m ciu_m qlist_m1 if qte_m1!=., color(gs12)) (line tau_bc_m1 qlist_m1 if qte_m1!=., lcolor(blue) lwidth(medium)) (function y=0 if qte_m1!=., range(qlist_m1) clcolor(red) clpattern(dash) clwidth(medium)), ylabel(-3(1)3) title(Log(manufacturing output)) xtitle("Quantiles") ytitle("Bias corrected Q-LATE") legend(off) graphregion(color(white) ) plotregion(color(white)  style(none)) 
graph save "JASA_Tables&Figures\QLATE_bc_Manufacturing.gph", replace 

cd JASA_Tables&Figures

gr combine "QLATE_bc_BankAge.gph" "QLATE_bc_BlackPop.gph" "QLATE_bc_Farmland.gph" "QLATE_bc_Manufacturing.gph", graphregion(color(white) margin(0) style(none)) plotregion(color(white) margin(0 0 0 0) style(none)) imargin(0 0 0 0)  
graph save "Fig5_falsification_test_1st_moment.gph", replace

