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
