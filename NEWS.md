# TapeR 0.5.3

* in case of non-availability of apVar, stop Taper_FIT_LME.f() with appropriate
  error message, i.e. the value hold by element 'apVar' of the fitted lme-object

# TapeR 0.5.2

* improved calculation of approximate confidence intervals for diamter
  estimation in E_DHx_HmDm_HT.f() in case of sHT > 0. Especially changed 
  smooth.spline estimation of confidence bands from fixed degrees of freedom 
  to all.knots=TRUE to more precisely follow the data, especially at the 
  lower stem part. Also updated the return list to reflect the approximate MSE
  in all elements, not only in the confidence bounds (before, list elements 
  hold the information of expected taper curve, without the extra uncertainty
  due to measurement error in tree height). 

# TapeR 0.5.1

* added KOV_Mean and KOV_Pred to output of E_DHx_HmDm_HT.f()

# TapeR 0.5.0

* minor adjustments for submission to CRAN

# TapeR 0.4.3

* added function resVar() for more flexibility in residual error 
  assumptions, hence, changed and renamed the R0 interface: now called Rfn
  awaiting a list instead of logical

# TapeR 0.4.2

* included 'R0'-interpolation to E_HDx_* and Hx_root.f
* improved root-finding in E_Hdx_* so that error is avoided when taper
  curve is not monotonically decreasing, which might happen rarely in rather
  abnormal taper forms. The maximum height of searched diameter is returned.

# TapeR 0.4.1

* improved the XZ_BSPLINE.f()-function to make it faster. This 
  especially means that it does not return both X and Z matrices but only X or
  Z as they are only dependent on order and knots.

