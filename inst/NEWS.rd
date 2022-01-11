\name{NEWS}
\title{News for Package \pkg{TapeR}}

\section{Changes in version 0.5.0}{
  \itemize{
    \item minor adjustments for submission to CRAN
  }
}

\section{Changes in version 0.4.3}{
  \itemize{
    \item added function resVar() for more flexibility in residual error 
    assumptions, hence, changed and renamed the R0 interface: now called Rfn
    awaiting a list instead of logical
  }
}
\section{Changes in version 0.4.2}{
  \itemize{
    \item included 'R0'-interpolation to E_HDx_* and Hx_root.f
    \item improved root-finding in E_Hdx_* so that error is avoided when taper
    curve is not monotonically decreasing, which might happen rarely in rather
    abnormal taper forms. The maximum height of searched diameter is returned.
  }
}
\section{Changes in version 0.4.1}{
  \itemize{
    \item improved the XZ_BSPLINE.f()-function to make it faster. This 
    especially means that it does not return both X and Z matrices but only X or
    Z as they are only dependent on order and knots.
  }
}



