;Linear fit y=mx+b for LMFIT called from dispersion_fitter_for_microbursts.pro

FUNCTION myfunct_linear, X, A

;  return, [A[0] + X*A[1], A[0], X*A[1]]
  return, [A[0] + X*A[1], 1., 1.]
END
