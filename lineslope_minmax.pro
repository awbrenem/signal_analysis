;Calculate the best fit linear line and minimum and maximum slopes based on 
;LINFIT and the input error bars. The min and max slopes are calculated based on the 
;returned standard deviation of the slope and intercept
;Returns useful values

;*********
;******SHOULD UPDATE THIS TO USE fitexy.pro WHICH ALLOWS UNCERTAINTIES IN BOTH X, Y DIRECTIONS
;*********


function lineslope_minmax,xvals,yvals,yerr,xerr=xerror


  coefs = LINFIT(xvals,yvals,MEASURE_ERRORS=yerr,/DOUBLE,chisqr=chs,sigma=siggy,prob=prob)

stop

  fitexy,xvals,yvals,a,b,x_sig=xerror,y_sig=yerr


;       FITEXY, x, y, A, B, X_SIG= , Y_SIG= , [sigma_A_B, chi_sq, q, TOL=]
;
; INPUTS:
;       x = array of values for independent variable.
;       y = array of data values assumed to be linearly dependent on x.
;
; REQUIRED INPUT KEYWORDS:
;       X_SIGMA = scalar or array specifying the standard deviation of x data.
;       Y_SIGMA = scalar or array specifying the standard deviation of y data.
;
; OPTIONAL INPUT KEYWORD:
;       TOLERANCE = desired accuracy of minimum & zero location, default=1.e-3.
;
; OUTPUTS:
;       A_intercept = constant parameter result of linear fit,
;       B_slope = slope parameter, so that:
;                       ( A_intercept + B_slope * x ) approximates the y data.











  ;Calculate reduced chi-sqr
  nparams = 2.  ;y=mx + b
  ndegfreedom = n_elements(xvals)
  v = ndegfreedom - nparams
  chs_reduced = chs/v
  print,'*******************'
  print,'Chisqr = ',chs
  print,'Chisqr reduced = ',chs_reduced
  print,'P-value = ',prob
  print,'Prob of good fit = ',1-prob
  print,'Coeffs = ',coefs
  print,'Sigma = ',siggy
  print,'*******************'



  ;Best fit line
  fitline = coefs[1]*xvals + coefs[0]
  ;Max slope line 
  fitlinemax = (coefs[1]+siggy[1])*xvals + (coefs[0]-siggy[0])
  ;Min slope line 
  fitlinemin = (coefs[1]-siggy[1])*xvals + (coefs[0]+siggy[0])



  return,{xvals:xvals,yvals:yvals,yerr:yerr,$
          chisqr:chs,chisqr_reduced:chs_reduced,$
          p_value:prob,prob_of_good_fit:1-prob,$
          coeff:coefs,sigma:siggy,$
          fitline:fitline,fitlinemax:fitlinemax,fitlinemin:fitlinemin}


end












;**************************
;Below code is obsolete
;**************************


;;delta-Intercept values
;dI = 0.1*coefs[0]
;baddat = 0
;bbq = 0.
;
;while not baddat do begin
;
;  A=[coefs[0]+(dI*bbq), coefs[1]]
;  fita=[0,1]
;
;  coef_line = LMFIT(xvals,yvals,A,measure_errors=err,/double, $
;  FITA=fita,FUNCTION_NAME = 'myfunct_linear',CHISQ=chisqrtmp)
;
;  if bbq eq 0 then bestfitY = coef_line;
;
;  oplot,xvals,coef_line
;  ;oplot,xvals,fitline,color=250
;
;  yv = coefs[0] + xvals*coefs[1]
;  oplot,xvals,yv,color=50,linestyle=2,thick=3
;
;
;  ;Stirling's approximation
;  n = v/2.
;  gamma_stirling = sqrt(2.*!pi)*exp(-n)*n^(n-0.5)*(1.0+0.0833/n);
;
;  x2 = chisqrtmp
;  x2_red = x2/(ndegfreedom-nparams)
;
;
;  ;Stop when this value > 0.05
;  print,'chisq = ', x2
;  print,'chisq reduced = ',x2_red
;  print,'p-value = ',1 - CHISQR_PDF(x2, v)
;  print,'prob of good fit = ',CHISQR_PDF(x2, v)
;  print,'*********'
;
;  if chisqr_pdf(x2,v) ge 0.95 then begin
;    if dI gt 0 then begin
;      minfitY = coef_line
;      ;change to decrementing the y-intercept
;      dI = -0.1*coefs[0]
;      bbq = 0
;    endif else begin
;      maxfitY = coef_line
;      baddat = 1
;    endelse
;  endif
;  bbq++
;
;endwhile

;stop
;

;plot,xvals,bestfitY
;oplot,xvals,maxfitY
;oplot,xvals,minfitY

;stop
