;Calculate the p-value. Low p-values indicate greater confidence that
;deviation from the null hypothesis is significant.
;https://en.wikipedia.org/wiki/Chi-squared_distribution
;However, this works in a weird way when considering lines slopes.
;If you have large error bars, then the best fit line may return a very
;low reduced chisqr. This will give a low prob of good fit b/c of the
;low chisq. So, keep adjusting the line slope until you have a chisq
;near(ish) unity. Eventually the prob of good fit metric will get to
;95%, which means that the adjusted line fit is straddling the edges
;of the error bars. This is what I'll take as my max/min line slope.


;err is the error bars on yvals

pro lineslope_minmax,xvals,yvals,err


;Use LINFIT to get an initial guess as to the slope and intercept, which
;are needed for LMFIT
coefs = LINFIT(xvals,yvals,MEASURE_ERRORS=err, $
/DOUBLE,chisqr=chs,sigma=siggy,prob=prob)

fitline = coefs[1]*xvals + coefs[0]
oplot,xvals,fitline,color=250

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
print,'*******************'


;Now we'll determine max and min slopes. Max(min) slope will have min(max)
;y-intercept. Using LMFIT, hold the intercept constant and find new best
;fit slope and chi-sq. Increment the slope until the chi-sq value is too low.


stop

;Calculate the p-value. Low p-values indicate greater confidence that
;deviation from the null hypothesis is significant.
;https://en.wikipedia.org/wiki/Chi-squared_distribution
;However, this works in a weird way when considering lines slopes.
;If you have large error bars, then the best fit line may return a very
;low reduced chisqr. This will give a low prob of good fit b/c of the
;low chisq. So, keep adjusting the line slope until you have a chisq
;near(ish) unity. Eventually the prob of good fit metric will get to
;95%, which means that the adjusted line fit is straddling the edges
;of the error bars. This is what I'll take as my max/min line slope.

;delta-Intercept values
dI = 0.1*coefs[0]
baddat = 0
bbq = 0.

while not baddat do begin

  A=[coefs[0]+(dI*bbq), coefs[1]]
  fita=[0,1]

  coef_line = LMFIT(xvals,yvals,A,measure_errors=err,/double, $
  FITA=fita,FUNCTION_NAME = 'myfunct_linear',CHISQ=chisqrtmp)

  if bbq eq 0 then bestfitY = coef_line

  oplot,xvals,coef_line
  ;oplot,xvals,fitline,color=250

  yv = coefs[0] + xvals*coefs[1]
  oplot,xvals,yv,color=50,linestyle=2,thick=3


  ;Stirling's approximation
  n = v/2.
  gamma_stirling = sqrt(2.*!pi)*exp(-n)*n^(n-0.5)*(1.0+0.0833/n)

  x2 = chisqrtmp
  x2_red = x2/(ndegfreedom-nparams)


  ;Stop when this value > 0.05
  print,'chisq = ', x2
  print,'chisq reduced = ',x2_red
  print,'p-value = ',1 - CHISQR_PDF(x2, v)
  print,'prob of good fit = ',CHISQR_PDF(x2, v)
  print,'*********'

  if chisqr_pdf(x2,v) ge 0.95 then begin
    if dI gt 0 then begin
      minfitY = coef_line
      ;change to decrementing the y-intercept
      dI = -0.1*coefs[0]
      bbq = 0
    endif else begin
      maxfitY = coef_line
      baddat = 1
    endelse
  endif
  bbq++

endwhile

stop


plot,xvals,bestfitY
oplot,xvals,maxfitY
oplot,xvals,minfitY

stop

end
