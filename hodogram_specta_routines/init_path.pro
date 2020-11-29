;+
;*****************************************************************************************
;
;  PROCEDURE :  init_path 
;  PURPOSE  :   creates the !data sys variable that contains the paths to various
;				data
;				
;
;  CALLED BY:   
;               
;
;  CALLS:
;               
;
;  REQUIRES:    
;               
;
;  INPUT:
;               
;
;  EXAMPLES:    
;               
;
;  KEYWORDS:    
;               
;
;   CHANGED:  1)  NA [MM/DD/YYYY   v1.0.0]
;
;   NOTES:      
;               
;
;   CREATED:  01/21/2011
;   CREATED BY:  Aaron W. Breneman
;    LAST MODIFIED:  MM/DD/YYYY   v1.0.0
;    MODIFIED BY: Aaron W. Breneman
;
;*****************************************************************************************
;-
pro init_path

pre = '/Users/aaronbreneman/Desktop/'



;GENERAL DATA PATH
general = pre + 'code/Aaron/datafiles/'


;STEREO----------------------------------------------

pre2 = 'code/Aaron/datafiles/stereo/'
pre2 = pre + pre2


mag_l1 = pre2 + 'impact/level1/mag/'
plastic_l1 = pre2 + 'plastic/level1/'
plastic_l2 = pre2 + 'plastic/level2/'
eventdens = pre2 + 'event_density'
reducedlist = pre2 + 'reduced_list'
wavequality = pre2 + 'wave_quality'
tdsmax = pre2 + 'tdsmax'
tplot = pre2 + 'tplot/' 
burstwaves = pre2 + 'burstwaveforms/'

stereo = {mag_l1:mag_l1,$
			plastic_l1:plastic_l1,$
			plastic_l2:plastic_l2,$
			eventdens:eventdens,$
			reducedlist:reducedlist,$
			wavequality:wavequality,$
			tdsmax:tdsmax,$
			tplot:tplot,$
			burstwaves:burstwaves}


;----------------------------------------------------
;LIGHTNING BOLT -------------------------------------

pre2 = 'code/Aaron/datafiles/lbolt/'
pre2 = pre + pre2

TDSascii = pre2 + 'TDS/LEVEL_1/ASCII/'
TDSidl = pre2 + 'TDS/LEVEL_1/IDL/'
flashdat = pre2 + 'lightning_flash_data/'
link1 = pre2 + 'TM/LINK1/'
link2 = pre2 + 'TM/LINK2/'
tplot = pre2 + 'tplot/'
calibration = pre + 'Research/UMN_rocket/rowland/LBOLT/FINALCALS/DATA/RAWCALDATA/PROGRAMS/'
ionosonde = pre + 'Research/UMN_rocket/IONOSONDE/'
nldn = pre2 + 'lightning_flash_data/'



lbolt = {TDSascii:TDSascii,$
		 TDSidl:TDSidl,$
		 flashdat:flashdat,$
		 link1:link1,$
		 link2:link2,$
		 tplot:tplot,$
		 calibration:calibration,$
		 ionosonde:ionosonde,$
		 nldn:nldn}

;-----------------------------------------------------
;LANL DATA -------------------------------------------

pre2 = 'code/Aaron/datafiles/lanl/'
pre2 = pre + pre2

lanl = pre2



;-----------------------------------------------------
;WIND ------------------------------------------------

pre2 = 'code/Aaron/datafiles/wind/'
pre2 = pre + pre2

wind = pre2

;-----------------------------------------------------

;SCRIPTS

pre2 = 'code/Aaron/scripts/'
pre2 = pre + pre2
scripts = pre2


finalstr = {general:general,$
			stereo:stereo,$
			lbolt:lbolt,$
			lanl:lanl,$
			wind:wind,$
			scripts:scripts}

defsysv,'!data', finalstr
print,'!data system variable set'




end