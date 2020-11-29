; test_tplot_rbsp_spec_xspec crib
;
; IDL crib for manually calculating XSPEC from efw l1_int b1 CIT data
; using rbsp_spec.pro and rbsp_xspec.pro
;
; version 1.0, Kris Kersten, UMN, June 2012
;			email: kris.kersten@gmail.com


; initialize RBSP environment
rbsp_efw_init


;-------------------------------------------------------------------------------
;
;  RBSPA
;
;-------------------------------------------------------------------------------

probe = 'a'

date = '2012-01-19'
duration = 1
timespan, date, duration

integration=1

datatype = ['eb1', 'mscb1']

get_support_data = 0

rbsp_load_efw_waveform, probe=probe, datatype=datatype, type=type, $
		get_support_data=get_support_data, integration=integration


split_vec, 'rbsp?_efw_eb1', suffix='_'+['E12', 'E34', 'E56']
split_vec, 'rbsp?_efw_mscb1', suffix='_'+['U', 'V', 'W']


tplot_options, 'xmargin', [ 20., 15.] ; increase margins so we don't lose labels
options,'*','labflag',-1


rbsp_xspec,'rbspa_efw_mscb1_W','rbspa_efw_eb1_E12', $
		/median_subtract, /nan_fill_gaps
rbsp_spec,'rbspa_efw_mscb1_W', $
		/median_subtract, /nan_fill_gaps
rbsp_spec,'rbspa_efw_eb1_E12', $
		/median_subtract, /nan_fill_gaps


tplot,[ 'rbspa_efw_mscb1_W','rbspa_efw_eb1_E12', $
	'rbspa_efw_mscb1_W_SPEC', 'rbspa_efw_eb1_E12_SPEC', $
	'rbspa_efw_mscb1_W_X_rbspa_efw_eb1_E12_coh', $
	'rbspa_efw_mscb1_W_X_rbspa_efw_eb1_E12_pha' ]


; set up interesting time limits
tlimit,'2012-01-19/14:53:30','2012-01-19/14:57:00' ; white noise?
tlimit,'2012-01-19/16:53:30','2012-01-19/17:27:00' ; dynamic range
tlimit,'2012-01-19/17:35:30','2012-01-19/17:41:10' ; frequency sweeps
tlimit,'2012-01-19/18:41:15','2012-01-19/18:53:15' ; plasma simulator
tlimit,'2012-01-19/20:50:14','2012-01-19/20:50:20' ; 1kHz burst?




;-------------------------------------------------------------------------------
;
;  RBSPB
;
;-------------------------------------------------------------------------------

probe = 'b'


date = '2012-01-26'
duration = 1
timespan, date, duration

integration=1
datatype = ['eb1', 'mscb1']
get_support_data = 0

rbsp_load_efw_waveform, probe=probe, datatype=datatype, type=type, $
		get_support_data=get_support_data, integration=integration

split_vec, 'rbspb_efw_eb1', suffix='_'+['E12', 'E34', 'E56']
split_vec, 'rbspb_efw_mscb1', suffix='_'+['U', 'V', 'W']

tplot_options, 'xmargin', [ 20., 15.]
options,'*','labflag',-1


rbsp_xspec,'rbspb_efw_mscb1_W','rbspb_efw_eb1_E12', $
		/median_subtract, /nan_fill_gaps
rbsp_spec,'rbspb_efw_mscb1_W', $
		/median_subtract, /nan_fill_gaps
rbsp_spec,'rbspb_efw_eb1_E12', $
		/median_subtract, /nan_fill_gaps

tplot,[ 'rbspb_efw_mscb1_W','rbspb_efw_eb1_E12', $
	'rbspb_efw_mscb1_W_SPEC', 'rbspb_efw_eb1_E12_SPEC', $
	'rbspb_efw_mscb1_W_X_rbspb_efw_eb1_E12_coh', $
	'rbspb_efw_mscb1_W_X_rbspb_efw_eb1_E12_pha' ]


; set up interesting time limits
tlimit,'2012-01-26/16:06:10','2012-01-26/16:31:00' ; dynamic range
tlimit,'2012-01-26/19:43:00','2012-01-26/19:51:20' ; frequency sweeps
tlimit,'2012-01-26/16:55:30','2012-01-26/17:06:30' ; plasma simulator



;detailed look
tlimit,'2012-01-26/19:46:11.100','2012-01-26/19:46:11.5' ; frequency sweeps
tlimit,'2012-01-26/19:46:11.000','2012-01-26/19:46:13.0' ; frequency sweeps
tlimit,'2012-01-26/19:49:00.000','2012-01-26/19:51:00.0' ; frequency sweeps



get_data,'rbspb_efw_mscb1_W_X_rbspb_efw_eb1_E12_coh',data=dat2
goo2 = where(dat2.x ge time_double('2012-01-26/19:46:10'))
plot,dat2.v,dat2.y[goo2[0],*]

get_data,'rbspb_efw_mscb1_W_X_rbspb_efw_eb1_E12_pha',data=dat
goo = where(dat.x ge time_double('2012-01-26/19:46:11.450'))
plot,dat.v,dat.y[goo[61],*],psym=2


ylim,'rbspb_efw_mscb1_W_X_rbspb_efw_eb1_E12_coh',0,8000
ylim,'rbspb_efw_mscb1_W_X_rbspb_efw_eb1_E12_pha',0,8000





end ; test_rbsp_efw_phase.pro





