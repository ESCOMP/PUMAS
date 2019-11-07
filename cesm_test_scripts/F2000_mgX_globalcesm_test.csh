#!/bin/csh -vx

# A simplified script for running F2000climo global test cases in 
# **Github CAM/CESM**

# Edit the variables below this comment:

set CESMDIR    = /glade/work/katec/mg3work/ESCOMP_cam_devolpment_ktfork_pumas
set MG_vers    = mg3
set COMPSET    = F2000climo
set CASE       = gittest.$COMPSET.$MG_vers.002
set CASEDIR    = /glade/work/$USER/mg3work/cesm_cases/$CASE

set RESUBMIT   = 0
set STOP_N     = 6
set STOP_OPTION = nmonths
set REST_N     = 6
set REST_OPTION = nmonths
set DEBUG_TF   = TRUE
set ARCHIVE_TF = TRUE

set PBSACC     = P93300606
set TIMEWALL   = 12:00:00
set QUEUE      = regular

set STARTDATE  = 0001-01-01
set GETREF_TF  = TRUE
set REFCASE    = f.e21.f09_f09_F2000climo.cesm2_1.004
set REFDATE    = 0011-01-01-00000
set REFDIR     = /glade/u/home/katec/archive/{$REFCASE}/rest/{$REFDATE}

###! Add all user specific namelist changes below in the form of 
###! namelist_var = new_namelist_value 
cat <<EOF >! user_nl_cam.TMP
&camexp
  history_budget = .true.
  print_energy_errors = .false.
  macrop_scheme = 'CLUBB_SGS'
  cld_macmic_num_steps = 2
  micro_mg_do_hail = .true.
  micro_mg_do_graupel = .false.
  fincl1  = 'CDNUMC','ICLDIWP','AQRAIN','AQSNOW','ANRAIN','ANSNOW','FREQS','FREQR',
           'NUMICE','NUMLIQ','TMNUMICE','VTRMI','UMS','UMR',
           'LWCF','SWCF', 'ICLDTWP', 'CONCLD','ADRAIN','ADSNOW',
           'QCSEVAP', 'QISEVAP', 'CMELIQ', 'CMEIOUT', 'EVAPPREC', 'EVAPSNOW',
           'PRCO','PRAO','PRECT','PRECC','PRECL','NUMLIQ','GCLDLWP','TGCLDLWP',
           'CME','REL','QRAIN','QSNOW','REI','LWC','IWC','PE','PEFRAC',
	   'Q_qneg3','CLDICE_qneg3'
  nhtfrq  = 0
  mfilt   = 1
  avgflag_pertape = 'A'
  dust_emis_fact         = 0.7D0
  soil_erod_file         = '/glade/u/home/mahowald/dst_source2x2tunedcam6-2x2-04062017.nc'
  ext_frc_cycle_yr               = 2000
  srf_emis_cycle_yr              = 2000
/
EOF

cat <<EOF >! user_nl_clm.TMP
&clm_inparm
  finidat = '/glade/u/home/katec/archive/${REFCASE}/rest/${REFDATE}/${REFCASE}.clm2.r.${REFDATE}.nc'
  fsurdat = '/glade/work/dlawren/cesm_code/surfacedatasets/surfdata_0.9x1.25_78pfts_simyr1850_c170412.FWS.nc'
  use_init_interp = .true.
  suplnitro = 'NONE'
/
EOF

# You probably shouldn't need to edit anything below this line:

set RESOLUTION = f09_f09_mg17
set DEF_TASKS  = 1080

## create new case
${CESMDIR}/cime/scripts/create_newcase --case $CASEDIR --res $RESOLUTION --compset $COMPSET --run-unsupported

mv user_nl_cam.TMP $CASEDIR
mv user_nl_clm.TMP $CASEDIR
cd $CASEDIR

## set PE Layout
./xmlchange NTASKS_ATM=$DEF_TASKS,NTASKS_LND=$DEF_TASKS,NTASKS_ICE=$DEF_TASKS,NTASKS_OCN=$DEF_TASKS,NTASKS_CPL=$DEF_TASKS,NTASKS_GLC=1,NTASKS_ROF=$DEF_TASKS,NTASKS_WAV=$DEF_TASKS

## set other options
./xmlchange RUN_STARTDATE=$STARTDATE
./xmlchange DEBUG=$DEBUG_TF 
./xmlchange DOUT_S=$ARCHIVE_TF
./xmlchange PROJECT=$PBSACC
./xmlchange RESUBMIT=$RESUBMIT
./xmlchange STOP_OPTION=$STOP_OPTION
./xmlchange STOP_N=$STOP_N
./xmlchange REST_OPTION=$REST_OPTION
./xmlchange REST_N=$REST_N
./xmlchange JOB_WALLCLOCK_TIME=$TIMEWALL
./xmlchange JOB_QUEUE=$QUEUE
./xmlchange GET_REFCASE=$GETREF_TF
./xmlchange RUN_REFCASE=$REFCASE
./xmlchange RUN_REFDATE=$REFDATE
./xmlchange RUN_REFDIR=$REFDIR
./xmlchange RUN_REFTOD=00000

## Change CAM Configure here. 
./xmlchange --append CAM_CONFIG_OPTS="-microphys "$MG_vers" "
## Change CLM configuration
./xmlchange CLM_BLDNML_OPTS="-bgc bgc -crop"


## call case.setup (reset just in case)
./case.setup --reset

# CLOBBER USER_NL_CAM
mv user_nl_cam.TMP user_nl_cam
# CLOBBER USER_NL_CLM
mv user_nl_clm.TMP user_nl_clm

## call case.build
setenv PBS_ACCOUNT $PBSACC
qcmd -- ./case.build

## Submit job
./case.submit
