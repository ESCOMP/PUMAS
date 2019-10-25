#!/bin/csh -vx

# A simplified script for running SCAM tests in **Github CAM/CESM**

# Edit the variables below this comment:

set CESMDIR = /glade/work/katec/mg3work/ESCOMP_cam_devolpment_ktfork_pumas
set CASEDIR = /glade/work/katec/mg3work/scam6
set CASETITLE = cleanscript_ESCOMP
#available: mpace arm97 arm95 togaII gateIII twp06 sparticus
set IOP = arm97 
set DEBUG_TF = 'TRUE'
set MG_vers = mg3
set ARCHIVE_TF = 'FALSE'
set PBSACC = P93300642

###! Add all user specific namelist changes below in the form of 
###! namelist_var = new_namelist_value 
cat >> user_nl_cam.TMP << EOF
  scm_relaxation               = .true.
  ncdata='/glade/p/cesmdata/cseg/inputdata/atm/cam/inic/gaus/cami_0000-01-01_64x128_L32_c170510.nc' 
  scm_relaxation         = .false.
  print_energy_errors = .true.
  scm_use_obs_uv = .true.
  scm_use_obs_T  = .true.
  use_gw_rdg_beta = .false.
  mfilt          = 3000
  nhtfrq         = 1
  history_budget = .true.
  macrop_scheme  = 'CLUBB_SGS'
  cld_macmic_num_steps = 2
  micro_mg_do_hail = .false.
  micro_mg_do_graupel = .true.
  fincl1  = 'CDNUMC','ICLDIWP','AQRAIN','AQSNOW','ANRAIN','ANSNOW','FREQS','FREQR',
           'AQGRAU','ANGRAU','NUMICE','NUMLIQ','TMNUMICE','VTRMI','UMS','UMR',
           'LWCF','SWCF', 'ICLDTWP', 'CONCLD','ADRAIN','ADSNOW',
           'QCSEVAP', 'QISEVAP', 'CMELIQ', 'CMEIOUT', 'EVAPPREC', 'EVAPSNOW',
           'PRCO','PRAO','PRECT','PRECC','PRECL','NUMLIQ','GCLDLWP','TGCLDLWP',
           'CME','REL','QRAIN','QSNOW','REI','LWC','IWC','PE','PEFRAC',
           'QGSEDTEN','PRACGO','PGRACSO','PRDGO','PSACRO','MNUCCRO','PSACWGO',
           'PGSACWO','VTRMI','VTRMC','UMR','UMS'
EOF


# You probably shouldn't need to edit anything below this line:

set COMPSET=FSCAM
set COMPILER=intel
set RES=T42_T42
set IOPNAME = scam_$IOP
set CASENAME=${CASETITLE}.${MG_vers}.${COMPSET}.${IOP}
set MODSDIR = $CESMDIR/cime_config/usermods_dirs

###create and run case

$CESMDIR/cime/scripts/create_newcase --compset $COMPSET  --res $RES --compiler $COMPILER --case $CASEDIR/$CASENAME  --user-mods-dir ${MODSDIR}/${IOPNAME} --run-unsupported

mv user_nl_cam.TMP $CASEDIR/$CASENAME
cd  $CASEDIR/$CASENAME

./xmlchange DEBUG=$DEBUG_TF
./xmlchange --append CAM_CONFIG_OPTS="-microphys "$MG_vers" "
./xmlchange DOUT_S=$ARCHIVE_TF
./xmlchange PROJECT=$PBSACC

./case.setup

# CLOBBER USER_NL_CAM
mv user_nl_cam.TMP user_nl_cam

setenv PBS_ACCOUNT $PBSACC
# Always clean over and build
./case.build --clean-all
qcmd -- ./case.build 

### Submit to Queue
./case.submit



 
