! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to deal with namelists for input
!
module namelists

  Use parameters
  Use header_data, only : mphys_id
  Use switches
  Use switches_bin

#if SHIPWAY_MICRO == 1
  ! Temporary for adding in 4a switches
  Use mphys_switches, only: iopt_act, option, aerosol_option        &
     , l_aaut, l_aacc, l_aevp, l_ased, iopt_tidy, l_warm            &
     , l_inuc, iopt_rcrit, iopt_inuc                                &
     , l_evaporation, l_rain, l_sed, l_boussinesq, diag_mu_option   &
     , l_sed_3mdiff, l_cons, l_abelshipway, l_condensation          &
     , l_active_inarg2000, l_oneway, l_newoptions 
  Use mphys_parameters, only: p1, p2, p3, sp1, sp2, sp3             &
     , max_step_length, max_sed_length
#endif

  implicit none
  
  namelist/mphys/num_h_moments, num_h_bins, mom_init, &
       h_names, mom_names, mom_units,num_aero_moments,num_aero_bins, &
       aero_mom_init, aero_N_init, aero_sig_init, aero_rd_init, aero_names
  
  namelist/control/dt, dg_dt, mphys_scheme, mphys_var &
       , wctrl, zctrl, tctrl, pctrl_z, pctrl_v, pctrl_T, ipctrl &
       , xctrl, lhf_ctrl, shf_ctrl, diaglevel, dgstart  

  namelist/case/input_file, l_input_file, ifiletype, icase

  namelist/switch/l_mphys, l_advect, l_diverge, l_pupdate &
       , l_fix_qv, l_nomphys_qv, l_noadv_qv, l_posadv_qv &
       , l_fix_theta, l_nomphys_theta, l_noadv_theta  &
       , l_noadv_hydrometeors, l_nodiv_hydrometeors, l_sediment &
       , isurface, l_noadv_aerosols, l_nodiv_aerosols, l_fix_aerosols &
       , l_sed_ult, l_diverge_advection, l_periodic_bound  &
       , l_force_positive

  logical :: iiwarm=.false.
  character(200) :: KiD_outdir=''
  character(200) :: KiD_outfile=''

  namelist/addcontrol/iiwarm, KiD_outdir, KiD_outfile  &
#if SHIPWAY_MICRO == 1
     ! Shipway 4A ...
     , option, l_evap, l_sed_3mdiff &
     , max_step_length, max_sed_length, diag_mu, l_rsource, l_raut &
     , l_evaporation, l_rain, l_sed, l_boussinesq, diag_mu_option   &
     , p1, p2, p3 &
     , sp1, sp2, sp3 &
     , l_abelshipway, l_cons &
     , l_coll_coal, l_cond, l_condensation, iopt_act &
     , aerosol_option, l_aaut, l_aacc, l_aevp, l_ased &
     , iopt_tidy, l_warm, l_inuc, iopt_rcrit   &
     , l_active_inarg2000, iopt_inuc, l_cu_cold, l_oneway, l_newoptions &
#endif
     ! Thompson 09...
     , l_reuse_thompson_lookup

  ! Namelist input...

  character(200) :: fileName=''
  character(200) :: fileNameIn=''
  character(200) :: fileNameOut=''
  character(200) :: namelistIn='namelists/input.nml'
  character(200) :: fileIn=''
  character(200) :: fileOut=''
  logical :: fexist

  namelist/namelistToUse/fileIn, fileOut

contains

  subroutine read_namelist
    !
    ! Read the namelists from file
    !

#if COMMANDLINE == 1
    ! This bit is F2003 - If your compiler doesnt support 
    ! it you need to comment the line out you can then specify 
    ! which namelist to use throught namelists/input.nml file 
    ! or else use command line processing that works with your 
    ! compiler (nearly all will do it but not in a portable way).
    write(*,*) 'Querying command line'
    CALL GET_COMMAND_ARGUMENT(1,fileNameIn)
    CALL GET_COMMAND_ARGUMENT(2,fileNameOut)
#endif

    if (trim(fileNameIn)=='')then  ! Not input at the command line 
                                   ! so use input.nml
#ifdef DEF_CASE
      write(namelistIn, '(A,A,A)') 'namelists/', DEF_CASE, '_input.nml'
#endif

      write(*,*) 'Unable to determine input file from command line, so querying ', trim(namelistIn), ' instead...'
      inquire(file=namelistIn, exist=fexist)
      if (fexist) then
        open(2, file=namelistIn)
        read(2, namelistToUse)
        close(2)
        write(*, namelistToUse)
      end if
      fileNameIn  = fileIn
      if (trim(fileOut)/='')fileNameOut = fileOut
    end if
    
    if (trim(fileNameIn)/='')fileName=fileNameIn

    write(6,*) 'Using namelist: ', trim(fileName)

    open(1, file=fileName)
!    rewind(1)
    read(1,mphys) 
!    rewind(1)
    read(1,case)
!    rewind(1)
    read(1,control) 
!    rewind(1)
    read(1,switch) 
!    rewind(1)
    read(1,addcontrol) 
    close(1)

    select case(mphys_scheme)
    case('lem2.4')
       imphys=imphys_lem2_4
       mphys_id='LEM2.4'
    case('tau_bin')
       imphys=imphys_tau_bin
       mphys_id='TAU_bin'
    case('thompson') ! NB same as thompson09
       imphys=imphys_thompson09
       mphys_id='thompson09'
    case('thompson09')
       imphys=imphys_thompson09
       mphys_id='thompson09'
    case('thompson06')
       imphys=imphys_thompson06
       mphys_id='thompson06'
    case('thompson07')
       imphys=imphys_thompson07
       mphys_id='thompson07' 
    case('morr_two_moment')
       imphys=imphys_morr_two_moment
       mphys_id='morr_two_moment'
    case('um7_3')
       imphys=imphys_um7_3
       mphys_id='um7_3'
    case('wsm6')
       imphys=imphys_wsm6
       mphys_id='wsm6'
    case('wdm6')
       imphys=imphys_wdm6
       mphys_id='wdm6'
    case('4A')
       imphys=imphys_4A
       mphys_id='4A'
    case('mg1_5')
       imphys=imphys_mg1_5
       mphys_id='mg1_5'
    case('mg2')
       imphys=imphys_mg2
       mphys_id='mg2'
    case('mg3')
       imphys=imphys_mg3
       mphys_id='mg3'
     case default
       print*, 'Mphys scheme not recognized: ', mphys_scheme
       print*, 'Did you mean:' 
       print*, '   lem2.4?'
       print*, '   um7_3?'
       print*, '   tau_bin?'
       print*, '   thompson09?'
       print*, '   thompson07?'
       print*, '   morr_two_moment?'
       print*, '   wsm6?'
       print*, '   4A?'
       print*, '   mg1_5?'
       print*, '   mg2?'
       print*, '(NB not all available in release version)'
       stop
    end select

    if (trim(input_file)=='')l_input_file=.False.

    if (.not. l_input_file)ifiletype=itest_case

  end subroutine read_namelist

end module namelists

