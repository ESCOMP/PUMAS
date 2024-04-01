module BOSS_utils

!--------------------------------------------------------------------------
!
! This module constaints parameter values to calculate BOSS process rates
! for warm rain
!
!--------------------------------------------------------------------------

implicit none
private
save

 public :: &
      BOSS_init, &
      BOSS_liq_autoconversion, &
      BOSS_accrete_cloud_water_rain, &
      BOSS_self_collection_rain, &
      BOSS_self_collection_cloud

 integer, parameter, public :: r8 = selected_real_kind(12)

! alternate threshold used for some in-cloud mmr
 real(r8), parameter :: icsmall = 1.e-8_r8
! Smallest mixing ratio considered in microphysics.
 real(r8), parameter :: qsmall = 1.e-18_r8

contains

subroutine BOSS_init(iautoq,pautoq,pautoN,paccq,psccN,pscrN,&
     log_a_auto_t1, b_auto_t1, log_mc_auto_inv, log_mr_auto_inv, &
     log_a_acc, b_acc_mc, b_acc_mr, log_a_sc_c, b_sc_c, log_a_sc_r, b_sc_r, &
     b_auto_t2_mc, b_auto_t2_mr, b_auto_t2_n, log_a_auto_t2, &
     pi,iulog)

  integer, intent(in) :: iautoq

  real(r8), intent(in) ::   log_a_auto_t1
  real(r8), intent(in) ::   b_auto_t1
  real(r8), intent(in) ::   log_mc_auto_inv
  real(r8), intent(in) ::   log_mr_auto_inv
  real(r8), intent(in) ::   log_a_acc
  real(r8), intent(in) ::   b_acc_mc
  real(r8), intent(in) ::   b_acc_mr
  real(r8), intent(in) ::   log_a_sc_c
  real(r8), intent(in) ::   b_sc_c
  real(r8), intent(in) ::   log_a_sc_r
  real(r8), intent(in) ::   b_sc_r

  !placeholder for iautoq = 2 and 3
  real(r8), intent(in) :: b_auto_t2_mc
  real(r8), intent(in) :: b_auto_t2_mr
  real(r8), intent(in) :: b_auto_t2_n
  real(r8), intent(in) :: log_a_auto_t2
  ! BOSS parameters
  real(r8), intent(out), allocatable, dimension(:) :: pautoq
  ! parameters for autoconversion (N)
  real(r8), intent(out), dimension(2) :: pautoN
  ! parameters for accretion (q)
  real(r8), intent(out), dimension(5) :: paccq
  ! parameters for cloud self collection (N)
  real(r8), intent(out), dimension(3) :: psccN
  ! parameters for rain self collection (N)
  real(r8), intent(out), dimension(3) :: pscrN
  real(r8), intent(in) :: pi
  real(r8) :: sxth,rhow,f_M32q
  integer, intent(in) :: iulog

  sxth  = 1./6.
  rhow   = 1000.

  f_M32q = pi*sxth*rhow

! set autoconversion mass parameters
! updated way of normalizing
if (iautoq.eq.1) then
   allocate(pautoq(3))
   pautoq(1) = 1D1**log_a_auto_t1
   pautoq(2) = b_auto_t1
   pautoq(3) = 2. - b_auto_t1
elseif (iautoq.eq.2) then
   allocate(pautoq(6))
   pautoq(1) = 1D1**log_a_auto_t1
   pautoq(2) = b_auto_t1
   pautoq(3) = 2. - b_auto_t1
   pautoq(4) = 1D1**log_a_auto_t2
   pautoq(5) = b_auto_t2_mc
   pautoq(6) = 2. - b_auto_t2_mc

elseif (iautoq.eq.3) then
   allocate(pautoq(8))
   pautoq(1) = 1D1**log_a_auto_t1
   pautoq(2) = b_auto_t1
   pautoq(3) = 2. - b_auto_t1
   pautoq(4) = 1D1**log_a_auto_t2
   pautoq(5) = b_auto_t2_mc
   pautoq(6) = 2. - b_auto_t2_mc - b_auto_t2_n
   pautoq(7) = b_auto_t2_mr
   pautoq(8) = b_auto_t2_n - b_auto_t2_mr
endif

   pautoN(1) = 1D1**log_mc_auto_inv
   pautoN(2) = 0.5*pautoN(1)

   paccq(1) = 1D1**log_a_acc
   paccq(2) = b_acc_mc
   paccq(3) = 1. - b_acc_mc
   paccq(4) = b_acc_mr
   paccq(5) = 1. - b_acc_mr

   psccN(1) = -1D1**log_a_sc_c
   psccN(2) = b_sc_c
   psccN(3) = 2. - b_sc_c

   pscrN(1) = 1D1**log_a_sc_r
   pscrN(2) = b_sc_r
   pscrN(3) = 2. - b_sc_r

end subroutine BOSS_init

subroutine BOSS_liq_autoconversion(microp_uniform,iautoq,relvar,qcic,ncic,qric,nric,pautoq,pautoN,prc,nprc,nprc1,vlen)

  use micro_pumas_utils, only: var_coef_r8_vect

  integer, intent(in) :: iautoq
  integer, intent(in) :: vlen
  logical, intent(in) :: microp_uniform
  real(r8), dimension(vlen), intent(in) :: relvar

  real(r8), dimension(vlen), intent(in) :: qcic
  real(r8), dimension(vlen), intent(in) :: ncic
  real(r8), dimension(vlen), intent(in) :: qric
  real(r8), dimension(vlen), intent(in) :: nric
  real(r8), dimension(vlen), intent(out) :: prc
  real(r8), dimension(vlen), intent(out) :: nprc
  real(r8), dimension(vlen), intent(out) :: nprc1

!  real(r8), intent(in), dimension(3) :: pautoq
  real(r8), intent(in), allocatable, dimension(:) :: pautoq
  real(r8), intent(in), dimension(2) :: pautoN

  real(r8), dimension(vlen) :: prc_coef, prc_coef2
  integer :: i

  !$acc data create (prc_coef)
  ! Take variance into account, or use uniform value.

  if (.not. microp_uniform) then
     call var_coef_r8_vect(relvar, pautoq(2), prc_coef,vlen)
     if ((iautoq.eq.2).or.(iautoq.eq.3)) then
        call var_coef_r8_vect(relvar, pautoq(4), prc_coef2,vlen)
     endif

  else

     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector
     do i = 1,vlen
        prc_coef(i) = 1._r8
        prc_coef2(i) = 1._r8
     end do
     !$acc end parallel
  end if

  do i=1,vlen
     if (qcic(i) >= icsmall) then
        if (iautoq.eq.1) then ! 1-term auto wo rain dependence
           prc(i) = prc_coef(i)*pautoq(1)*(qcic(i)**pautoq(2)*ncic(i)**pautoq(3))
        elseif (iautoq.eq.2) then ! 2-term auto wo rain dependence
           prc(i) = prc_coef(i)*pautoq(1)*(qcic(i)**pautoq(2)*ncic(i)**pautoq(3)) + &
                prc_coef2(i)*pautoq(4)*qcic(i)**pautoq(5)*ncic(i)**pautoq(6)
        elseif (iautoq.eq.3) then ! 2-term auto w rain dependence
           if ((qric(i).gt.0) .and. (nric(i).gt.0)) then
              prc(i) = prc_coef(i)*pautoq(1)*(qcic(i)**pautoq(2)*ncic(i)**pautoq(3)) + &
                   prc_coef2(i)*pautoq(4)*qcic(i)**pautoq(5)*ncic(i)**pautoq(6)*qric(i)**pautoq(7)*nric(i)**pautoq(8)
           else
              prc(i) = prc_coef(i)*pautoq(1)*(qcic(i)**pautoq(2)*ncic(i)**pautoq(3))
           endif
        endif

        ! handle overflow (from P3)
        if(isnan(prc(i)) .or. (prc(i).ge.1e38)) then
           prc = 1e38
        endif

        nprc(i) = pautoN(2)*prc(i)
        nprc1(i)= pautoN(1)*prc(i)

        ! handle overflow
        if(isnan(nprc(i)) .or. (nprc(i).ge.1e38)) then
           nprc(i) = 1e38
        endif
        if(isnan(nprc1(i)) .or. (nprc1(i).ge.1e38)) then
           nprc1(i) = 1e38
        endif

     else
        prc(i)   = 0._r8
        nprc(i)  = 0._r8
        nprc1(i) = 0._r8
     end if
  end do
  !$acc end parallel

end subroutine BOSS_liq_autoconversion

subroutine BOSS_accrete_cloud_water_rain(microp_uniform,relvar,qcic,ncic,qric,nric,paccq,pra,npra,vlen)

  use micro_pumas_utils, only: var_coef_r8_vect

  integer, intent(in) :: vlen
  logical, intent(in) :: microp_uniform
  real(r8), dimension(vlen), intent(in) :: relvar

  ! In-cloud rain                                 
  real(r8), dimension(vlen), intent(in) :: qric ! MMR
  real(r8), dimension(vlen), intent(in) :: nric ! Number
  ! Cloud droplets                          
  real(r8), dimension(vlen), intent(in) :: qcic ! MMR
  real(r8), dimension(vlen), intent(in) :: ncic ! Number
  ! Output tendencies                          
  real(r8), dimension(vlen), intent(out) :: pra  ! MMR
  real(r8), dimension(vlen), intent(out) :: npra ! Number

  real(r8), intent(in), dimension(5) :: paccq

  real(r8), dimension(vlen) :: pra_coef
  integer :: i                 



  !$acc data create (prc_coef)
  ! Take variance into account, or use uniform value.
  if (.not. microp_uniform) then                                                                      
     call var_coef_r8_vect(relvar, paccq(2), pra_coef,vlen)                                      
  else                                                                                                   
     !$acc parallel vector_length(VLENS) default(present)                                                
     !$acc loop gang vector                                                                              
     do i = 1,vlen                                                                                       
        pra_coef(i) = 1._r8                                                                              
     end do                                                                                              
     !$acc end parallel                                                                                  
  end if                   

  !$acc data create (pra_coef)                                                    
  !$acc parallel vector_length(VLENS) default(present)                            
  !$acc loop gang vector                                                          
  do i=1,vlen                                                                     
     if (qric(i) >= qsmall .and. qcic(i) >= qsmall) then                         
        ! trude test
!        pra_coef(i)=1._r8
        ! include sub-grid distribution of cloud water
        pra(i) = pra_coef(i)*paccq(1)*qcic(i)**paccq(2)*ncic(i)**paccq(3)*qric(i)**paccq(4)*nric(i)**paccq(5)
!        pra(i) = paccq(1)*qcic(i)**paccq(2)*ncic(i)**paccq(3)*qric(i)**paccq(4)*nric(i)**paccq(5)
        !    ncacc = qcacc*ncic(i,k)/qcic(i,k)
        npra(i) = pra(i)*ncic(i)/qcic(i)
     else                                    
        pra(i)  = 0._r8
        npra(i) = 0._r8                                                           
     end if                                                                       
  end do                                                                          
  !$acc end parallel
  !$acc end data       
end subroutine BOSS_accrete_cloud_water_rain

subroutine BOSS_self_collection_rain(qric, nric, pscrN, nragg, vlen)
  integer,                   intent(in) :: vlen
  real(r8), dimension(vlen), intent(in) :: qric ! MMR
  real(r8), dimension(vlen), intent(in) :: nric ! Number
  ! Output number tendency
  real(r8), dimension(vlen), intent(out) :: nragg
  real(r8), intent(in), dimension(3) :: pscrN
  integer :: i

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen                                                                  
     if (qric(i) >= qsmall) then
!        nragg(i) = -8._r8*nric(i)*qric(i)*rho(i)
        nragg(i) = -pscrN(1)*qric(i)**pscrN(2)*nric(i)**pscrN(3)
     else
        nragg(i) = 0._r8
     end if
  end do
  !$acc end parallel
end subroutine BOSS_self_collection_rain


subroutine BOSS_self_collection_cloud(microp_uniform, relvar, qcic, ncic, psccN, ncagg, vlen)

  use micro_pumas_utils, only: var_coef_r8_vect
  
  integer,                   intent(in) :: vlen
  logical, intent(in) :: microp_uniform
  real(r8), dimension(vlen), intent(in) :: relvar

  real(r8), dimension(vlen), intent(in) :: qcic ! MMR
  real(r8), dimension(vlen), intent(in) :: ncic ! Number
  ! Output number tendency
  real(r8), dimension(vlen), intent(out) :: ncagg
  real(r8), intent(in), dimension(3) :: psccN

  real(r8), dimension(vlen) :: ncagg_coef 
  integer :: i

  !$acc data create (prc_coef)                                                                           
                                                                                                         
  ! Take variance into account, or use uniform value.                                                    
  if (.not. microp_uniform) then                                                                         
     call var_coef_r8_vect(relvar, psccN(2), ncagg_coef,vlen)                                      
  else                                                                                                   
     !$acc parallel vector_length(VLENS) default(present)                                                
     !$acc loop gang vector                                                                              
     do i = 1,vlen                                                                                       
        ncagg_coef(i) = 1._r8                                                                              
     end do                                                                                              
     !$acc end parallel                                                                                  
  end if                   

  
  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen                                                                  
     if (qcic(i) >= qsmall) then
!        nragg(i) = -8._r8*nric(i)*qric(i)*rho(i)
        ! trude test
!        ncagg_coef(i)=1._r8
        ncagg(i) = ncagg_coef(i)*(-psccN(1)*qcic(i)**psccN(2)*ncic(i)**psccN(3))
!        ncagg(i) = (-psccN(1)*qcic(i)**psccN(2)*ncic(i)**psccN(3))

        ! handle overflow                                     
        if(isnan(ncagg(i)) .or. (ncagg(i).le.-1e38)) then
           ncagg(i) = -1e38
        endif

     else
        ncagg(i) = 0._r8
     end if
  end do
  !$acc end parallel
end subroutine BOSS_self_collection_cloud

end module BOSS_utils
