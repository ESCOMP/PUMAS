module pumas_stochastic_collect_tau
! From Morrison (Lebo, originally TAU bin code)
! Gettelman and Chen 2018
!the subroutines take in air density, air temperature, and the bin mass boundaries, and 
!output the mass and number mixing ratio tendencies in each bin directly.
!this is then wrapped for CAM. 

! note, this is now coded locally. Want the CAM interface to be over i,k I think.

#ifndef HAVE_GAMMA_INTRINSICS
use shr_spfn_mod, only: gamma => shr_spfn_gamma
#endif

use shr_kind_mod,      only: r8=>shr_kind_r8
use cam_history,       only: addfld
use micro_pumas_utils, only: pi, rhow, qsmall, VLENS
use cam_logfile,       only: iulog

implicit none
private
save

! Subroutines
public :: pumas_stochastic_kernel_init, pumas_stochastic_collect_tau_tend

!In the module top, declare the following so that these can be used throughout the module:

integer, parameter, public  :: ncd = 35
integer, parameter, public  :: ncdp = ncd + 1
integer, parameter, public  :: ncdl = ncd
integer, parameter, public  :: ncdpl = ncdl+1

! for Zach's collision-coalescence code

real(r8), private :: knn(ncd,ncd)

real(r8), private :: mmean(ncd), diammean(ncd)       ! kg & m at bin mid-points
real(r8), private :: medge(ncdp), diamedge(ncdp)     ! kg & m at bin edges 
integer, private  :: cutoff_id                       ! cutoff between cloud water and rain drop, D = 40 microns

! Assume 6 microns for each...
real(r8), parameter :: m1 = 4._r8/3._r8*pi*rhow*(6.e-6_r8)**3

!$acc declare create(knn,cutoff_id,mmean,diammean,medge,diamedge)

!===============================================================================
contains
!===============================================================================

subroutine calc_bins    

  real(r8) :: DIAM(ncdp)
  real(r8) :: X(ncdp)
  real(r8) :: radsl(ncdp)
  real(r8) :: radl(ncd)
  integer  :: L, lcl  
  real(r8) :: kkfac

!Then before doing any calculations you'll need to calculate the bin mass grid 
! (note this code could be cleaned up, I'm just taking it as it's used in our bin scheme). 
! This only needs to be done once, since we'll use the same bin mass grid for all calculations. 

! use mass doubling bins from Graham Feingold (note cgs units)

  DIAM(1)=1.5625*2.E-04_r8                ! cm
  X(1)=PI/6._r8*DIAM(1)**3*rhow/1000._r8  ! rhow kg/m3 --> g/cm3 
  radsl(1) = X(1)                         ! grams 

  DO l=2,ncdp
     X(l)=2._r8*X(l-1)
     DIAM(l)=(6._r8/pi*X(l)*1000._r8/rhow)**(1._r8/3._r8)  ! cm
     radsl(l)=X(l)             
  ENDDO

! now get bin mid-points

  do l=1,ncd
     radl(l)=(radsl(l)+radsl(l+1))/2._r8         ! grams   
     diammean(l) = (6._r8/pi*radl(l)*1000._r8/rhow)**(1._r8/3._r8) ! cm
  end do

! set bin grid for method of moments

  ! for method of moments

  do lcl = 1,ncd+1
     medge(lcl) = radsl(lcl)               ! grams
     diamedge(lcl) = DIAM(lcl)             ! cm
  enddo

  do lcl = 1,ncd
     mmean(lcl) = radl(lcl)  
     diammean(lcl) = diammean(lcl)
  enddo

  do lcl = ncdp,1,-1
     if( diamedge(lcl).ge.40.e-4_r8 ) cutoff_id = lcl
  end do  

end subroutine calc_bins

subroutine pumas_stochastic_kernel_init(kernel_filename)

  use cam_history_support, only: add_hist_coord

  character(len=*), intent(in) :: kernel_filename  ! Full pathname to kernel file

  integer :: iunit ! unit number of opened file for collection kernel code from a lookup table.

  integer :: idd, jdd
  real(r8) :: kkfac

  call calc_bins

! Read in the collection kernel code from a lookup table. Again, this only needs to be done once.
! use kernel from Zach (who got it from Jerry)

  KNN(:,:)=0._r8 ! initialize values
  kkfac=1.5_r8   ! from Zach

  open(newunit=iunit,file=kernel_filename,status='old')

  do idd=1,ncd
     do jdd=1,idd
        READ(iunit,941) KNN(IDD,JDD)
941     FORMAT(2X,E12.5)

        KNN(IDD,JDD)=(mmean(IDD)*kkfac+mmean(JDD)*kkfac)*KNN(IDD,JDD)
        if (knn(idd,jdd) < 0._r8) knn(idd,jdd)=0._r8
     end do
  end do

  !$acc update device(knn,cutoff_id,mmean,diammean,medge,diamedge)

end subroutine pumas_stochastic_kernel_init

!main driver routine
!needs to pull in i,k fields (so might need dimensions here too)

subroutine pumas_stochastic_collect_tau_tend(deltatin,t,rho,qc,qr,qcin,     &
                        ncin,qrin,nrin,lcldm,precip_frac,mu_c,lambda_c,     &
                        n0r,lambda_r,qcin_new,ncin_new,qrin_new,nrin_new,   &
                        qctend,nctend,qrtend,nrtend,qctend_TAU,nctend_TAU,  &
                        qrtend_TAU,nrtend_TAU,scale_qc,scale_nc,scale_qr,   &
                        scale_nr,amk_c,ank_c,amk_r,ank_r,amk,ank,amk_out,   &
                        ank_out,gmnnn_lmnnn_TAU,mgncol,nlev)

  !inputs: mgncol,nlev,t,rho,qcin,ncin,qrin,nrin
  !outputs: qctend,nctend,qrtend,nrtend
  !not sure if we want to output bins (extra dimension). Good for testing?  
  
  integer, intent(in) :: mgncol,nlev
  real(r8), intent(in) :: deltatin
  real(r8), intent(in) :: t(mgncol,nlev)
  real(r8), intent(in) :: rho(mgncol,nlev)
  real(r8), intent(in) :: qc(mgncol,nlev)
  real(r8), intent(in) :: qr(mgncol,nlev)
  real(r8), intent(in) :: qcin(mgncol,nlev)
  real(r8), intent(in) :: ncin(mgncol,nlev)
  real(r8), intent(in) :: qrin(mgncol,nlev)
  real(r8), intent(in) :: nrin(mgncol,nlev)
  real(r8), intent(in) :: lcldm(mgncol,nlev)
  real(r8), intent(in) :: precip_frac(mgncol,nlev)
  real(r8), intent(inout) :: qctend(mgncol,nlev)
  real(r8), intent(inout) :: nctend(mgncol,nlev)
  real(r8), intent(inout) :: qrtend(mgncol,nlev)
  real(r8), intent(inout) :: nrtend(mgncol,nlev)
  real(r8), intent(out) :: qctend_TAU(mgncol,nlev)
  real(r8), intent(out) :: nctend_TAU(mgncol,nlev)
  real(r8), intent(out) :: qrtend_TAU(mgncol,nlev)
  real(r8), intent(out) :: nrtend_TAU(mgncol,nlev)
  
  real(r8), intent(out) :: scale_qc(mgncol,nlev)
  real(r8), intent(out) :: scale_nc(mgncol,nlev)
  real(r8), intent(out) :: scale_qr(mgncol,nlev)
  real(r8), intent(out) :: scale_nr(mgncol,nlev)
  
  real(r8), intent(out) :: amk_c(mgncol,nlev,ncd)
  real(r8), intent(out) :: ank_c(mgncol,nlev,ncd)
  real(r8), intent(out) :: amk_r(mgncol,nlev,ncd)
  real(r8), intent(out) :: ank_r(mgncol,nlev,ncd)
  real(r8), intent(out) :: amk(mgncol,nlev,ncd)
  real(r8), intent(out) :: ank(mgncol,nlev,ncd)
  real(r8), intent(out) :: amk_out(mgncol,nlev,ncd)
  real(r8), intent(out) :: ank_out(mgncol,nlev,ncd)
  
  real(r8), intent(out) :: mu_c(mgncol,nlev)
  real(r8), intent(out) :: lambda_c(mgncol,nlev)
  real(r8), intent(out) :: lambda_r(mgncol,nlev)
  real(r8), intent(out) :: n0r(mgncol,nlev)
  
  real(r8), intent(out) :: qcin_new(mgncol,nlev)
  real(r8), intent(out) :: ncin_new(mgncol,nlev)
  real(r8), intent(out) :: qrin_new(mgncol,nlev)
  real(r8), intent(out) :: nrin_new(mgncol,nlev)
  real(r8), intent(out) :: gmnnn_lmnnn_TAU(mgncol,nlev)

  ! Local variables
  
  integer :: i,k,n,lcl
  integer :: cutoff_amk(mgncol,nlev),cutoff(mgncol,nlev)
  
  real(r8) :: all_gmnnn,all_lmnnn
  real(r8) :: qscl
  
  real(r8) :: qcin_old(mgncol,nlev)
  real(r8) :: ncin_old(mgncol,nlev)
  real(r8) :: qrin_old(mgncol,nlev)
  real(r8) :: nrin_old(mgncol,nlev)
  
  real(r8) :: amk0(mgncol,nlev,ncd)
  real(r8) :: ank0(mgncol,nlev,ncd)
  real(r8) :: gnnnn(mgncol,nlev,ncd)
  real(r8) :: gmnnn(mgncol,nlev,ncd)
  real(r8) :: lnnnn(mgncol,nlev,ncd)
  real(r8) :: lmnnn(mgncol,nlev,ncd)
  real(r8) :: gnnnn0(mgncol,nlev,ncd)
  real(r8) :: gmnnn0(mgncol,nlev,ncd)
  real(r8) :: lnnnn0(mgncol,nlev,ncd)
  real(r8) :: lmnnn0(mgncol,nlev,ncd)
  
  integer, parameter :: sub_step = 60

  !$acc data create (cutoff_amk,cutoff,qcin_old,ncin_old,qrin_old, &
  !$acc              nrin_old,amk0,ank0,gnnnn,gmnnn,lnnnn,lmnnn,   &
  !$acc              gnnnn0,gmnnn0,lnnnn0,lmnnn0)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)  
  do k=1,nlev
     do i=1,mgncol
        cutoff(i,k) = cutoff_id - 1
     end do
  end do
  !$acc end parallel

  ! First make bins from cam size distribution (bins are diagnostic)
  
  call cam_bin_distribute(qc,qr,qcin,ncin,qrin,nrin,mu_c,lambda_c, &
                          lambda_r,n0r,lcldm,precip_frac,scale_qc, &
                          scale_nc,scale_qr,scale_nr,amk_c,ank_c,  &
                          amk_r,ank_r,amk,ank,cutoff_amk,mgncol,nlev)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)  
  do k=1,nlev
     do i=1,mgncol
        if ( cutoff_amk(i,k) > 0 ) then
           cutoff(i,k) = cutoff_amk(i,k)
        end if
     end do
  end do
  !$acc end parallel
 
!Then call the subroutines that actually do the calculations. The inputs/ouputs are described in comments below. 

!This part of the code needs to be called each time for each process rate calculation 
! (i.e., for each sampled cloud/rain gamma distribution):

! note: variables passed to compute_column_params are all inputs,
! outputs from this subroutine are stored as global variables

! inputs: t --> input air temperature (K)
!         rho --> input air density (kg/m^3)
!         medge --> bin mass boundary (g) 
!         amk --> array of bin mass mixing ratio, i.e., the input drop mass distribution (kg/kg)
!         ank --> array of bin number mixing ratio, i.e., the input drop number distribution (kg^-1)

! inputs: medge --> bin mass boundary (g), same as above

! outputs: gnnnn --> bin number mass mixing tendency gain, array in bins (#/cm^3/s)
!          lnnnn --> bin number mass mixing tendency loss, array in bins (#/cm^3/s)
!          gmnnn --> bin mass mixing ratio tendency gain, array in bins (g/cm^3/s) 
!          lmnnn --> bin mass mixing ratio tendency loss, array in bins (g/cm^3/s)


! Call Kernel

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)  
  do k=1,nlev
     do i=1,mgncol
        qcin_new(i,k) = 0._r8
        ncin_new(i,k) = 0._r8
        qrin_new(i,k) = 0._r8
        nrin_new(i,k) = 0._r8
        
        qcin_old(i,k) = 0._r8
        ncin_old(i,k) = 0._r8
        qrin_old(i,k) = 0._r8
        nrin_old(i,k) = 0._r8
        
        qctend_TAU(i,k) = 0._r8
        nctend_TAU(i,k) = 0._r8
        qrtend_TAU(i,k) = 0._r8
        nrtend_TAU(i,k) = 0._r8
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(3)
  do lcl=1,ncd
     do k=1,nlev
        do i=1,mgncol
           gnnnn(i,k,lcl) = 0._r8
           gmnnn(i,k,lcl) = 0._r8
           lnnnn(i,k,lcl) = 0._r8
           lmnnn(i,k,lcl) = 0._r8
        end do
     end do
  end do
  !$acc end parallel

! update qc, nc, qr, nr

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)  
  do k=1,nlev
     do i=1,mgncol
        !$acc loop seq
        do lcl=1,ncd
           amk0(i,k,lcl) = amk(i,k,lcl)
           ank0(i,k,lcl) = ank(i,k,lcl)
        end do
        ! substep bin code
        !$acc loop seq
        do n=1,sub_step
           call compute_coll_params(rho(i,k),medge,amk0(i,k,1:ncd),ank0(i,k,1:ncd),gnnnn0(i,k,1:ncd),gmnnn0(i,k,1:ncd),lnnnn0(i,k,1:ncd),lmnnn0(i,k,1:ncd))

           all_gmnnn=0._r8
           all_lmnnn=0._r8
           !scaling gmnnn, lmnnn
           !$acc loop seq
           do lcl=1,ncd
              all_gmnnn = all_gmnnn+gmnnn0(i,k,lcl)
              all_lmnnn = all_lmnnn+lmnnn0(i,k,lcl)
           end do
 
           if ( (all_gmnnn == 0._r8) .or. (all_lmnnn == 0._r8) ) then
              !$acc loop seq
              do lcl=1,ncd
                 gmnnn0(i,k,lcl) = 0._r8
                 lmnnn0(i,k,lcl) = 0._r8
              end do
           else
              !$acc loop seq
              do lcl=1,ncd
                 lmnnn0(i,k,lcl) = lmnnn0(i,k,lcl)*(all_gmnnn/all_lmnnn)
              end do
           end if

           !$acc loop seq
           do lcl=1,ncd
              amk0(i,k,lcl) = amk0(i,k,lcl)+(gmnnn0(i,k,lcl)-lmnnn0(i,k,lcl))*1.e3_r8/ &
                          rho(i,k)*deltatin/dble(sub_step)
              ank0(i,k,lcl) = ank0(i,k,lcl)+(gnnnn0(i,k,lcl)-lnnnn0(i,k,lcl))*1.e6_r8/ &
                          rho(i,k)*deltatin/dble(sub_step)
              gmnnn(i,k,lcl) = gmnnn(i,k,lcl)+gmnnn0(i,k,lcl)/sub_step
              gnnnn(i,k,lcl) = gnnnn(i,k,lcl)+gnnnn0(i,k,lcl)/sub_step
              lmnnn(i,k,lcl) = lmnnn(i,k,lcl)+lmnnn0(i,k,lcl)/sub_step
              lnnnn(i,k,lcl) = lnnnn(i,k,lcl)+lnnnn0(i,k,lcl)/sub_step
           end do
        end do ! end of loop "sub_step"

        ! cloud water
        !$acc loop seq
        do lcl=1,cutoff(i,k)
           qcin_old(i,k) = qcin_old(i,k)+amk(i,k,lcl)
           ncin_old(i,k) = ncin_old(i,k)+ank(i,k,lcl)
           qcin_new(i,k) = qcin_new(i,k)+(gmnnn(i,k,lcl)-lmnnn(i,k,lcl))*1.e3_r8/rho(i,k)*deltatin
           ncin_new(i,k) = ncin_new(i,k)+(gnnnn(i,k,lcl)-lnnnn(i,k,lcl))*1.e6_r8/rho(i,k)*deltatin
           qctend_TAU(i,k) = qctend_TAU(i,k)+(amk0(i,k,lcl)-amk(i,k,lcl))/deltatin
           nctend_TAU(i,k) = nctend_TAU(i,k)+(ank0(i,k,lcl)-ank(i,k,lcl))/deltatin
           gmnnn_lmnnn_TAU(i,k) = gmnnn_lmnnn_TAU(i,k)+gmnnn(i,k,lcl)-lmnnn(i,k,lcl)
        end do

        ! rain
        !$acc loop seq
        do lcl=cutoff(i,k)+1,ncd
           qrin_old(i,k) = qrin_old(i,k)+amk(i,k,lcl)
           nrin_old(i,k) = nrin_old(i,k)+ank(i,k,lcl)
           qrin_new(i,k) = qrin_new(i,k)+(gmnnn(i,k,lcl)-lmnnn(i,k,lcl))*1.e3_r8/rho(i,k)*deltatin
           nrin_new(i,k) = nrin_new(i,k)+(gnnnn(i,k,lcl)-lnnnn(i,k,lcl))*1.e6_r8/rho(i,k)*deltatin
           qrtend_TAU(i,k) = qrtend_TAU(i,k)+(amk0(i,k,lcl)-amk(i,k,lcl))/deltatin
           nrtend_TAU(i,k) = nrtend_TAU(i,k)+(ank0(i,k,lcl)-ank(i,k,lcl))/deltatin
           gmnnn_lmnnn_TAU(i,k) = gmnnn_lmnnn_TAU(i,k)+gmnnn(i,k,lcl)-lmnnn(i,k,lcl)
        end do

        !$acc loop seq     
        do lcl=1,ncd
           amk_out(i,k,lcl) = amk(i,k,lcl) + (gmnnn(i,k,lcl)-lmnnn(i,k,lcl))*1.e3_r8/rho(i,k)*deltatin
           ank_out(i,k,lcl) = ank(i,k,lcl) + (gnnnn(i,k,lcl)-lnnnn(i,k,lcl))*1.e6_r8/rho(i,k)*deltatin
        end do
     
        qcin_new(i,k) = qcin_new(i,k)+qcin_old(i,k)
        ncin_new(i,k) = ncin_new(i,k)+ncin_old(i,k)
        qrin_new(i,k) = qrin_new(i,k)+qrin_old(i,k)
        nrin_new(i,k) = nrin_new(i,k)+nrin_old(i,k)
     end do
  end do
  !$acc end parallel

! Conservation checks 
! AG: Added May 2023

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)  
  do k=1,nlev
     do i=1,mgncol

        ! First make sure all not negative
        qcin_new(i,k)=max(qcin_new(i,k),0._r8)
        ncin_new(i,k)=max(ncin_new(i,k),0._r8)
        qrin_new(i,k)=max(qrin_new(i,k),0._r8)
        nrin_new(i,k)=max(nrin_new(i,k),0._r8)

! Now adjust so that sign is correct. qc_new,nc_new <= input, qr_new >= input
! NOte that due to self collection nr can be larger or smaller than input....
! Makes above check redundant I think.

        qcin_new(i,k)=min(qcin_new(i,k),qcin(i,k))
        ncin_new(i,k)=min(ncin_new(i,k),ncin(i,k))
        qrin_new(i,k)=max(qrin_new(i,k),qrin(i,k))

! Next scale mass...so output qc+qr is the same as input

        if ( (qcin_new(i,k)+qrin_new(i,k)) > 0._r8 ) then
           qscl = (qcin(i,k)+qrin(i,k))/(qcin_new(i,k)+qrin_new(i,k))
        else
           qscl = 0._r8
        end if
        qcin_new(i,k) = qcin_new(i,k) * qscl
        qrin_new(i,k) = qrin_new(i,k) * qscl

! Now zero nr,nc if either small or no mass?

        if ( qcin_new(i,k) < qsmall ) then
           ncin_new(i,k) = 0._r8
        end if
        
        if ( qrin_new(i,k) < qsmall ) then
           nrin_new(i,k) = 0._r8
        end if

!Finally add number if mass but no (or super small) number

        if ( qcin_new(i,k) > qsmall .and. ncin_new(i,k) < qsmall ) then
           ncin_new(i,k) = qcin_new(i,k)/m1
        end if
        
        if ( qrin_new(i,k) > qsmall .and. nrin_new(i,k) < qsmall) then
           nrin_new(i,k) = qrin_new(i,k)/m1
        end if

! Then recalculate tendencies based on difference
! Clip tendencies for cloud (qc,nc) to be <= 0. 
! Qrtend is not used in pumas (-qctend is used) but clip that too). 
! Nr tend can be muliply signed. 

        qctend_TAU(i,k)= min((qcin_new(i,k) - qcin(i,k)) / deltatin,0._r8)
        nctend_TAU(i,k)= min((ncin_new(i,k) - ncin(i,k)) / deltatin,0._r8)
        qrtend_TAU(i,k)= max((qrin_new(i,k) - qrin(i,k)) / deltatin,0._r8)
        nrtend_TAU(i,k)= (nrin_new(i,k) - nrin(i,k)) / deltatin

     end do
  end do
  !$acc end parallel

  !$acc end data

end subroutine pumas_stochastic_collect_tau_tend

subroutine cam_bin_distribute(qc_all,qr_all,qc,nc,qr,nr,mu_c,lambda_c, &
                              lambda_r,n0r,lcldm,precip_frac,scale_qc, &
                              scale_nc,scale_qr,scale_nr,amk_c,ank_c,  &
                              amk_r,ank_r,amk,ank,cutoff_amk,mgncol,nlev)

  implicit none

  integer, intent(in) :: mgncol,nlev
  real(r8), dimension(mgncol,nlev), intent(in) :: qc_all,qr_all,qc,nc,qr,nr,mu_c, &
                                                  lambda_c,lambda_r,n0r,lcldm,    &
                                                  precip_frac
  real(r8), dimension(mgncol,nlev,ncd), intent(out) :: amk_c,ank_c,amk_r,ank_r,amk,ank
  real(r8), dimension(mgncol,nlev), intent(out) :: scale_nc,scale_qc,scale_nr,scale_qr 
  integer, dimension(mgncol,nlev), intent(out) :: cutoff_amk

  ! Local variables

  integer  :: i,j,k
  real(r8) :: phi
  integer  :: id_max_qc, id_max_qr
  real(r8) :: max_qc, max_qr, min_amk 

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(3)
  do j=1,ncd
     do k=1,nlev
        do i=1,mgncol
           ank_c(i,k,j) = 0._r8
           amk_c(i,k,j) = 0._r8
           ank_r(i,k,j) = 0._r8
           amk_r(i,k,j) = 0._r8
           ank(i,k,j) = 0._r8
           amk(i,k,j) = 0._r8
        end do
     end do
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k=1,nlev
     do i=1,mgncol
        scale_nc(i,k) = 0._r8
        scale_qc(i,k) = 0._r8
        scale_nr(i,k) = 0._r8
        scale_qr(i,k) = 0._r8
        cutoff_amk(i,k) = 0

        id_max_qc = 0
        id_max_qr = 0
        max_qc = 0._r8
        max_qr = 0._r8

        ! cloud water, nc in #/m3 --> #/cm3
        if ( (qc_all(i,k) > qsmall) .and. (qc(i,k) > qsmall) ) then
           !$acc loop seq
           do j=1,ncd
              phi = nc(i,k)*lambda_c(i,k)**(mu_c(i,k)+1._r8)/ &
                    gamma(mu_c(i,k)+1._r8)*(diammean(j)*1.e-2_r8)**mu_c(i,k)* &
                    exp(-lambda_c(i,k)*diammean(j)*1.e-2_r8)  ! D cm --> m
              ank_c(i,k,j) = phi*(diamedge(j+1)-diamedge(j))*1.e-2_r8   ! D cm --> m
              amk_c(i,k,j) = phi*(diamedge(j+1)-diamedge(j))*1.e-2_r8*mmean(j)*1.e-3_r8  ! mass in bin g --> kg
              scale_nc(i,k) = scale_nc(i,k)+ank_c(i,k,j)
              scale_qc(i,k) = scale_qc(i,k)+amk_c(i,k,j) 
           end do
           scale_nc(i,k) = scale_nc(i,k)/nc(i,k)
           scale_qc(i,k) = scale_qc(i,k)/qc(i,k)

           !$acc loop seq
           do j=1,ncd
              ank_c(i,k,j) = ank_c(i,k,j)/scale_nc(i,k)*lcldm(i,k)
              amk_c(i,k,j) = amk_c(i,k,j)/scale_qc(i,k)*lcldm(i,k)
              if ( amk_c(i,k,j) > max_qc ) then
                 id_max_qc = j
                 max_qc = amk_c(i,k,j)
              end if
           end do
        end if

        ! rain drop
        if ( (qr_all(i,k) > qsmall) .and. (qr(i,k) > qsmall) ) then
           !$acc loop seq
           do j=1,ncd
              phi = n0r(i,k)*exp(-lambda_r(i,k)*diammean(j)*1.e-2_r8)   ! D cm --> m
              ank_r(i,k,j) = phi*(diamedge(j+1)-diamedge(j))*1.e-2_r8   ! D cm --> m  
              amk_r(i,k,j) = phi*(diamedge(j+1)-diamedge(j))*1.e-2_r8*mmean(j)*1.e-3_r8
              scale_nr(i,k) = scale_nr(i,k) + ank_r(i,k,j)
              scale_qr(i,k) = scale_qr(i,k) + amk_r(i,k,j)
           end do
           scale_nr(i,k) = scale_nr(i,k)/nr(i,k)
           scale_qr(i,k) = scale_qr(i,k)/qr(i,k)

           !$acc loop seq
           do j=1,ncd
              ank_r(i,k,j) = ank_r(i,k,j)/scale_nr(i,k)*precip_frac(i,k)
              amk_r(i,k,j) = amk_r(i,k,j)/scale_qr(i,k)*precip_frac(i,k)
              if ( amk_r(i,k,j) > max_qr ) then
                 id_max_qr = j
                 max_qr = amk_r(i,k,j)
              end if
           end do
        end if

        !$acc loop seq
        do j=1,ncd
           amk(i,k,j) = amk_c(i,k,j) + amk_r(i,k,j)
           ank(i,k,j) = ank_c(i,k,j) + ank_r(i,k,j)
        end do

        if ( (id_max_qc > 0) .and. (id_max_qr > 0) ) then
           if ( (max_qc/max_qr < 10._r8) .or. (max_qc/max_qr > 0.1_r8) ) then
              min_amk = amk(i,k,id_max_qc)
              !$acc loop seq
              do j=id_max_qc,id_max_qr
                 if ( amk(i,k,j) <= min_amk ) then
                    cutoff_amk(i,k) = j
                    min_amk = amk(i,k,j)
                 end if
              end do
           end if
        end if
     end do  ! end of loop "mgncol"
  end do     ! end of loop "nlev"
  !$acc end parallel

!input: qc,nc,qr,nr, medge (bin edges). May also need # bins?
!output: amk, ank (mixing ratio and number in each bin)

!this part will take a bit of thinking about.
!use size distribution parameters (mu, lambda) to generate the values at discrete size points
!need to also ensure mass conservation  

end subroutine cam_bin_distribute

! here are the subroutines called above that actually do the collision-coalescence calculations:

! The Kernel is from Jerry from many moons ago (included)

! I read in the file data and multiply by the summed mass of the individual bins 
! (with a factor of 1.5 so that the values represent the middle of the bin

! 941 FORMAT(2X,E12.5)
!     READ(iunit,941) KNN(IDD,JDD)
!     KNN(IDD,JDD)=(XK_GR(IDD)*kkfac+XK_GR(JDD)*kkfac)*KNN(IDD,JDD)

!where idd and jdd are the indexes for the bins and xk_gr is the mass of drops in a bin in grams
!

!************************************************************************************
! Setup variables needed for collection
! Either pass in or define globally the following variables
! tbase(height) - temperature in K as a function of height
! rhon(height) - air density as a function of height in kg/m^3
! xk_gr(bins) - mass of single drop in each bin in grams
! lsmall - small number
! QC - mass mixing ratio in kg/kg
! QN - number mixing ratio in #/kg
! All parameters are defined to be global in my version so that they are readily available throughout the code:
! SMN0,SNN0,SMCN,APN,AMN2,AMN3,PSIN,FN,FPSIN,XPSIN,HPSIN,FN2,XXPSIN (all arrays of drop bins)
!************************************************************************************

!AG: Global arrays need to be passed around I think? Right now at the module level. Is that okay?

SUBROUTINE COMPUTE_COLL_PARAMS(rhon,xk_gr,qc,qn,gnnnn,gmnnn,lnnnn,lmnnn)

  !$acc routine seq

  IMPLICIT NONE

! variable declarations (added by hm, 020118)
! note: vertical array looping is stripped out, this subroutine operates
! only on LOCAL values

  real(r8), dimension(ncd) :: qc,qn
  real(r8), dimension(ncdp) :: xk_gr
  real(r8) :: tbase,rhon
  integer :: lk
  integer :: l
  real(r8), parameter :: lsmall = 1.e-12_r8
  real(r8), dimension(ncd) :: smn0,snn0,smcn,amn2,amn3,psin,fn,fpsin, &
                               xpsin,hpsin,fn2,xxpsin
  real(r8) :: apn

  real(r8), dimension(ncd) :: gnnnn,gmnnn,lnnnn,lmnnn
  integer :: lm1,ll

  lk=ncd

  DO L=1,LK
     SMN0(L)=QC(L)*RHON/1.E3_r8
     SNN0(L)=QN(L)*RHON/1.E6_r8

     IF(SMN0(L).LT.lsmall.OR.SNN0(L).LT.lsmall)THEN
        SMN0(L)=0.0_r8
        SNN0(L)=0.0_r8
     ENDIF
  ENDDO

  DO L=1,LK
     IF(SMN0(L) .gt. 0._r8.AND.SNN0(L) .gt. 0._r8)THEN
        SMCN(L)=SMN0(L)/SNN0(L)
        IF((SMCN(L) .GT. 2._r8*XK_GR(L)))THEN
           SMCN(L) = (2._r8*XK_GR(L))
        ENDIF
        IF((SMCN(L) .LT. XK_GR(L)))THEN
           SMCN(L) = XK_GR(L)
        ENDIF
     ELSE
        SMCN(L)=0._r8
     ENDIF
     IF (SMCN(L).LT.XK_GR(L).OR.SMCN(L).GT.(2._r8*XK_GR(L)).OR.SMCN(L).EQ.0.0_r8)THEN
        APN=1.0_r8
     ELSE
        APN=0.5_r8*(1._r8+3._r8*(XK_GR(L)/SMCN(L))-2._r8*((XK_GR(L)/SMCN(L))**2._r8))
     ENDIF

     IF(SNN0(L) .GT. LSMALL)THEN
        AMN2(L)=APN*SMN0(L)*SMN0(L)/SNN0(L)
        AMN3(L)=APN*APN*APN*SMN0(L)*SMN0(L)*SMN0(L)/(SNN0(L)*SNN0(L))
     ELSE
        AMN2(L)=0._r8
        AMN3(L)=0._r8
     ENDIF

     IF(SMCN(L).LT.XK_GR(L))THEN
        PSIN(L)=0.0_r8
        FN(L)=2._r8*SNN0(L)/XK_GR(L)
     ELSE
        IF(SMCN(L).GT.(2._r8*XK_GR(L)))THEN
           FN(L)=0.0_r8
           PSIN(L)=2._r8*SNN0(L)/XK_GR(L)
        ELSE
           PSIN(L)=2._r8/XK_GR(L)*(SMN0(L)/XK_GR(L)-SNN0(L))
           FN(L)=2._r8/XK_GR(L)*(2._r8*SNN0(L)-SMN0(L)/XK_GR(L))
        ENDIF
     ENDIF

     IF(SNN0(L).LT.LSMALL.OR.SMN0(L).LT.LSMALL)THEN
        PSIN(L)=0.0_r8
        FN(L)=0.0_r8
     ENDIF

     FPSIN(L)=0.5_r8/XK_GR(L)*(PSIN(L)-FN(L))
     XPSIN(L)=2._r8*XK_GR(L)*PSIN(L)
     HPSIN(L)=PSIN(L)-0.5_r8*FN(L)
     FN2(L)=FN(L)/2._r8

     IF(L.GT.1)THEN
        XXPSIN(L)=XK_GR(L)*PSIN(L-1)
     ENDIF
  ENDDO

!************************************************************************************
! Compute collision coalescence
! Either pass in or define globally the following variables
! Gain terms begin with G, loss terms begin with L
! Second letter defines mass (M) or number (N)
! Third and fourth letters define the types of particles colling, i.e., NN means drops colliding with drops
! Last letter defines the category the new particles go into, in this case just N for liquid drops
! The resulting rates are in units of #/cm^3/s and g/cm^3/s
! Relies on predefined kernel array KNN(bins,bins) - see top of this file
!************************************************************************************

   GMNNN = 0._r8
   GNNNN = 0._r8
   LMNNN = 0._r8
   LNNNN = 0._r8
! remove verical array index, calculate gain/loss terms locally

  DO L=3,LK-1
     LM1=L-1
     DO LL=1,L-2
        GNNNN(L)=GNNNN(L)+(PSIN(LM1)*SMN0(LL)-FPSIN(LM1)*AMN2(LL))*KNN(LM1,LL)
        GMNNN(L)=GMNNN(L)+(XK_GR(L)*PSIN(LM1)*SMN0(LL)+FN2(LM1)*AMN2(LL)-FPSIN(LM1)*AMN3(LL))*KNN(LM1,LL)
     ENDDO
  ENDDO

  DO L=2,LK-1
     LM1=L-1
     GNNNN(L)=GNNNN(L)+0.5_r8*SNN0(LM1)*SNN0(LM1)*KNN(LM1,LM1)
     GMNNN(L)=GMNNN(L)+0.5_r8*(SNN0(LM1)*SMN0(LM1)+SMN0(LM1)*SNN0(LM1))*KNN(LM1,LM1)
     DO LL=1,L-1
        LNNNN(L)=LNNNN(L)+(PSIN(L)*SMN0(LL)-FPSIN(L)*AMN2(LL))*KNN(L,LL)
        GMNNN(L)=GMNNN(L)+(SMN0(LL)*SNN0(L)-PSIN(L)*AMN2(LL)+FPSIN(L)*AMN3(LL))*KNN(L,LL)
        LMNNN(L)=LMNNN(L)+(XPSIN(L)*SMN0(LL)-HPSIN(L)*AMN2(LL))*KNN(L,LL)
     ENDDO
  ENDDO

  DO L=1,LK-1
     DO LL=L,LK-1
        LNNNN(L)=LNNNN(L)+(SNN0(LL)*SNN0(L))*KNN(LL,L)
        LMNNN(L)=LMNNN(L)+(SNN0(LL)*SMN0(L))*KNN(LL,L)
     ENDDO
  ENDDO

END SUBROUTINE COMPUTE_COLL_PARAMS


end module pumas_stochastic_collect_tau


