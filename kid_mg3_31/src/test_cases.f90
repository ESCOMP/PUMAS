! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Some test cases
!
!

module test_cases

  Use typeKind
  Use parameters, only : dt, num_h_moments, num_h_bins &
     , aero_N_init, aero_rd_init, aero_sig_init
  Use runtime, only : n_times
  Use common_physics
  Use column_variables
  Use input_variables
  Use derived_fields, only : calc_derived_fields
  Use physconst, only : Ru=>R, r_on_cp, g, p0, pi, cp
  Use interpolation, only : interpolate, smooth1D
  Use switches
  Use class_species, only: species_allocate

  Use aerosols, only: moment_logn, awp=>wp, set_aerosol

  Implicit none

  !local variables
  integer, private :: k, itime, ih, j
  
  ! Levels for interpolation
  real(wp), allocatable :: &
        pHeight(:)         & ! height
       ,pTheta(:)          & ! theta
       ,pqv(:)             & ! qv
       ,pRH(:)               ! RH

  ! 1-D column for interpolation. This is assigned to 
  ! 2-D arrays 
  real(wp), allocatable :: &
        Theta_1d(:)          & ! theta
       ,qv_1d(:)               ! qv  

  real(wp), allocatable :: RH(:) ! RH

  ! Local variables/parameters for GCSS warm rain case 1
  real(wp), private :: A, B, C, utop
  ! Local variables/parameters for GCSS  deep case 1
  real(wp), private :: Dw, z0w, Dq, z0q, zinv, zbase
  ! Local variables/parameters for GCSS  deep case 2
  real(wp), private :: w0, w1,  tup, zup, s, sigma, eta
  ! Local variables/parameters for GCSS 2D_Cu case
  real(wp), private :: hx, x0, hz, z0, alpha, beta & 
       , H2, AMP1, AMP2, a10, a1, a2, xc, f1, f2   &
       , t1

  
  real(wp) :: &
        maxT  & ! max time
       ,maxZ  & ! max height
       ,maxX    ! max horizontal extent


  
contains

  subroutine standard_cases(icase)
    
    integer, intent(in) :: icase

    !local variables
    real(wp) :: t     ! temporary local time variable

    ! for setting aerosol
    integer, parameter :: Ninit=naerosol ! Number of aerosol we can intialize
    integer :: indices(Ninit)
    logical :: lainits(Ninit) = .false.
    real(wp) :: Nds(Ninit)
    real(wp) :: sigmas(Ninit)
    real(wp) :: rds(Ninit)
    real(wp) :: densitys(Ninit)
    real(wp) :: fscale(nz)

    select case (icase)

    case(ipassive_sed_bulk1)
       !==============================================
       ! GCSS microphysics intercomparison Warm Rain 1
       !==============================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1:3)=(/ 3000., 2500., 2800. /)
       if (all(wctrl==0.))wctrl(1)=2.
       if (all(tctrl==0.))tctrl(1)=3600.
       if (all(pctrl_v==0.))pctrl_v(1:3)=(/ 2., .001, 100.e6 /)
       if (ipctrl==0)ipctrl=1

       maxZ=zctrl(1)
       maxT=tctrl(1)
       n_times=int(maxT/dt)

       call set_standard_profile(ipctrl,maxZ)

       call allocate_forcing(nz,nx,n_times)

       !initialize some rain (species 2?)
       ih=int(pctrl_v(1))
       l_hinit(ih)=.true.
       do k=1,nz
          call species_allocate(hydrometeors_init(k,ih) &
               , num_h_moments(ih), num_h_bins(ih), ih)
       end do
       do k=1,nz
          if (z(k)>zctrl(2) .and. z(k)<zctrl(3))then
             hydrometeors_init(k, ih)%moments(1,1)=pctrl_v(2)
             if (num_h_moments(ih)>1)then
                hydrometeors_init(k, ih)%moments(1,2)=pctrl_v(3)
             end if
          end if
       end do
    
    case(igcss_warm1)
       !==============================================
       ! GCSS microphysics intercomparison Warm Rain 1
       !==============================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1)=3000.
       if (all(wctrl==0.))wctrl(1)=2.
       if (all(tctrl==0.))tctrl(1:2)=(/3600., 600./)
       if (ipctrl==0)ipctrl=1

       maxZ=zctrl(1)
       maxT=tctrl(1)
       n_times=int(maxT/dt)

       call set_standard_profile(ipctrl,maxZ)

       call allocate_forcing(nz,nx,n_times)

       do itime=1,n_times
          t=itime*dt
          time_in(itime)=t
          if (t<tctrl(2))then
             w_t(:,:,itime)= wctrl(1)*sin(pi*t/tctrl(2))
          else
             w_t(:,:,itime)=0
          end if
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
       end do

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale)
 
    case(igcss_warm1A)
       !==============================================
       ! GCSS microphysics intercomparison Warm Rain 1
       !==============================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1)=3000.
       if (all(wctrl==0.))wctrl(1)=0.3
       if (all(tctrl==0.))tctrl(1:2)=(/3600., 600./)
       if (ipctrl==0)ipctrl=1

       maxZ=zctrl(1)
       maxT=tctrl(1)
       n_times=int(maxT/dt)

       call set_standard_profile(ipctrl,maxZ)

       call allocate_forcing(nz,nx,n_times)

       do itime=1,n_times
          t=itime*dt
          time_in(itime)=t
          if (t<tctrl(2))then
             w_t(:,:,itime)= wctrl(1)*sin(pi*t/tctrl(2))
          else
             w_t(:,:,itime)=0
          end if
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
       end do

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale)

    case(igcss_warmX)
       !==============================================
       ! GCSS microphysics intercomparison Warm Rain 1: slow updraft
       !==============================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1)=3000.
       if (all(wctrl==0.))wctrl(1)=0.3
       if (all(tctrl==0.))tctrl(1:2)=(/9600., 600./)
       if (ipctrl==0)ipctrl=8  !use different profile...

       maxZ=zctrl(1)
       maxT=tctrl(1)
       n_times=int(maxT/dt)

       call set_standard_profile(ipctrl,maxZ)

       call allocate_forcing(nz,nx,n_times)

       do itime=1,n_times
          t=itime*dt
          time_in(itime)=t
          if (t<tctrl(2))then
             w_t(:,:,itime)= wctrl(1)*sin(pi*t/tctrl(2))
          else
             w_t(:,:,itime)=0
          end if
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
       end do

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 
    case(igcss_warm0)
       !==============================================
       ! GCSS microphysics intercomparison Warm Rain 0: NO updraft
       !==============================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1)=3000.
       if (all(wctrl==0.))wctrl(1)=0.0
       if (all(tctrl==0.))tctrl(1:2)=(/21600., 600./)
       if (ipctrl==0)ipctrl=1  !use different profile...

       maxZ=zctrl(1)
       maxT=tctrl(1)
       n_times=int(maxT/dt)

       call set_standard_profile(ipctrl,maxZ)

       call allocate_forcing(nz,nx,n_times)

       do itime=1,n_times
          t=itime*dt
          time_in(itime)=t
          w_t(:,:,itime)=0
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
       end do

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale)

    case(igcss_warm00)
       !==============================================
       ! GCSS microphysics intercomparison Warm Rain 0: NO updraft
       !==============================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1)=3000.
       if (all(wctrl==0.))wctrl(1)=0.0
       if (all(tctrl==0.))tctrl(1:2)=(/28800., 600./)
       if (ipctrl==0)ipctrl=1  !use different profile...

       maxZ=zctrl(1)
       maxT=tctrl(1)
       n_times=int(maxT/dt)

       call set_standard_profile(ipctrl,maxZ)

       call allocate_forcing(nz,nx,n_times)

       do itime=1,n_times
          t=itime*dt
          time_in(itime)=t
          w_t(:,:,itime)=0
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
       end do

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale)

    case(igcss_warm2)
       !==============================================
       ! GCSS microphysics intercomparison Warm Rain 2
       !==============================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1)=3000.
       if (all(wctrl==0.))wctrl(1)=2.
       if (all(tctrl==0.))tctrl(1:2)=(/7200., 600./)
       if (ipctrl==0)ipctrl=1

       maxZ=zctrl(1)
       maxT=tctrl(1)
       n_times=int(maxT/dt)

       call set_standard_profile(ipctrl,maxZ)

       call allocate_forcing(nz,nx,n_times)

       do itime=1,n_times
          t=itime*dt
          time_in(itime)=t
          w_t(:,:,itime)=wctrl(1)*sin(pi*t/tctrl(2))
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
        end do

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 
    case(igcss_warm3)
       !==============================================
       ! GCSS microphysics intercomparison Warm Rain 3
       !==============================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1)=3000.
       if (all(wctrl==0.))wctrl(1)=2.
       if (all(tctrl==0.))tctrl(1:3)=(/3600., 600., 1200./)
       if (ipctrl==0)ipctrl=1

       maxZ=zctrl(1)
       maxT=tctrl(1)
       n_times=int(maxT/dt)

       call set_standard_profile(ipctrl,maxZ)

       call allocate_forcing(nz,nx,n_times)

       do itime=1,n_times
          t=itime*dt
          time_in(itime)=t
          w_t(:,:,itime)=wctrl(1)*sin(pi*t/tctrl(2))*exp(-t/tctrl(3))
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
       end do       


       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 
!    case(igcss_warm4)
!       !==============================================
!       ! GCSS microphysics intercomparison Warm Rain 4
!       !==============================================
!       !removed

    case(igcss_warm5)
       !==============================================
       ! GCSS microphysics intercomparison Warm Rain 5
       !==============================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1:6)=(/ 7000., 650., 1000., 2000., 3000.,&
            & 4000. /)
       if (all(wctrl==0.))wctrl(1)=7.
       if (all(tctrl==0.))tctrl(1:4)=(/5400., 1800., 3240., 3600./)
       if (ipctrl==0)ipctrl=3

       maxZ=zctrl(1)
       maxT=tctrl(1)
       n_times=int(maxT/dt)

       call set_standard_profile(ipctrl,maxZ)

       call allocate_forcing(nz,nx,n_times)

       do itime=1,n_times
          t=itime*dt
          time_in(itime)=t
          if (t < tctrl(4))then
             A=wctrl(1)*(1-((t-tctrl(2))/tctrl(2))**4)
             if (t <  tctrl(2)) then
                B=zctrl(2)+(zctrl(4)-zctrl(2))*t/tctrl(2)
                utop=zctrl(3)+(zctrl(6)-zctrl(3))*t/tctrl(2)
             else
                B=zctrl(4)
                if (t<tctrl(3))then
                   utop=zctrl(6)
                else
                   utop=zctrl(6)+(zctrl(5)-zctrl(6))*(t-tctrl(3))/(tctrl(4)-tctrl(3))
             end if
             end if
          else
             A=0
             B=1
             utop=1
          end if
          C=max(utop-b,b)
          do k=1, nz
             w_t(k,:,itime)=A*exp(-((z(k)-B)/C)**2)
          end do
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
       end do

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 
    case(igcss_warm6)
       !==============================================
       ! GCSS microphysics intercomparison Warm Rain 6
       !==============================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1:6)=(/ 3000., 300., 500., 1000., 1500.,&
            & 2000. /)
       if (all(wctrl==0.))wctrl(1)=3.
       if (all(tctrl==0.))tctrl(1:4)=(/4500., 1800., 3240., 3600./)
       if (ipctrl==0)ipctrl=2

       maxZ=zctrl(1)
       maxT=tctrl(1)
       n_times=int(maxT/dt)

       call set_standard_profile(ipctrl,maxZ)

       call allocate_forcing(nz,nx,n_times)

       do itime=1,n_times
          t=itime*dt
          time_in(itime)=t
          if (t < tctrl(4))then
             A=wctrl(1)*(1-((t-tctrl(2))/tctrl(2))**4)
             if (t <  tctrl(2)) then
                B=zctrl(2)+(zctrl(4)-zctrl(2))*t/tctrl(2)
                utop=zctrl(3)+(zctrl(6)-zctrl(3))*t/tctrl(2)
             else
                B=zctrl(4)
                if (t<tctrl(3))then
                   utop=zctrl(6)
                else
                   utop=zctrl(6)+(zctrl(5)-zctrl(6))*(t-tctrl(3))/(tctrl(4)-tctrl(3))
             end if
             end if
          else
             A=0
             B=1
             utop=1
          end if
          C=max(utop-b,b)
          do k=1, nz
             w_t(k,:,itime)=A*exp(-((z(k)-B)/C)**2)
          end do
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
       end do

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 


    case(igcss_warm7)
       !==============================================
       ! GCSS microphysics intercomparison Warm Rain 7
       !==============================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1:6)=(/ 1000., 200., 400., 170., 630.,&
            & 800. /)
       if (all(wctrl==0.))wctrl(1)=1.
       if (all(tctrl==0.))tctrl(1:2)=(/ 8*3600., 600./)
       if (all(pctrl_v==0.))pctrl_v(1)=.05 ! proxy for surface rain rate
       if (ipctrl==0)ipctrl=2
       maxZ=zctrl(1)
       maxT=tctrl(1)
       n_times=int(maxT/dt)

       Dw=zctrl(2)
       z0w=zctrl(3)
       Dq=zctrl(4)
       z0q=zctrl(5)
       zinv=zctrl(6)

       call set_DYCOMS_profile(maxZ)

       call allocate_forcing(nz,nx,n_times)

       do itime=1,n_times
          t=itime*dt
          time_in(itime)=t
          do j = 1, nx
             where(z<zinv)
                w_t(:,j,itime)=wctrl(1)*(z/(zinv))*(1.-exp(-((z-zinv)/Dw)**2))*sin(pi*t/tctrl(2))
             elsewhere
                w_t(:,j,itime)=0.
             end where
             where(z< z0q-1.25*Dq)
                qforce_in(:,j,itime)=cos(1.25*pi/2)
             elsewhere(z<z0q+Dq)
                qforce_in(:,j,itime)=cos(pi/2.*((z-z0q)/(Dq))) !&
             elsewhere
                qforce_in(:,j,itime)=0
             end where
             qforce_in(:,j,itime)=qforce_in(:,j,itime)*(pctrl_v(1)/3600.0) &
                  /sum(qforce_in(:,j,itime)*rho(:)*dz_half(:)) 

          end do
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
       end do

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 
    case(igcss_deep1)
       !================================================
       ! GCSS microphysics intercomparison Deep 1
       !================================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1:5)=(/ 15000., 5000., 5000., 7000., 7000. /)
       if (all(wctrl==0.))wctrl(1)=10.
       if (all(tctrl==0.))tctrl(1)=12*3600.
       if (all(pctrl_v==0.))pctrl_v(1)=5. ! proxy for surface rain rate
       if (ipctrl==0)ipctrl=7
       maxZ=zctrl(1)
       maxT=tctrl(1)
       n_times=int(maxT/dt)
       nt_in=max(100,n_times)
       dt_for=maxT/nt_in

       Dw=zctrl(2)
       z0w=zctrl(3)
       Dq=zctrl(4)
       z0q=zctrl(5)

       call set_standard_profile(ipctrl, maxZ)

       call allocate_forcing(nz,nx,nt_in)

       do itime=1,nt_in
          t=itime*dt_for
          time_in(itime)=t
          do j = 1,nx
             where(abs(z-z0w) < Dw)
                w_t(:,j,itime)=wctrl(1)*cos(pi/2.*((z(:)-z0w)/(Dw)))**4
             elsewhere
                w_t(:,j,itime)=0.
             end where
             where(abs(z-z0q) < Dq)
                qforce_in(:,j,itime)=cos(pi/2.*((z(:)-z0q)/(Dq)))**2
             elsewhere
                qforce_in(:,j,itime)=0
             end where
             qforce_in(:,j,itime)=qforce_in(:,j,itime)*(pctrl_v(1)/3600.0) &
                  /sum(qforce_in(:,j,itime)*rho(:)*dz_half(:))
          end do
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
       end do  
 
       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 
    case(igcss_deep2)
       !================================================
       ! GCSS microphysics intercomparison Deep 2
       !================================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1:3)=(/ 16000., 10000., .6 /)
       if (all(wctrl==0.))wctrl(1:2)=(/ 16., 6. /)
       if (all(tctrl==0.))tctrl(1:3)=(/ 7200., 1200., .2 /)
       if (ipctrl==0)ipctrl=7
       maxZ=zctrl(1)
       maxT=tctrl(1)
       n_times=int(maxT/dt)
       
       w0=wctrl(1)
       w1=wctrl(2)
       A=zctrl(3)
       B=tctrl(3)
       tup=tctrl(2)
       zup=zctrl(2)

       sigma=sqrt(log(w0/w1))

       call set_standard_profile(ipctrl, maxZ)

       call allocate_forcing(nz,nx,n_times)

       do itime=1,n_times
          t=itime*dt
          time_in(itime)=t
          do k=1,nz
             eta=z(k)/zup
             s=sqrt((t/tup+eta-1.)**2/(2.*A**2)+(eta-t/tup)**2/(2.*B**2))
             w_t(k,:,itime)=w0*exp(-(s/sigma)**2)
          end do
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
       end do 
 
       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 
    case(igcss_mixed1)
       !================================================
       ! GCSS microphysics intercomparison Mixed Phase 1
       !================================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1:2)=(/ 1000., 400. /)
       if (all(wctrl==0.))wctrl(1)=.3
       if (all(tctrl==0.))tctrl(1:2)=(/6*3600., 600./)
       if (ipctrl==0)ipctrl=5

       maxZ=zctrl(1)
       maxT=tctrl(1)
       n_times=int(maxT/dt)

       call set_standard_profile(ipctrl, maxZ)

       call allocate_forcing(nz,nx,n_times)

       do itime=1,n_times
          t=itime*dt
          time_in(itime)=t
          do j = 1,nx
             w_t(:,j,itime)=wctrl(1)*sin(pi*t/tctrl(2))*z(:)/zctrl(2)
          end do
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
       end do  

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 


    case(igcss_mixed3)
       !================================================
       ! GCSS microphysics intercomparison Mixed Phase 3
       !================================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1:6)=(/ 1000., 200., 400., 100., 350., 450. /)
       if (all(wctrl==0.))wctrl(1)=.3
       if (all(tctrl==0.))tctrl(1:2)=(/ 6*3600., 600. /)
       if (all(pctrl_v==0.))pctrl_v(1)=.05 ! proxy for surface rain rate
       if (ipctrl==0)ipctrl=5
       maxZ=zctrl(1)
       maxT=tctrl(1)
       n_times=int(maxT/dt)
       nt_in=max(100,n_times)
       dt_for=maxT/nt_in

       Dw=zctrl(2)
       z0w=zctrl(3)
       Dq=zctrl(4)
       z0q=zctrl(5)
       zinv=zctrl(6)

       call set_standard_profile(ipctrl, maxZ)

       call allocate_forcing(nz,nx,n_times)

       do itime=1,n_times
          t=itime*dt
          time_in(itime)=t
          do j = 1,nx
             where(z<zinv)
                w_t(:,j,itime)=wctrl(1)*(z/(zinv)) &
                     *(1.-exp(-((z-zinv)/Dw)**2))*sin(pi*t/tctrl(2))
             elsewhere
                w_t(:,j,itime)=0.
             end where
             where(z< z0q-1.25*Dq)
                qforce_in(:,j,itime)=cos(1.25*pi/2)
             elsewhere(z<z0q+Dq)
                qforce_in(:,j,itime)=cos(pi/2.*((z-z0q)/(Dq))) 
             elsewhere
                qforce_in(:,j,itime)=0
             end where
             qforce_in(:,j,itime)=qforce_in(:,j,itime)*(pctrl_v(1)/3600.0) &
                  /sum(qforce_in(:,j,itime)*rho(:)*dz_half(:)) 
          end do
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
    
       end do

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 
    case(iicel)
       !================================================
       ! ICE-L like simulation
       !================================================
       maxZ=5000.
       maxT=3600.
       n_times=int(maxT/dt)
       p_surf=50000.
       z_surf=(840.-0.01*p_surf)/0.06
       do k=1,nz
          z(k)=maxZ*k/float(nz)-25.
       end do
       do k=1,nz
          qv(k,:)=max(1e-8_wp,0.0052-5.9e-7*(z(k)+z_surf))
          theta(k,:)=(273.15+38.2-0.0091*(z(k)+z_surf))* &
               (1000./(840.-0.06*(z(k)+z_surf)))**r_on_cp
       end do
       
       call z2exner

       k=nz
       qv(k,:)=0.02*exp(-z(k)/1000.)

       call allocate_forcing(nz,nx,n_times)

       do itime=1,n_times
          t=itime*dt
          time_in(itime)=t
          do j=1,nx
             w_t(:,j,itime)=3.2*sin(2.*pi*t/800.)*exp(-t/800.)*rho(int(nz/2))/rho(:)
          end do
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
       end do 

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 

    case(ice1)
       !================================================
       ! Cirrus like simulation: colder and higher...
       !================================================
       ipctrl=9
       maxZ=12000.
       maxT=7200.
       n_times=int(maxT/dt)

       call set_standard_profile(ipctrl, maxZ)

       call allocate_forcing(nz,nx,n_times)

       do itime=1,n_times
          t=itime*dt
          time_in(itime)=t
          do j=1,nx
             w_t(:,j,itime)=2.*sin(2.*pi*t/1600.)*exp(-t/1600.)*rho(int(nz/2))/rho(:)
          end do
          do j = 0, nx+1 
            call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
          enddo
       end do 

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 


    case default
      print*, ''
      print*, '================================='
      print*, 'Warning: 1D case chosen'
      print*, 'Test case value not recognized - ',icase
      print*, 'Was a 2D case intended?'
      print*, 'Exitting...'
      print*, '================================='
      print*, ''
      stop
         
    end select

  end subroutine standard_cases
  
  subroutine z2p
    ! Calculate pressure and density from z and theta

    do j=1,nx
       call smooth1D(theta(:,j), theta_ref(:))
    end do
    
    rho(1)=p_surf/(Ru*theta_ref(1))
    exner(1,:)=(p_surf/p0)**r_on_cp-g*z(1)/(cp*theta_ref(1))
    
    do k=1,nz-1
       exner(k+1,:)=exner(k,:)- &
            g*(z(k+1)-z(k))/(cp*.5*(theta_ref(k)+theta_ref(k+1)))
       rho(k+1)=p0*exner(k+1,1)**(1./r_on_cp-1)/(Ru*theta_ref(k+1))
    end do
    
  end subroutine z2p

  subroutine z2exner
    ! Calculate anelasic exner profile from z grid
    ! and also density on full and half levels
    ! and generate a reference theta profile
    
    real(wp) :: dexner(nz)    ! delta exner across grid levels
    real(wp) :: dexner_column ! delta exner across column

    ! Save smoothed reference theta
    do j = 1, nx
       call smooth1D(theta(:,j), theta_ref(:))
    end do

    dexner(1)=g*z(1)/(cp*theta_ref(1))
    do k=2,nz
       ! calculate exner difference
       dexner(k)=g*(z(k)-z(k-1))/(cp*.5*(theta_ref(k)+theta_ref(k-1)))
    end do

    dexner_column=sum(dexner)
    
    exner(nz,:)=(p_surf/p0)**r_on_cp-dexner_column
    do k=nz-1,1,-1
       exner(k,:)=exner(k+1,:)+dexner(k)
    end do
    
    do j = 1,nx
       rho(:)=p0*exner(:,j)**(1./r_on_cp-1.)/(Ru*theta_ref(:))
    end do
    
  end subroutine z2exner

  subroutine set_standard_profile(ip, maxZ)

    integer, intent(in) :: ip ! case switch
    real(wp), intent(in) :: maxZ

    select case(ip)
    case (1)
       call set_RICO_profile(maxZ)
    case (2)
       call set_DYCOMS_profile(maxZ)
    case (3)
       call set_RICO_AS_profile(maxZ)
    case (4)
       call set_AS_updraught_profile(maxZ)       
    case (5)
       call set_SHEBA_profile(maxZ)
    case (6)
       call set_MPACE_profile(maxZ)
    case (7)
       call set_deep_profile(maxZ)
!++ag
    case (8)
       call set_low_updraught_profile(maxZ)   
    case (9)
       call set_cirrus_profile(maxZ)    
!--ag
    case default
       call set_RICO_profile(maxZ)
    end select

  end subroutine set_standard_profile

  subroutine set_RICO_profile(maxZ)
    !
    ! Set up the profile according to the RICO GCSS intercomparison 
    ! 
    real(wp), intent(in) :: maxZ

    allocate(pHeight(3))
    allocate(pTheta(3))
    allocate(pqv(3))
    allocate(theta_1d(nz))
    allocate(qv_1d(nz))
    
    pheight=(/ 0., 740., 3260./)
    ptheta=(/ 297.9, 297.9, 312.66 /)
    pqv=(/ .015, .0138, .0024 /)
    do k=1,nz
       z(k)=maxZ*k/float(nz)
    end do

    call interpolate(pHeight, ptheta, z, theta_1d)
    call interpolate(pHeight, pqv, z, qv_1d)   
    
    do j = 0,nx+1
       theta(:,j) = theta_1d(:)
       qv(:,j) = qv_1d(:)
    enddo

    where (qv<1e-6) qv=1e-6

    p_surf=p0
    call z2exner
    
    call calc_derived_fields
    
    deallocate(pqv)
    deallocate(pTheta)
    deallocate(pHeight)
    deallocate(theta_1d)
    deallocate(qv_1d)

  end subroutine set_RICO_profile

  subroutine set_DYCOMS_profile(maxZ)
    !
    ! Set up the profile according to the RICO GCSS intercomparison 
    ! 
    real(wp), intent(in) :: maxZ
    real(wp) :: zinv=800.  ! inversion height
    
    do k=1,nz
       z(k)=maxZ*k/float(nz)
    end do

    do k=1,nz
       if (z(k)<=zinv)then
          theta(k,:)=288.3
          qv(k,:)=.00945
       else
          theta(k,:)=295+(z(k)-zinv)**(1./3.)
          qv(k,:)=0.005 - 0.003*(1.-exp( (zinv-z(k)) /500.))
       end if
    end do

    p_surf=101780.
    call z2exner
    
    call calc_derived_fields
    

  end subroutine set_DYCOMS_profile

  subroutine set_RICO_AS_profile(maxZ)
    !
    ! Set up the profile according to the RICO 
    ! Abel & Shipway runs
    ! 
    real(wp), intent(in) :: maxZ

    allocate(pHeight(5))
    allocate(pTheta(5))
    allocate(pqv(5))
    allocate(qv_1d(nz))
    allocate(theta_1d(nz))
    pheight=(/ 0., 1000., 3500., 4625., 10000./)
    ptheta=(/ 298., 298., 310., 322., 337. /)
    pqv=(/ .014, .014, .0006, 0.0003, 0.00005 /)
    do k=1,nz
       z(k)=maxZ*k/float(nz)
    end do

    call interpolate(pHeight, ptheta, z, theta_1d)
    call interpolate(pHeight, pqv, z, qv_1d)

    do j=1,nx
       qv(:,j) = qv_1d(:)
       theta(:,j)=theta_1d(:)
    end do
    
    p_surf=p0
    call z2exner
    
    call calc_derived_fields
    
    deallocate(pqv)
    deallocate(pTheta)
    deallocate(pHeight)
    deallocate(qv_1d)
    deallocate(theta_1d)

  end subroutine set_RICO_AS_profile

  subroutine set_AS_updraught_profile(maxZ)
    !
    ! Set up the profile according to
    ! a saturated updraught from the RICO 
    ! Abel & Shipway runs
    ! 
    real(wp), intent(in) :: maxZ

    allocate(pHeight(5))
    allocate(pTheta(5))
    allocate(pRH(5))
    allocate(RH(nz))
    allocate(theta_1d(nz))

    pheight=(/ 0., 500.,  3500., 4000., 10000. /)
    ptheta=(/ 297.2, 297.2, 311.3, 315., 340. /)
    pRH=(/ .86, 1., 1., .5, .25 /)
    do k=1,nz
       z(k)=maxZ*k/float(nz)
    end do

    call interpolate(pHeight, ptheta, z, theta_1d)
    call interpolate(pHeight, pRH, z, RH)
    
    do j = 1,nx
       theta(:,j) = theta_1d(:)
    end do
    
    p_surf=p0
    call z2exner

    call calc_derived_fields

    do j=1,nx
       do k=1,nz
          qv(k,j)=RH(k)*qsaturation(TdegK(k,j),pmb(k,j))
       end do
    end do

    deallocate(RH)
    deallocate(theta_1d)
    deallocate(pRH)
    deallocate(pTheta)
    deallocate(pHeight)

  end subroutine set_AS_updraught_profile

!++ag
  subroutine set_low_updraught_profile(maxZ)
    !
    ! Set up the profile according to
    ! a saturated updraught...
    ! Andrew Gettelman
    ! 
    real(wp), intent(in) :: maxZ

    allocate(pHeight(3))
    allocate(pTheta(3))
    allocate(pRH(3))
    allocate(RH(nz))
    allocate(theta_1d(nz))

! 
!    pheight=(/ 0., 800.,  1500.,3000./)
!    ptheta=(/ 297.2, 297.2, 311.3, 315./)
!USE RICO values, but control RH

    pheight=(/ 0., 740., 3260./)
    ptheta=(/ 297.9, 297.9, 312.66 /)
    pRH=(/ .90, 1., .5 /)
    do k=1,nz
       z(k)=maxZ*k/float(nz)
    end do

    call interpolate(pHeight, ptheta, z, theta_1d)
    call interpolate(pHeight, pRH, z, RH)
    
    do j = 1,nx
       theta(:,j) = theta_1d(:)
    end do
    
    p_surf=p0
    call z2exner

    call calc_derived_fields

    do j=1,nx
       do k=1,nz
          qv(k,j)=RH(k)*qsaturation(TdegK(k,j),pmb(k,j))
       end do
    end do

    deallocate(RH)
    deallocate(theta_1d)
    deallocate(pRH)
    deallocate(pTheta)
    deallocate(pHeight)

  end subroutine set_low_updraught_profile

  subroutine set_cirrus_profile(maxZ)
    !
    ! Set up the profile according to
    ! a saturated updraught from the RICO 
    ! Abel & Shipway runs
    ! 
    real(wp), intent(in) :: maxZ

    allocate(pHeight(5))
    allocate(pTheta(5))
    allocate(pRH(5))
    allocate(RH(nz))
    allocate(theta_1d(nz))

    pheight=(/ 0., 5000.,  6000., 8000., 12000. /)
    ptheta=(/ 294., 280., 270., 270., 325. /)  !280, 300, 325.
    pRH=(/ .2, .2, 1., 1., .25 /)
    do k=1,nz
       z(k)=maxZ*k/float(nz)
    end do

    call interpolate(pHeight, ptheta, z, theta_1d)
    call interpolate(pHeight, pRH, z, RH)
    
    do j = 1,nx
       theta(:,j) = theta_1d(:)
    end do
    
    p_surf=p0
    call z2exner

    call calc_derived_fields

    do j=1,nx
       do k=1,nz
          qv(k,j)=RH(k)*qsaturation(TdegK(k,j),pmb(k,j))
       end do
    end do

    deallocate(RH)
    deallocate(theta_1d)
    deallocate(pRH)
    deallocate(pTheta)
    deallocate(pHeight)

  end subroutine set_cirrus_profile

!--ag

  subroutine set_SHEBA_profile(maxZ)
    !
    ! Set up the profile according to the RICO GCSS intercomparison 
    ! 
    real(wp), intent(in) :: maxZ

    allocate(pHeight(4))
    allocate(pTheta(4))
    allocate(pqv(4))
    allocate(theta_1d(nz))
    allocate(qv_1d(nz))
    
    pheight=(/ 0., 450., 480., 2000./)
    ptheta=(/ 257., 257.,  262.5, 272. /)
    pqv=(/ .000915, .000915, .0008, .00055 /)
    do k=1,nz
       z(k)=maxZ*k/float(nz)
    end do

    call interpolate(pHeight, ptheta, z, theta_1d)
    call interpolate(pHeight, pqv, z, qv_1d)

    do j=1,nx
       theta(:,j) = theta_1d(:)
       qv(:,j) = qv_1d(:)
    end do
    
    p_surf=p0
    call z2exner
    
    call calc_derived_fields
    
    deallocate(pqv)
    deallocate(pTheta)
    deallocate(pHeight)
    deallocate(theta_1d)
    deallocate(qv_1d)

  end subroutine set_SHEBA_profile

  subroutine set_MPACE_profile(maxZ)
    !
    ! Set up the profile according to the RICO GCSS intercomparison 
    ! 
    real(wp), intent(in) :: maxZ

    allocate(pHeight(4))
    allocate(pTheta(4))
    allocate(pqv(4))
    allocate(theta_1d(nz))
    allocate(qv_1d(nz))
    
    pheight=(/ 0., 1350., 1380., 4000./)
    ptheta=(/ 269., 269.,  273., 293. /)
    pqv=(/ .00192, .00192, .0008, .0003 /)
    do k=1,nz
       z(k)=maxZ*k/float(nz)
    end do

    call interpolate(pHeight, ptheta, z, theta_1d)
    call interpolate(pHeight, pqv, z, qv_1d)
    
    do j=1,nx
       theta(:,j) = theta_1d(:)
       qv(:,j) = qv_1d(:)
    end do

    p_surf=p0
    call z2exner
    
    call calc_derived_fields
    
    deallocate(pqv)
    deallocate(pTheta)
    deallocate(pHeight)
    deallocate(theta_1d)
    deallocate(qv_1d)

  end subroutine set_MPACE_profile


  subroutine set_deep_profile(maxZ)
    !
    ! Set up the profile roughly in accordance with KWAJEX data
    ! 
    real(wp), intent(in) :: maxZ
    real(wp) :: t_surf, q_surf

    allocate(theta_1d(nz))

    t_surf=300.
    p_surf=p0    
    q_surf=0.018

    do k=1,nz
       z(k)=maxZ*k/float(nz)
    end do
    
    call salr(z,t_surf,p_surf,q_surf,theta_1d)

    do j=1,nx
       theta(:,j)=theta_1d(:)
    end do

    call z2exner
    
    call calc_derived_fields
    
    do j=1,nx
       do k=1,nz
          if (q_surf < qsaturation(TdegK(k,j),pmb(k,j)))then
             qv(k,j)=q_surf
          else
             if (TdegK(k,j) > tK0c)then
                qv(k,j)=qsaturation(TdegK(k,j),pmb(k,j))
             else
                qv(k,j)=qisaturation(TdegK(k,j),pmb(k,j))
             end if
          end if
       end do
    end do

  end subroutine set_deep_profile

  subroutine allocate_forcing(nz,nx,ntimes)
    integer, intent(in) :: nz, nx, ntimes

    allocate(w_t(nz,0:nx+1,ntimes))
    allocate(w_t_half(nz,0:nx+1,ntimes))  
    allocate(wth_surf_t(0:nx+1,ntimes))
    allocate(wqv_surf_t(0:nx+1,ntimes))
    allocate(Tforce_in(nz,0:nx+1,ntimes))
    allocate(qforce_in(nz,0:nx+1,ntimes))
    allocate(time_in(ntimes))
    allocate(v_t(nz,0:nx+1,ntimes))
    allocate(v_t_half(nz,0:nx+1,ntimes))  

    w_t=0.
    v_t=0.
    w_t_half=0.
    v_t_half=0.
    wth_surf_t=0.
    wqv_surf_t=0.
    Tforce_in=0.
    qforce_in=0.
    
  end subroutine allocate_forcing

end module test_cases
 


