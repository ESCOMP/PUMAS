! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Fields
! 
! AH 2010 - Modified from 1-D KiD fields to be 2-D 
!
  Module column_variables

    Use typeKind

    Use parameters, only : nz, nx, nspecies, naerosol &
         , unset_real
  
    Use class_species, only : species

    Implicit None

    ! fields carried on each grid level
    real(wp) ::        &
         theta_ref(nz)    ! Smoothed initial theta profile used as 
                          ! a reference profile
    
    real(kind=wp) ::   &
          z(nz)            &   ! height (m)
         ,x(0:nx+1)            &   ! horizontal dimension (m)
         ,exner(nz, 0:nx+1)    &   ! exner pressure
         ,w(nz, 0:nx+1)    =0  &   ! vertical velocity (m/s)
         ,v(nz, 0:nx+1)    =0  &   ! horizontal velocity (m/s)
         ,rho(nz)          &   ! density (kg/m3)
         ,theta(nz, 0:nx+1)    &   ! potential temperature(K)
         ,qv(nz, 0:nx+1)       &   ! water vapour mixing ratio (kg/kg)
         ,ss(nz, 0:nx+1)   =0      ! supersaturation (transported as a 
                           ! passive scalar)

    real(kind=wp) ::     &
          z_half(nz)     &   ! height on half levels (m)
         ,w_half(nz,0:nx+1)  =0 &   ! vertical velocity on half levels(m/)
         ,rho_half(nz)=0  &   ! density on half levels(kg/m3)
         ,x_half(0:nx+1) & ! horizontal half levels (m)
         ,v_half(nz,0:nx+1)   ! horizontal winds on half horiz. grid (m/s)
 
!++ag
    ! extra interface level variables    
     real(kind=wp) ::                          &
       exner_half(0:nz, 0:nx+1)               &   ! exner pressure
       ,pmb_half(0:nz, 0:nx+1)  = unset_real    ! Pressure in mb on interface
!--ag
         
    !grid spacing
    real(kind=wp) ::      &
          dz(nz)          &   ! dz (m)
         ,dz_half(nz)     &   ! dz on half levels (m)
         ,dx(0:nx+1)          &   ! dx (m)
         ,dx_half(0:nx+1)


    real(wp) ::                  &
          pmb(nz, 0:nx+1)  =unset_real   & ! Pressure in mb
         ,TdegK(nz, 0:nx+1)=unset_real     ! Temperature in K

    type(species) ::                &
          hydrometeors(nz, 0:nx+1, nspecies) & ! hydrometeors
         ,aerosol(nz, 0:nx+1, naerosol)        ! aerosols

    ! Initial profiles
    real(wp) :: thinit(nz, 0:nx+1)   &   ! initial potential temperature 
               ,qvinit(nz, 0:nx+1)       ! initial water vapour  

    ! Applied forcing tendencies (if used)
    real(wp) ::          &    
          Tforce(nz,0:nx+1)=0   &   ! (Horizontal) Advective forcing:
                             ! temperature(K)
         ,qforce(nz,0:nx+1)=0       ! (Horizontal) Advective forcing:
                             ! vapour(kg/kg)

    real(wp) :: Trelax(nz)=1. ! relaxation timescale (if using relaxation)

    ! Tendencies calculated by the model (advection/microphysics)
    real(wp) ::            &
          dtheta_adv(nz, 0:nx+1)   & ! Advective tendency: Theta (K/s)
         ,dqv_adv(nz, 0:nx+1)      & ! Advective tendency: qv (kg/kg/s)
         ,dss_adv(nz, 0:nx+1)        ! Advective super sat tendency for bin micro
        
    type(species) ::       &
          dhydrometeors_adv(nz, 0:nx+1, nspecies) & ! Advective tendency: 
                                                ! hydrometeors (?/s)
         ,daerosol_adv(nz, 0:nx+1, naerosol)        ! Advective tendency: 
                                                ! aerosol (?/s)

    real(wp) ::              &
          dTheta_mphys(nz, 0:nx+1)   &   ! Mphys tendency: Theta (K/s)
         ,dqv_mphys(nz, 0:nx+1)      &   ! Mphys tendency: qv (kg/kg/s)
         ,dss_mphys(nz, 0:nx+1)          ! Mphys tendency: ss
    
    type(species) ::       &
          dhydrometeors_mphys(nz, 0:nx+1, nspecies) & ! Mphys tendency: 
                                             ! hydrometeors (?/s)
         ,daerosol_mphys(nz, 0:nx+1, naerosol)        ! Mphys tendency: 
                                             ! aerosol (?/s)

    real(wp) ::            &
          dTheta_div(nz, 0:nx+1)   &   ! Divergence tendency: Theta (K/s)
         ,dqv_div(nz, 0:nx+1)      &   ! Divergence tendency: qv (kg/kg/s)
         ,dss_div(nz, 0:nx+1)          ! Divergence tendency: ss
    
    type(species) ::       &
          dhydrometeors_div(nz, 0:nx+1, nspecies) & ! Divergence tendency: 
                                           ! hydrometeors (?/s)
         ,daerosol_div(nz, 0:nx+1, naerosol)        ! Divergence tendency: 
                                           ! aerosol (?/s)

    ! Some test cases use initial values for hydrometeors as follows
    type(species) ::       &
          hydrometeors_init(nz,nspecies) & ! Initial values for
                                           ! hydrometeors 
         ,aerosol_init(nz,naerosol)        ! Initial values for
                                           ! aerosol (?/s)

! NB Forcing of variables not yet full functional
    type(species) ::       &
          dhydrometeors_force(nz,0:nx+1,nspecies) & ! Imposed forcing: 
                                           ! hydrometeors (?/s)
         ,daerosol_force(nz,0:nx+1,naerosol)        ! Imposed forcing: 
                                           ! aerosol (?/s)

    ! scalars 
    real(wp) ::      &
          p_surf     &   ! surface pressure (Pa)
         ,z_surf=0   &   ! surface height (m) 
         ,wth_surf(0:nx+1)=0 &   ! surface potential temperature flux (Km/s)
         ,wqv_surf(0:nx+1)=0     ! surface vapour flux (kg/kg m/s)

    ! Mask for applying increments. Increments are not applied where 
    ! field is zero (and are multiplied by weight otherwise).

    real(wp) :: field_mask(nz,0:nx+1)=1.  ! Mask for increments
    
  end Module column_variables
    
