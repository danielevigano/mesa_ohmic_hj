! ***********************************************************************
!
!   Copyright (C) 2011  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff
      
      implicit none
      
      include "test_suite_extras_def.inc"

      contains

      include "test_suite_extras.inc"
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         ! Point to the new heat source subroutine (also in this file).
         s% other_energy => ohmic_heating
!         s% other_energy_implicit => ohmic_heating
      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step

	   ! 1d Interpolation function
      real(dp) function interp(a,b,c,d,e)
         implicit none
         real(dp) :: slope,intercept
         real(dp) :: a,b,c,d,e

         slope= (d-b)/(c-a)
         intercept = b - slope*a
         interp = intercept + slope*e
      end function interp
  		
  	   ! 2d interpolation function
      real(dp) function interp2d(x1,x2,y1,y2,z1,z2,z3,z4,xx,yy)
         implicit none
         real(dp) :: x1,x2,y1,y2,z1,z2,z3,z4,xx,yy
         real(dp) :: z_int1, z_int2
		
         z_int1 = interp(x1,z1,x2,z2,xx)
         z_int2 = interp(x1,z3,x2,z4,xx)
         interp2d = interp(y1,z_int1,y2,z_int2,yy)
      end function interp2d

      subroutine get_j(ndim,n_pieces,J0cgs,sigma,katm,rr_dynamo,rad,lmin,lmax,sigma_fit,jind,beta,alpha)

         implicit none
         integer, intent(in) :: ndim, n_pieces, katm, lmin, lmax
         real(dp), intent(in) :: J0cgs, rr_dynamo
         real(dp), intent(in), dimension(ndim) :: sigma, rad
         real(dp), intent(inout), dimension(ndim) :: sigma_fit, jind

         integer :: count, ll, p, k, kmin, kmax
         integer, dimension(0:n_pieces) :: k_pieces
         real(dp), dimension(n_pieces), intent(out) :: beta
         real(dp), dimension(:,:), allocatable, intent(out) :: alpha
         real(dp) :: avx, avy, num, den, jnorm, sigmanorm, wj

         allocate(alpha(n_pieces,lmax-lmin+1))

         beta(:) = 0d0
         alpha(:,:) = 0d0
         sigma_fit(:) = 0d0
         jind(:) = 0d0
         
         jnorm = J0cgs

         do k = 1,katm
            jind(k) = jnorm
            sigma_fit(k) = sigma(k)
         enddo

         k_pieces(0) = katm
         k_pieces(n_pieces) = ndim
         do p = 1,n_pieces
            k_pieces(p) = katm + (p-1)*(ndim-katm)/(n_pieces-1)
         enddo

         do p = 1,n_pieces-1
            avx = 0d0
            avy = 0d0
            count = 0
            num = 0d0
            den = 0d0
            kmin = k_pieces(p)
            kmax = k_pieces(p+1)

            ! compute best-fit slope of the trend
            ! y = -beta*x + q
            ! log10(sigma) = -beta*log10(rad) + q
            do k = kmin, kmax
               count = count + 1
               avx = avx + dlog10(rad(k))
               avy = avy + dlog10(sigma(k))
            enddo
            avx = avx/count
            avy = avy/count

            ! Standard linear regression fit
            do k = kmin, kmax
               num = num + (dlog10(rad(k)) - avx)*(dlog10(sigma(k)) - avy)
               den = den + (dlog10(rad(k)) - avx)**2
            enddo
            beta(p) = - num/den
            sigmanorm = 10d0**(avy + beta(p)*avx)

            do k = kmin, kmax
               ! sigma_fit itself is not used, but we print it for output purposes
               ! sigma_fit = 10^q*rad^(-beta)
               sigma_fit(k) = sigmanorm * rad(k)**(-beta(p))
               jind(k) = 0d0

               do ll = lmin,lmax
                  ! alpha is the power law of the analytic solution to the Laplace equation for the current field amplitude
                  alpha(p,ll-lmin+1) = 0.5d0*(beta(p) - 1 + dsqrt((1d0-beta(p))**2 + 4d0*ll*(ll+1)))
   
                  ! We consider that at r_dynamo the magnetic spectrum is flat, so that in the atmosphere, where the normalization is taken,
                  ! the lowest multipoles l_b are the largest:  B \propto (r_dynamo/r_atm)^(l_b+2)
                  ! The v x B contribution to J follow the magnetic spectrum weight, shifted in l by the wind contribution (usually l_wind=1):
                  ! l_vxb \sim l_b + l_wind
                  ! lmin = min(l_b) + l_wind
                  ! Therefore, the relative contributions of each multipole, for a fixed profile goes with (rr_dynamo/rad(katm))^ll
                  wj = (rr_dynamo/rad(katm))**ll
                  jind(k) = jind(k) + ( wj * rad(k)**(- beta(p) - 1 + alpha(p,ll-lmin+1)) )

               enddo
   
            enddo
         
            ! Renormalize to guarantee continuity with the previous range
            jind(kmin:kmax) = jind(kmin:kmax) * jnorm/jind(kmin)
            jnorm = jind(kmax)
         enddo

      end subroutine get_j

!-----------------------------------------------------------------------
! Joule heating.
!
! Notes:
! - Follow the instructions provided in the README file located in the
! MESA folder /star/other/.
! - The module const_def is defined in /const/public/const_def.f90 and
! contains standard definitions of constants.
! - The derived type star_info (of the variable s) is defined in the
! file /star_data/public/star_data_def.f90 (as well as the subroutine
! star_ptr). Most declarations of the derived type are in the included
! file star_data.inc in the same folder.
! - Surface irradiation from the host star is controlled by the input
! variables irradiation_flux and column_depth_for_irradiation (defined
! in inlist_evolution). The input variables are processed by the file
! /star/private/ctrls_io.f90 where the "controls" namelist is defined.
!-----------------------------------------------------------------------
      subroutine ohmic_heating(id, ierr)
            use const_def, only: Lsun,Msun,Rsun,m_jupiter,r_jupiter, &
                              boltz_sigma,clight,me,dp,pi 

         ! Default definitions.
         implicit none
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ! Constants
         real(dp), parameter :: m_u = 1.660539d-24  ! Atomic mass unit (in g).
         real(dp), parameter :: e_e = 4.8032d-10 ! Elementary charge in cgs units (g^0.5 cm^1.5 s^-1).

         ! Conversion factors
         real(dp), parameter :: sigma_SI_to_cgs = clight**2d0 / 1d11 ! Conversion multiplicative factor for conductivity from cgs (s^-1) to S/m (SI).
         real(dp), parameter :: J_SI_to_cgs = clight/1e5
         real(dp), parameter :: p_cgs_to_bar = 1d-6 ! Pressure (1 barye = 1e-6 bar).

         ! Custom (local) definitions.
         integer :: k,ndim
         real(dp) :: age,density,pressure,temperature,entropy,gravity
         real(dp) :: radius,rr,mass,dmass,luminosity       ! Variables of the structure equations
         real(dp) :: radius_p,luminosity_p                 ! Planetary radius (cm) quantities (in cgs).
         real(dp), save :: X_frac, Y_frac, Z_frac, mass_p ! Composition and mass (fixed)

         ! Irradiation and Ohmic terms
         real(dp) :: eps_irrad ! Specific heat from irradiation and Joule effect (erg/g/s)
         real(dp) :: lum_irrad,lum_joule,lum_joule_above,lum_joule_below ! Integrated luminosity for irradiation, Joule, Joule above/below RCB (erg/s)
         real(dp), save :: irradiation,column_depth_for_irradiation,Teq
         real(dp) :: x_e,sv_e,f_time_ohm

         ! Parameters read from the input (see below for their description)
         real(dp), save :: v_avg, p_atm, p_dyn, pcmin, s0
         real(dp), save :: age_ohmic_on, age_ohmic_on_full
         real(dp), save :: fdip, fohm, f_avg, w_age_max
         integer, save :: lmin, lmax
         character(len=1), save :: scaling_b

         ! Joule heating quantities
         integer :: katm, kdyn                               ! Index of radius identifying the wind region and the dynamo surface
         real(dp), dimension(:), allocatable :: eps_joule, sigma_thermal, sigma_free_e, sigma_compl, sigma_fit, jind, scale_height, scale_height_t, q_c, dr
         integer, parameter :: nsig = 3000   ! Set it high enough (> ndim), it's used to save the previous Joule heating   
         real(dp), dimension(nsig), save :: joule_heating_prev
         real(dp), save :: mu_mean, age_prev, Bdips_prev
         real(dp) :: joule_heating, window, norm_J, Bdips
         real(dp) :: t_avg, w_age, sigma_atm
         real(dp) :: rho_mean, volume
         real(dp) :: r_dynamo, B_dynamo, B_flux, B_reiners, B_dip_surface       ! Inferred values of B dynamo and B_dip_surface
         real(dp) :: luminosity_int ! Magnetic field and internal luminosity just above the dynamo region (1e6 bar).
         real(dp) :: q_0, F_factor   ! Factors entering in the scaling laws

  		   ! For the low-density conductivity tables
         integer :: i,j,p1,p2,q1,q2
         integer, parameter :: ndens_ku = 6, ntemp_ku = 700
         real(dp), dimension(ndens_ku),save :: rho_ku = (/1e-7,1e-6,1e-5,1e-4,1e-3,1e-2/)
         real(dp), dimension(ntemp_ku),save :: T_ku
         real(dp), dimension(ndens_ku,ntemp_ku), save :: sigma_ku
         real(dp) :: z1,z2,z3,z4,sigma_ku_intp,sigma_bo_intp
         character (len=10):: rho_str(ndens_ku)=(/"1e-7.dat","1e-6.dat","1e-5.dat","1e-4.dat","1e-3.dat","1e-2.dat"/)
         real(dp), parameter :: a_K = 1d-7
         
		   ! For the high-density conductivity tables (Bonitz et al. 2024)
         integer, parameter :: ndens_bo = 10, ntemp_bo = 6
         real(dp), dimension(ndens_bo), save :: rho_bo
         real(dp), dimension(ntemp_bo) :: T_bo = (/5000.0,10000.0,15000.0,20000.0,30000.0,50000.0/)
         real(dp), dimension(ndens_bo,ntemp_bo), save :: sigma_bo
         character (len=10):: temp_str(ntemp_bo)=(/"5000.K","10000K","15000K","20000K","30000K","50000K"/)
		
         ! Pieces for the piecewise power-law reconstruction of the conductivity profile
         integer, parameter :: n_pieces = 50
         ! Magnetic multipoles considered in the induction solution
         real(dp), dimension(n_pieces) :: beta
         real(dp), dimension(:,:), allocatable :: alpha

         ! Pressures and radii of the main (1) and superficial (2, if any)
         ! radiative-convective boundaries (RCB).
         real(dp) :: r_top_rcb1,r_top_rcb2,r_bot_rcb2
         real(dp) :: p_top_rcb1,p_top_rcb2,p_bot_rcb2
         real(dp) :: r_top_rcb1_old,r_top_rcb2_old,r_bot_rcb2_old
         real(dp) :: p_top_rcb1_old,p_top_rcb2_old,p_bot_rcb2_old
         real(dp) :: p_surface
         ! index for the presence(1)/absence(0) of the convective region
         integer :: conv_sup

         ! Variables for the outputs
         logical, save :: first_call = .TRUE.
         logical, save :: write_profile = .FALSE.
         integer, save :: ia = 1
         integer :: iaa
         ! Maximum age of the output and number of outputs
         ! (set like this to have the same f3.1 format)
         integer, parameter :: nage = 30
         real(dp), parameter :: max_age_out = 9.9d9
         real(dp), dimension(nage), save :: age_output
         character(len = 7), parameter :: AGE_STR_FORMAT = "(f3.1)"
         character(len = 3) :: str
         character(len = 3) :: tmp
         integer,save :: iage
         
         !--------------------------------------------------------------
         ! Initialize the pointer s (of derived type star_info).
         ! Before the call to star_ptr the pointer s is uninitialized.
         ! Note: Assigning uninitialized pointers to variables can cause
         ! segmentation fault.
         !--------------------------------------------------------------
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! After the call to star_ptr the pointer s is initialized.
         ndim = s% nz ! Number of radial zones (variable).
         age = s% star_age ! Current age [yrs]

         if (first_call) then
            iage = 1
            do iaa = 1,nage
               age_output(iaa) = iaa * (max_age_out) / nage
            end do
            age_prev = 0d0
            Bdips_prev = 0d0
    
            mass_p = s% m(1)   ! Total mass of the planet (outermost value, constant in time)
            irradiation = s% irradiation_flux      ! Stellar irradiation flux [erg cm^-2 s^-1]
            column_depth_for_irradiation = s% column_depth_for_irradiation  ! Column depth for irradiation deposition [cm^2 g^-1]
            ! Composition (constant in time and space)
            X_frac = s% X(1)   ! Mass fraction of neutral Hydrogen
            Y_frac = s% Y(1)   ! Mass fraction of neutral Helium
            Z_frac = s% Z(1)   ! Mass fraction of neutral heavier elements
            ! Derived quantities
            ! Mean molecular weight, assuming for heavy elements Z/A=1/2 
            ! and negligible H-dissociation and negligible ionization
            mu_mean = 1d0/(0.5d0*X_frac + 0.25d0*Y_frac + Z_frac*0.5d0)
            ! Equilibrium temperature (K) assuming perfect redistribution
            Teq = (irradiation/4d0/boltz_sigma)**0.25d0
            ! Read Ohmic input parameters of the models.
            open(1,file = "inlist_ohmic")
            read(1,*) v_avg     ! average velocity in the wind layer [cm s^-1]
            read(1,*) lmin      ! minimum l of the (vxB) term decomposition
            read(1,*) lmax      ! maximum l of the (vxB) term decomposition
            read(1,*) p_atm     ! Pressure considered for the wind layers [bar]
            read(1,*) p_dyn     ! Pressure considered to evaluate the dynamo quantities [bar]
            read(1,*) scaling_b ! Flag for scaling law: Reiners et al. 2009 'R', Christensen et al. 2009 'C'
            read(1,*) fohm      ! Fraction of Ohmically dissipated power over the total, for magnetic field scaling laws
            read(1,*) fdip      ! Dipolarity at the dynamo surface, used to infer the surface field.
            read(1,*) age_ohmic_on       ! Age at which we gradually start switching on the Ohmic term - std: 1d5
            read(1,*) age_ohmic_on_full  ! Age from which the Ohmic term is fully on - std: 1e6
            read(1,*) f_avg     ! Fraction of age over which taking the time-average of the Ohmic terms, for stability purposes (0 for no averages in time)
            read(1,*) w_age_max ! Maximum weight of the current value of Ohmic heating to give in the averages (1 for no averages in time)
            read(1,*) pcmin     ! Mpressure at which the Ohmic heating is manually damped outward [bar]. Values <= 0 for not having it.
            read(1,*) s0        ! Width of the tanh of the deposition region, i.e., how smooth is the beginning and end of the region (log(p)) - std: 0.25
            close(1)

            ! Sanity checks for the Ohmic model input
            if (v_avg .lt. 0) then
               write(6,*) "Input ERROR: need v_avg >=0. vavg = ",v_avg
               stop
            endif
            if ((lmin .lt. 2) .or. (lmax .lt. lmin)) then
               write(6,*) "Input ERROR: need lmin >=2 and lmax>lmin. lmin, lmax = ",lmin,lmax
               stop
            endif
            if ((p_atm .lt. 1) .or. (p_atm .gt. 1e3) .or. (p_dyn .lt. p_atm)) then
               write(6,*) "Input ERROR: realistic values of p_atm \in [1,1000], and p_dyn > p_atm. p_atm, p_dyn = ",p_atm, p_dyn
               stop
            endif
            if ((scaling_b .ne. 'R') .and. (scaling_b .ne. 'C')) then
               write(6,*) "Input ERROR: need scaling_b = 'R' (Reiners et al. 2009) or 'C' (Christensen 2009). scaling_b = ",scaling_b
               stop
            endif
            if (age_ohmic_on_full .lt. age_ohmic_on) then
               write(6,*) "Input ERROR: need age_ohmic_on_full >= age_ohmic_on: ",age_ohmic_on,age_ohmic_on_full
               stop
            endif
            if ((f_avg .lt. 0) .or. (f_avg .ge. 1) .or. (w_age_max .le. 0) .or. (w_age_max .gt. 1)) then
               write(6,*) "Input ERROR: need f_avg >=0 and < 1; w_age_max > 0 and <=1. f_avg,w_age_max = ",f_avg,w_age_max
               stop
            endif
            if ((pcmin .gt. 0) .and. (s0 .le. 0)) then
               write(6,*) "Input ERROR: s0 needs to be > 0 if pcmin > 0. pcmin, s0 = ",pcmin,s0
               stop
            endif

            ! Convert pressures from bar to cgs
            p_atm = p_atm/p_cgs_to_bar
            p_dyn = p_dyn/p_cgs_to_bar
            pcmin = pcmin/p_cgs_to_bar

 	         ! Read the low-density conductivity data from the Kumar table (Kumar et al. 2021)
            do j=1,ndens_ku
               open(21,file="sigma_low_density/T_EC_"//rho_str(j))
               do i=1,ntemp_ku
                  read(21,*)T_ku(i),sigma_ku(j,i)
               enddo
               close(21)
            enddo
            ! Read the high-density conductivity data from Bonitz et al. 2024
            do j=1,ntemp_bo
               open(22,file="sigma_high_density/B_data_cond_"//temp_str(j))
               do i=1,ndens_bo
                  read(22,*)rho_bo(i),sigma_bo(i,j)
               enddo
               close(22)
            enddo

            joule_heating_prev(:) = -1d0

         endif ! Close the first-call loop
         
         ! Planetary luminosity, mass and radius (in cgs).
         ! The surface values are at index 1.
         luminosity_p = s% L(1)
         radius_p = s% r(1)

         ! Initialize variables
         joule_heating = 0d0
         allocate(eps_joule(ndim))
         allocate(dr(ndim))
         allocate(sigma_thermal(ndim))
         allocate(sigma_free_e(ndim))
         allocate(sigma_compl(ndim))
         allocate(sigma_fit(ndim))
         allocate(jind(ndim))
         allocate(q_c(ndim))
         allocate(scale_height(ndim))
         allocate(scale_height_t(ndim))
         dr(:) = 0d0
         sigma_thermal(:) = 0d0
         sigma_free_e(:) = 0d0
         sigma_compl(:) = 0d0
         sigma_fit(:) = 0d0
         jind(:) = 0d0
         q_c(:) = 0d0
         scale_height(:) = 0d0
         scale_height_t(:) = 0d0
         lum_irrad = 0d0
         lum_joule = 0d0
         lum_joule_below = 0d0
         lum_joule_above = 0d0
         luminosity_int = 0d0
         sigma_atm = 0d0

         !--------------------------------------------------------------
         ! Determine the pressure at the RCB (linear interpolation).
         ! Here we track the radius and pressure delimiting the main 
         ! and superficial (if any) convective regions
         ! The main region bottom is the deepest layer by definition
         !--------------------------------------------------------------
         ! Define the surface pressure (at outermost cell).
         p_surface = s% peos(1)
         ! Initial values.
         r_top_rcb1_old = radius_p
         r_top_rcb2_old = radius_p
         r_bot_rcb2_old = radius_p
         p_top_rcb1_old = p_surface
         p_top_rcb2_old = p_surface
         p_bot_rcb2_old = p_surface
         r_top_rcb1 = (s% conv_mx1_top_r) * Rsun   ! top of the main convective region, converted to cm
         r_top_rcb2 = (s% conv_mx2_top_r) * Rsun   ! top of the possible second convective region (cm)
         r_bot_rcb2 = (s% conv_mx2_bot_r) * Rsun   ! bottom of the possible second convective region (cm)
         p_top_rcb1 = 0d0   ! corresponding pressures
         p_top_rcb2 = 0d0
         p_bot_rcb2 = 0d0
         conv_sup = 0     ! 1 if there is a secondary convective region


         do k = 1,ndim
            ! Determine the luminosity and radius of the top of the dynamo region
            if (s% peos(k) .le. p_dyn) kdyn = k

            ! Identify the index of the pressure at which we consider no more relevant wind
            if (s% peos(k) .le. p_atm) katm = k
         enddo

         r_dynamo = s% r(kdyn)
         luminosity_int = s% L(kdyn)

         ! Loop to calculate the pressures of the convective regions
         ! Loop over radial index. (The radial dimension ndim is variable.)
         ! From outside inwards (k=1 corresponds to the surface).
         do k = 1,ndim
            radius = s% r(k)
            pressure = s% peos(k) ! Total pressure (P = P_rad + P_gas).
            luminosity = s% L(k)
            if (k .gt. 1 .and. k .lt. ndim) dr(k) = 0.5*(s%r(k-1) - s%r(k+1))

            scale_height(k) = min(s%scale_height(k), s%r(kdyn) - s%r(ndim-1))

            ! main convection region (extending from r_top_rcb1 down to the center).
            if(r_top_rcb1_old.gt.r_top_rcb1.and.radius.le.r_top_rcb1)then
               p_top_rcb1 = p_top_rcb1_old + (pressure - p_top_rcb1_old) * (r_top_rcb1 - r_top_rcb1_old) / (radius - r_top_rcb1_old)
            endif
            r_top_rcb1_old = radius
            p_top_rcb1_old = pressure

            ! Secondary convection region (from r_top_rcb2 down to r_bot_rcb2).
            if(r_top_rcb2.ne.0d0.and.r_top_rcb2_old.gt.r_top_rcb2.and.radius.le.r_top_rcb2)then
               p_top_rcb2 = p_top_rcb2_old + (pressure - p_top_rcb2_old) * (r_top_rcb2 - r_top_rcb2_old) / (radius - r_top_rcb2_old)
               conv_sup = 1 ! It's one if there is a secondary convective region
            endif
            r_top_rcb2_old = radius
            p_top_rcb2_old = pressure

            if(r_bot_rcb2.ne.0d0.and.r_bot_rcb2_old.gt.r_bot_rcb2.and.radius.le.r_bot_rcb2)then
               p_bot_rcb2 = p_bot_rcb2_old + (pressure - p_bot_rcb2_old) * (r_bot_rcb2 - r_bot_rcb2_old) / (radius - r_bot_rcb2_old)
            endif
            r_bot_rcb2_old = radius
            p_bot_rcb2_old = pressure

         enddo

         !--------------------------------------------------------------
         ! Magnetic field strength from scaling laws.
         !--------------------------------------------------------------


         q_c = 2 * s%cp * s%t * s%rho**2 * s%conv_vel**3 / s%peos / s%chiT * s%chiRho ! convective flux vector
         q_0 = q_c(kdyn) 
         volume = 4 * pi  * (s%r(kdyn)**3 - s%r(ndim-1)**3) / 3 ! dynamo shell volume
         rho_mean = (s%m(kdyn) -  s%m(ndim-1)) / volume
         scale_height_t = s%peos / (s%rho * s%grav * s%grada)

         ! F factor appearing in Christensen 2009
         F_factor = sum((scale_height(kdyn:ndim-1)/scale_height_t(kdyn:ndim-1)*q_c(kdyn:ndim-1)/q_0)**(2d0/3d0) * &
                         (s%rho(kdyn:ndim-1)/rho_mean)**(1d0/3d0) * s%r(kdyn:ndim-1)**2 * dr(kdyn:ndim-1))
         F_factor = sum((scale_height(kdyn:ndim-1)/scale_height_t(kdyn:ndim-1)*q_c(kdyn:ndim-1))**(2d0/3d0) * &
                         (s%rho(kdyn:ndim-1)/rho_mean)**(1d0/3d0) * s%r(kdyn:ndim-1)**2 * dr(kdyn:ndim-1))
         F_factor = (4 * pi * F_factor / volume)**(3d0/2d0)

         ! Christensen 2009 scaling law:
         ! <B>**2/(2 mu0) = c fohm <rho>**(1/3) (F q0)**(2/3)
         ! In cgs: <B>**2/(8 pi) = c fohm <rho>**(1/3) (F q0)**(2/3)
         ! where we set c = 0.63 
         B_flux = dsqrt(8d0 * pi * fohm * 0.63 * (rho_mean)**(1.0/3.0) * (F_factor * q_0)**(2.0/3.0))  ! In cgs
         B_flux = dsqrt(8d0 * pi * fohm * 0.63 * (rho_mean)**(1.0/3.0) * (F_factor)**(2.0/3.0))  ! In cgs

         ! Dynamo magnetic field [in Gauss] from Reiners & Christensen 2010, eq. 1
         ! We use the internal lumonisity at p=p_dyn, even though the original formula (for stars) uses the surface luminosity, not applicable for irradiated HJs
         B_reiners = 4.8d3*((mass_p/Msun)*(luminosity_int/Lsun)**2d0/(radius_p/Rsun)**7d0)**(1d0/6d0)

         ! Extrapolation of the dipolar component to the surface of the planet,
         ! assuming a dynamo surface at p_dyn
         ! B_dip_surface = fdip B_dynamo (r_dynamo/radius_p)**3

         if (scaling_b .eq. 'C') then
            B_dynamo = B_flux
         elseif (scaling_b .eq. 'R') then
            B_dynamo = B_reiners
         endif
         B_dip_surface = fdip * B_dynamo * (r_dynamo/radius_p)**3d0

         if (B_flux .ne. B_flux) then
            print*,B_flux,F_factor,q_0,rho_mean,sum(scale_height(kdyn:ndim-1)/scale_height_t(kdyn:ndim-1)),sum(q_c(kdyn:ndim-1)/q_0)
            print*,kdyn,s% peos(kdyn)*p_cgs_to_bar
            print*,s% chiT
            print*,"ERROR IN COMPUTING B_FLUX"
            stop
         endif

         !--------------------------------------------------------------
         ! Radial profiles header
         !--------------------------------------------------------------
         ! Writing the header
         if(ia .le. nage .and. age >= age_output(ia) .and. age <= max_age_out)then
            write(tmp, AGE_STR_FORMAT) age/1d9
            str = trim(tmp)
            open(2,file = "LOGS/radial_profile_"//str//"Gyr.txt")
 
            ! Header with general information               
            write(2,'(a1,6a18)')"#","Age[yr]","Teq[K]","Mass[Mj]","Radius[Rj]","R_rcb[Rj]","P_rcb[bar]"
            write(2,'(a1,6es18.9)')"#",age,Teq,mass_p/m_jupiter, &
            radius_p/r_jupiter,r_top_rcb1/r_jupiter,p_top_rcb1*p_cgs_to_bar
            
            write(2,'(a1,a6,25a18)')"#","1.Index","2.Radius","3.Radius", &
            "4.Density","5.Pressure","6.Temperature","7.Entropy","8.Gravity", &
            "9.Mass","10.dm","11.Chi_rho","12.Chi_T","13.Gamma_1","14.Gamma_3", &
            "15.Lum_int","16.Spec_Irr_heat","17.Spec_Ohm_heat","18.Vol_Ohm_heat", &
            "19.Sigma K+","20.Sigma complete","21.Sigma fit","22.J_induced"
            
            write(2,'(a1,a6,25a18)')"#","","[Rj]","[Rp]", &
            "[g cm^-3]","[bar]","[K]","Entropy","[cm s^-2]", &
            "[Mj]","[Mj]","Chi_rho","Chi_T","Gamma_1","Gamma_3", &
            "[erg s^-1]","[erg s^-1 g^-1]","[erg s^-1 g^-1]",&
            "[erg s^-1 cm^-3]","[S/m]","[S/m]","[S/m]","[A/m^2]"

            write_profile = .TRUE.
         endif

         ! Defining F(t): gradual increase of the Ohmic term to avoid
         ! unrealistic higher values of J at early age (depending on the used normalization)
         f_time_ohm = 1d0
         if (age .lt. age_ohmic_on_full .and. age .gt. age_ohmic_on) then
            f_time_ohm = (age - age_ohmic_on)/(age_ohmic_on_full-age_ohmic_on)
         elseif (age .lt. age_ohmic_on_full) then
            f_time_ohm = 0d0
         endif

         !--------------------------------------------------------------
         ! main loop for Ohmic heating terms.
         !--------------------------------------------------------------
         ! Loop over radial index. (The radial dimension ndim is variable.)
         ! From outside inwards (k=1 corresponds to the surface).
         do k = 1,ndim
            ! All quantities are in cgs units.
            radius = s% r(k)
            rr = radius/radius_p ! Dimensionless radial coordinate.
            mass = s% m(k) ! Mass enclosed [in g].
            dmass = s% dm(k) ! Baryonic mass of cell k.
            density = s% rho(k)
            ! There is also s% entropy(k), which seems to be something different.
            entropy = exp(s% lnS(k)) ! Log of S (specific entropy).
            ! Gravitational acceleration: g = standard_cgrav*mass/radius**2.
            gravity = s% grav(k)
            pressure = s% peos(k) ! Total pressure (P = P_rad + P_gas).
            temperature = s% T(k)
            luminosity = s% L(k) ! For output purpose only.

            !-----------------------------------------------------------
            ! Electrical conductivity profile sigma (s^-1).
            ! Recipe used in Perna et al. 2010 and following works
            !-----------------------------------------------------------
            
            ! Fraction of electrons assuming only contributions from Potassium, Perna et al. 2010 Eq.1
            x_e = 6.47d-13 * (a_K * 1d7)**0.5d0 * (temperature / 1d3)**0.75d0 &
                  * (2.4d15 / (density / mu_mean / m_u))**0.5d0 &
                  * exp(-25188d0 / temperature) / 1.15d-11

            ! Eq.2 (See eq.14 of Draine et al. 1983 for the correct units!)
            sv_e = 1d-15 * (128d0 * boltzm * temperature / (9d0 * pi * me))**0.5d0
            ! Conductivity (eq.3) (s^-1)
            sigma_thermal(k) = x_e * e_e**2 / me / sv_e

            ! Contribution from the fraction of electrons per nucleon coming from pressure ionization
            sigma_free_e(k) = exp(s% lnfree_e(k)) * e_e**2 / me / sv_e

            if (age .gt. age_ohmic_on) then
            ! Interpolation is done for the age above which the Ohmic heating is activated 

               !-----------------------------------------------------------
               ! Electrical conductivity profile sigma (s^-1).
               ! Data from Kumar 
               ! Useful in the following range:
               ! density 1e-7 to 1e-2
               ! temperature 50 K - 35000 K
               !-----------------------------------------------------------
!               if (density .lt. rho_ku(1) .or. temperature .lt. T_ku(1)) then
!                  print*,"Density: ",density,"Temperature: ",temperature
!                  print*,"WARNING: Density and/or temperature is less than the available tables (Kumar et al. 2021)"
!               endif
            
               if (density .le. rho_ku(ndens_ku)) then
               ! Use the table of Kumar
             
                  ! Locate closest densities on the grid
                  i=1
                  do while (density .gt. rho_ku(i+1))
                     i=i+1
                  end do
                  p1=i
                  p2=i+1

                  ! Locate closest temperatures on the grid           
                  i=1
                  do while(temperature .gt. T_ku(i+1) .and. (i+1) .lt. ntemp_ku)
                     i=i+1
                  end do
                  q1=i
                  q2=i+1

               	! 2d interpolation points.		
                  z1 = sigma_ku(p1,q1)
                  z2 = sigma_ku(p2,q1)
                  z3 = sigma_ku(p1,q2)
                  z4 = sigma_ku(p2,q2)
               
                  sigma_compl(k) = interp2d(dlog10(rho_ku(p1)),dlog10(rho_ku(p2)),dlog10(T_ku(q1)),dlog10(T_ku(q2)), &
                     dlog10(z1),dlog10(z2),dlog10(z3),dlog10(z4),dlog10(density),dlog10(temperature))

               else if (density .gt. rho_ku(ndens_ku) .and. density .lt. rho_bo(1)) then
               ! Use the blending between Kumar's and Bonitz's tables
               ! 2d interpolation in three steps

              	   ! Locate closest temperatures on Kumar grid           
                  i=1
                  do while (temperature .gt. T_ku(i+1) .and. (i+1) .lt. ntemp_ku)
                     i=i+1
                  end do
                  q1=i
                  q2=i+1
                  z1 = sigma_ku(ndens_ku,q1)
                  z3 = sigma_ku(ndens_ku,q2)
                  sigma_ku_intp = interp(dlog10(T_ku(q1)),dlog10(z1),dlog10(T_ku(q2)),dlog10(z3),dlog10(temperature))
               
   		        	! Locate closest temperatures on Bonitz grid           
                  i=1
                  do while (temperature .gt. T_bo(i+1) .and. (i+1) .lt. ntemp_bo)
                     i=i+1
                  end do
                  q1=i
                  q2=i+1
               
                  z2 = sigma_bo(1,q1)
                  z4 = sigma_bo(1,q2)
                  sigma_bo_intp = interp(dlog10(T_bo(q1)),dlog10(z2),dlog10(T_bo(q2)),dlog10(z4),dlog10(temperature))
                  sigma_compl(k)=interp(dlog10(rho_ku(ndens_ku)),sigma_ku_intp,dlog10(rho_bo(1)),sigma_bo_intp,dlog10(density))

               else if (density .ge. rho_bo(1)) then

                  ! Locate closest densities on Bonitz's grid
                  i=1
                  do while (density .gt. rho_bo(i+1) .and. (i+1) .lt. ndens_bo)
                     i=i+1
                  end do
                  p1=i       
                  p2=i+1                  
                  
                  ! Locate closest temperatures on Bonitz's grid
                  i=1
                  do while (temperature .gt. T_bo(i+1) .and. (i+1) .lt. ntemp_bo)
                     i=i+1
                  end do
                  q1=i
                  q2=i+1
               
                  z1 = sigma_bo(p1,q1)
                  z2 = sigma_bo(p2,q1)
                  z3 = sigma_bo(p1,q2)
                  z4 = sigma_bo(p2,q2)

                  sigma_compl(k) = interp2d(dlog10(rho_bo(p1)),dlog10(rho_bo(p2)),dlog10(T_bo(q1)),dlog10(T_bo(q2)), &
                     dlog10(z1),dlog10(z2),dlog10(z3),dlog10(z4),dlog10(density),dlog10(temperature))

               else
                  print*, "t, rho, T: ", age, density, temperature, "Conductivity tables: CRITICAL ERROR IN DENSITY, check run_star_extras.f90. STOP"
                  stop
               endif

               ! Convert to cgs
               sigma_compl(k) = (10**(sigma_compl(k)))*sigma_SI_to_cgs

               ! Volume-average of sigma in the atmospheric level, to normalize the current
               if (k .lt. katm) &
                  sigma_atm = sigma_atm + sigma_compl(k)*(s% r(k) - s% r(k+1))*0.75d0*(s% r(k) + s% r(k+1))**2/(s% r(1)**3 - s% r(katm)**3)
            endif ! end if block for age

         end do ! end k-block for the conductivity calculation

         ! Weights entering in the average of the Joule heating normalization, to stabilize the code
         ! w_age is the weight given to the current point J^2/sigma
         ! w_age_max limits it, in order to avoid instabilities
         ! t_avg is the timescale over which the average is done (but w_age is limited, so it has no big importance)
         ! Smaller values of w_age and larger values of t_avg help in stabilizing the code
         t_avg = f_avg*age
         w_age = (age - age_prev) / t_avg
         if (w_age .gt. w_age_max) w_age = w_age_max
         if (w_age .lt. 0.d0) w_age=0.d0

         ! J is in statamp/cm^2 = (g^1/2 cm^-1/2 s^-2), 1 statamp = 10/clight A
         ! Normalization of the current as J ~ (K0/clight)*<sigma>*Bdip, where K0 is the input parameter in cm/s

         if (Bdips_prev .gt. 0) then
            Bdips = w_age*B_dip_surface + (1-w_age)*Bdips_prev
         else
            Bdips = B_dip_surface
         endif

         Bdips = B_dip_surface

         norm_J = Bdips*sigma_atm*v_avg/clight
         Bdips_prev = Bdips

         ! Calculate the current density using the conductivity profile
         if(age > age_ohmic_on) call get_j(ndim,n_pieces,norm_J,sigma_compl, &
            katm,r_dynamo/radius_p,(s% r)/radius_p,lmin,lmax,sigma_fit,jind,beta,alpha)
        
         do k = 1,ndim

            radius = s% r(k)
            rr = radius/radius_p ! Dimensionless radial coordinate.
            mass = s% m(k) ! Mass enclosed [in g].
            dmass = s% dm(k) ! Baryonic mass of cell k.
            density = s% rho(k)
            ! There is also s% entropy(k), which seems to be something different.
            entropy = exp(s% lnS(k)) ! Log of S (specific entropy).
            ! Gravitational acceleration: g = standard_cgrav*mass/radius**2.
            gravity = s% grav(k)
            pressure = s% peos(k) ! Total pressure (P = P_rad + P_gas).
            temperature = s% T(k)
            luminosity = s% L(k) ! For output purpose only.

            ! Use the interpolated conductivity to calculate the joule heating
            ! Joule heating term: Q = J^2/sigma * F(t) * window
            if(age > age_ohmic_on) then
               joule_heating = jind(k)**2d0/sigma_compl(k)

               ! Trying to regularize time oscillations
               if (k .le. nsig) then
                  if (joule_heating_prev(k) .ge. 0 .and. f_avg .gt. 0. .and. w_age .lt. 1.) then
                     joule_heating = w_age*(jind(k)**2d0/sigma_compl(k)) + (1-w_age)*joule_heating_prev(k)
                  else
                     joule_heating = jind(k)**2d0/sigma_compl(k)
                  endif
                  joule_heating_prev(k) = joule_heating                    ! Joule heating from previous timesteps

               else
                  joule_heating = jind(k)**2d0/sigma_compl(k)
               endif

               if (pcmin .gt. 0.) then
               ! Multiplying by a smoothened radial window function and by the activation function F(t)
                  window = 0.5d0*(1+dtanh((dlog(pressure)-dlog(pcmin))/s0))
                  joule_heating = joule_heating * window
               endif

               joule_heating = joule_heating * f_time_ohm

            else
               joule_heating = 0d0
            endif

            ! Specific irradiation and Joule heating rates (erg g^-1 s^-1).
            eps_irrad = s% irradiation_heat(k)
            eps_joule(k) = joule_heating/density

            ! Update the component of s with the extra heat source.
            ! s% extra_heat is of type(auto_diff_real_star_order1) so includes partials
            s% extra_heat(k) = eps_joule(k)

            ! Luminosity due to irradiation and Joule heating (erg s^-1)
            ! obtained by integrating specific heating rates in mass
            lum_irrad = lum_irrad + eps_irrad*dmass
            lum_joule = lum_joule + eps_joule(k)*dmass

            ! Joule heating above and below the main RCB.
            if(radius.gt.r_top_rcb1)then
               lum_joule_above = lum_joule_above + eps_joule(k)*dmass
            else
               lum_joule_below = lum_joule_below + eps_joule(k)*dmass
            endif

            !-----------------------------------------------------------
            ! Write radial output.
            !-----------------------------------------------------------

            ! Note: Heat equation terms.
            ! eps_heat = eps_nuc - non_nuc_neu + extra_heat + irradiation_heat
            ! Here eps_nuc is the sum for all reactions (power per unit mass),
            ! and non_nuc_neu is for the non-nuclear reaction neutrino losses.
            ! The heat terms are not summed until after this subroutine.
            ! (So, the total term s% eps_heat is not updated at this point!)
            ! Other terms: s% eps_heat(k),s% eps_nuc(k),s% non_nuc_neu(k)
            if(write_profile)then
               write(2,'(i7,25es18.9)')k,radius/r_jupiter,rr,density,pressure*p_cgs_to_bar, &
               temperature,entropy,gravity,mass,dmass, &
               s% chiRho(k),s% chiT(k),s% gamma1(k),s% gamma3(k), &
               luminosity,s% irradiation_heat(k),eps_joule(k),joule_heating, &
               sigma_thermal(k)/sigma_SI_to_cgs,sigma_compl(k)/sigma_SI_to_cgs, &
               sigma_fit(k)/sigma_SI_to_cgs,jind(k)/J_SI_to_cgs,sigma_free_e(k)/sigma_SI_to_cgs
            endif
          
         enddo
         ! End of loop over radial index.
	
         ! Close the opened files (for the radial profiles).
         if (write_profile) then
            write_profile = .FALSE.
            close(2)
            ia = ia + 1
         endif

         !--------------------------------------------------------------
         ! Output of quantitites as a function of age.
         !--------------------------------------------------------------
         if(first_call)then

            open(11,file = "LOGS/evolution.txt")
            write(11,'(a20,2f8.2)') "# Mass [Mj], Teq [K]: ",mass_p/m_jupiter,Teq
            write(11,'(a1,22a16)')"#","1.Age","2.Radius","3.L_surf","4.Teff",&
            "5.L_int","6.R_dynamo","7.B_flux","8.B_dip_surf",&
            "9.L_irr","10.L_joule","11.L_joule>RCB","12.L_joule<RCB","13.Sigma_atm",&
            "14.B_reiners","15.HeatFlux q_0","16.F_factor","17.rho_mean"

            write(11,'(a1,22a16)')"#","[yr]","[Rj]","[erg s^-1]","[K]",&
            "[erg s^-1]","[Rj]","[G]","[G]","[erg s^-1]","[erg s^-1]",&
            "[erg s^-1]","[erg s^-1]","[S/m]","[G]","[erg s^-1 cm^-2]","-","[g cm^-3]"

            open(12,file = "LOGS/evolution_conv.txt")
            write(12,'(a1,12a16)')"#","1.Age","2.Radius",&
            "3.R_bot main RCB","4.P_bot main RCB","5.R_top main RCB","6.P_top main RCB", &
            "7.R_bot Sh.RCB","8.P_bot Sh.RCB","9.R_top Sh.RCB","10.P_top.Sh.RCB", &
            "11.Surf P","12.Sh.RCB?"

            write(12,'(a1,12a16)')"#","[yr]","[Rj]","[Rj]","[bar]", &
            "[Rj]","[bar]","[Rj]","[bar]","[Rj]","[bar]","[bar]","[1=yes]"


            open(13,file = "LOGS/evolution_beta.txt")
            open(14,file = "LOGS/evolution_alpha.txt")

         else
            open(11,file = "LOGS/evolution.txt",position = "append")
            write(11,'(22es16.7)') age,radius_p/r_jupiter, &
            luminosity_p,s% Teff,luminosity_int,r_dynamo/r_jupiter,&
            B_flux,B_dip_surface,lum_irrad,lum_joule,lum_joule_above,&
            lum_joule_below,sigma_atm/sigma_SI_to_cgs,B_reiners,q_0,F_factor,rho_mean

            open(12,file = "LOGS/evolution_conv.txt",position = "append")
            write(12,'(11es16.7,i8)') age,radius_p/r_jupiter, &
            (s% r(ndim))/r_jupiter,(s% peos(ndim))*p_cgs_to_bar, &
            r_top_rcb1/r_jupiter,p_top_rcb1*p_cgs_to_bar, &
            r_bot_rcb2/r_jupiter,p_bot_rcb2*p_cgs_to_bar, &
            r_top_rcb2/r_jupiter,p_top_rcb2*p_cgs_to_bar, &
            p_surface*p_cgs_to_bar,conv_sup

            if (age .gt. age_ohmic_on) then
               open(13,file = "LOGS/evolution_beta.txt",position = "append")
               open(14,file = "LOGS/evolution_alpha.txt",position = "append")
               write(13,'(1es16.7,100f9.4)') age,beta(1:n_pieces)
               write(14,'(1es16.7,100f9.4)') age,alpha(1:n_pieces,1),alpha(1:n_pieces,lmax-lmin+1)
            endif

         endif

         close(11)
         close(12)

         if (first_call) first_call = .FALSE.

         age_prev = age
         iage = iage + 1

         return

      end subroutine ohmic_heating

end module run_star_extras
