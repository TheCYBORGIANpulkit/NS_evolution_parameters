! ***********************************************************************
!
!   Copyright (C) 2012-2019  Bill Paxton & The MESA Team
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
!   gnu library general public license for more details.o;
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! *********************************************************************** 
      module run_binary_extras 

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      use math_lib
      !use binary_wind
      
      implicit none
      
      contains
      
     
      
      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         write(*,*) 'hello from extra_binary_controls'

         ! Set these function pointers to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
         b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% other_adjust_mdots => my_adjust_mdots

         b% other_mdot_edd => eval_mdot_edd_routine
         !b% other_jdot_ml  => default_jdot_ml  
         b% other_jdot_ml  => rcl_jdot_ml        ! take the SAM of CL when radio turn on (disk instable) Jia
         !b% other_jdot_mb  => cz_jdot_mb          ! considering convective envelop 
         !b% other_mdot_edd => my_mdot_edd

         b% extras_binary_startup=> extras_binary_startup
         b% extras_binary_start_step=> extras_binary_start_step
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
         ! b% warn_binary_extra =.false.
         
      end subroutine extras_binary_controls

!      subroutine my_mdot_edd(binary_id, mdot_edd, mdot_edd_eta, ierr)
!         use const_def, only: dp
!         integer, intent(in) :: binary_id
!         real(dp), intent(out) :: mdot_edd
!         real(dp), intent(out) :: mdot_edd_eta
!         integer, intent(out) :: ierr
!         type (binary_info), pointer :: b
!         ierr = 0
!         call binary_ptr(binary_id, b, ierr)
!         if (ierr /= 0) then
!            write(*,*) 'failed in binary_ptr'
!            return
!         end if
!         mdot_edd = 0d0 
!         mdot_edd_eta = 0d0 
!      end subroutine my_mdot_edd
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!written by Xiakun and Li and modified by me to accomodate for the gravitational and baryonic mass conversion factor

      subroutine eval_mdot_edd_routine(binary_id , mdot_edd , mdot_edd_eta, ierr)
         integer , intent (in) :: binary_id
         real ( dp ) , intent (out) :: mdot_edd
         real(dp), intent(out) :: mdot_edd_eta
         integer , intent ( out ) :: ierr
         real :: mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
                 cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
                 ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
                 rlo_start_age,rlo_age,radio_start_age,radio_age,vesc_sur_d,eva_wind,&
                 mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,hoc_accretor,hoc_donor
         !lp_chen: Chen,H-L et al.,2013 ApJ

         common /binary/   mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
                           cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
                           ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
                           rlo_start_age,rlo_age,radio_start_age,radio_age,vesc_sur_d,eva_wind,&
                           mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,&
                           hoc_accretor,hoc_donor
         type(binary_info), pointer :: b
         include 'formats.inc'
         ierr = 0
         
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         if (b% point_mass_i == 0) then
            if (b% limit_retention_by_mdot_edd) then
               write(*,*) "Default mdot_edd calculation cannot be used when evolving both stars"
               write(*,*) "Maybe you want to set limit_retention_by_mdot_edd=.false. in binary_controls?"
               write(*,*) "Setting mdot_edd to zero"
            end if
            mdot_edd = 0
            mdot_edd_eta = 0
            return
         end if
         
!         if (is_bad(mdot_edd_eta) .or. mdot_edd_eta<0) then
!            write(*,*) "ERROR while computing mdot_edd_eta"
!            ierr = -1
!            write(*,*) "mdot_edd_eta, b% m(b% a_i), b% eq_initial_bh_mass", &
!               mdot_edd_eta, b% m(b% a_i), b% eq_initial_bh_mass
!            stop
!         end if

!         if (is_bad(mdot_edd) .or. mdot_edd<0) then
!            write(*,*) "ERROR while computing mdot_edd"
!            ierr = -1
!            write(*,*) "mdot_edd, b% m(b% a_i), b% s_donor% opacity(1), b% s_donor% surface_h1", &
!               mdot_edd, b% m(b% a_i), b% s_donor% opacity(1), b% s_donor% surface_h1
!            stop
!         end if
         
         if (b% use_this_for_mdot_edd_eta > 0) then
            !mdot_edd_eta = b% use_this_for_mdot_edd_eta
         else
            ! eg., eq. (6) of Podsiadlowski, Rappaport & Han 2003, MNRAS, 341, 385
            !mdot_edd_eta = 1d0 - sqrt(1d0 - pow2(min(b% m(b% a_i),sqrt(6d0)*b% eq_initial_bh_mass)/(3d0*b% eq_initial_bh_mass)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           modification to the mdot_edd_eta, given by eta = 1 - dM_g/dM_b
            mdot_edd_eta = 0.1 * ((b% m(b% a_i)/(1.989d33)) + 0.065*(b% m(b% a_i)/(1.989d33))*(b% m(b% a_i)/(1.989d33)))
            b% xtra(1) = mdot_edd_eta
            write(*,*) "eta" , "..................", mdot_edd_eta , "........" , b% xtra(1)
            !write(*,*) mdot_edd_eta , " .,. " , (b% m(b% a_i)/(1.989d33))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         end if
         !mdot_edd_eta = 0

         if (b% use_this_for_mdot_edd > 0) then
            mdot_edd = b% use_this_for_mdot_edd*(Msun/secyer)
         else
            ! eg., eq. (9) of Podsiadlowski, Rappaport & Han 2003, MNRAS, 341, 385
            if (.not. b% use_es_opacity_for_mdot_edd) then
               mdot_edd = 4d0*pi*standard_cgrav*b% m(b% a_i) &
                  /(clight*b% s_donor% opacity(1)*mdot_edd_eta)
            else
               mdot_edd = 4d0*pi*standard_cgrav*b% m(b% a_i)&
                  /(clight*0.2d0*(1d0+b% s_donor% surface_h1)*mdot_edd_eta)
            end if
         end if

              
         mdot_edd = 3.6d-8*(b% m(2)/(1.4*Msun))*(1.0d5*(clight**2.0)/(6.67428d-8*b% m(2)))*&
                   (1.7/(1.0+b% s1% X(1)))*Msun/secyer

         if (mdot_act .lt. mdot_critical_dubus) then
            mdot_edd = 0.01 * mdot_edd
         end if
         
      end subroutine eval_mdot_edd_routine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine eval_wind_xfer_fractions(binary_id, ierr)
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in binary_ptr'
         return
      end if
      
      ! for the primary
      if (b% point_mass_i /= 1) then
         if (.not. b% do_wind_mass_transfer_1 .or. b% model_twins_flag) then
            b% wind_xfer_fraction(1) = 0d0
         else if(.not. b% use_other_binary_wind_transfer) then
            call Bondi_Hoyle_wind_transfer(b% binary_id, 1, ierr)
            if (ierr /=0) then
               write(*,*) "Error in Bondi_Hoyle_wind_transfer(b% binary_id, 1, ierr)"
               return
            end if
         else
            call b% other_binary_wind_transfer(b% binary_id, 1, ierr)
            if (ierr /=0) then
               write(*,*) "Error in other_binary_wind_transfer(b% binary_id, 1, ierr)"
               return
            end if
         end if
      end if
      
      ! check if secondary needs wind transfer
      if (b% point_mass_i /= 2) then
         if (.not. b% do_wind_mass_transfer_2) then
            b% wind_xfer_fraction(2) = 0d0
         else if(.not. b% use_other_binary_wind_transfer) then
            call Bondi_Hoyle_wind_transfer(b% binary_id, 2, ierr)
            if (ierr /=0) then
               write(*,*) "Error in Bondi_Hoyle_wind_transfer(b% binary_id, 2, ierr)"
               return
            end if
         else
            call b% other_binary_wind_transfer(b% binary_id, 2, ierr)
            if (ierr /=0) then
               write(*,*) "Error in other_binary_wind_transfer(b% binary_id, 2, ierr)"
               return
            end if
         end if
      end if
         
   end subroutine eval_wind_xfer_fractions

   subroutine Bondi_Hoyle_wind_transfer(binary_id, s_i, ierr)
      integer, intent(in) :: binary_id, s_i ! s_i is index of the wind mass losing star
      integer, intent(out) :: ierr

      ! wind transfer fraction based on Bondi-Hoyle mechanism as described in
      ! Hurley et al. 2002, MNRAS, 329, 897-928

      type(binary_info), pointer :: b
      type (star_info), pointer :: s
      real(dp) :: v_orb, v_wind, b_BH
      real(dp) :: alpha  ! Bondi-Hoyle alpha for that star
      real(dp) :: max_xfer  ! Maximum transfer fraction

      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in binary_ptr'
         return
      end if
      
      if (s_i == 1) then
         s => b% s1
         alpha = b% wind_BH_alpha_1
         max_xfer = b% max_wind_transfer_fraction_1
      else
         s => b% s2
         alpha = b% wind_BH_alpha_2
         max_xfer = b% max_wind_transfer_fraction_2
      end if
      
      ! orbital speed Hurley et al 2002 eq. 8
      v_orb = sqrt(standard_cgrav * b% m(s_i) / b% separation) !cm/s
      
      ! windspeed from Hurley et al 2002 eq. 9
      v_wind = sqrt( 2d0 / 8d0 *  standard_cgrav * b% m(s_i) / b% r(s_i) )
      
      ! Bondi-Hoyle transfer fraction Hurley et al. 2002 eq. 6
      b% wind_xfer_fraction(s_i) = alpha / pow2(b% separation) /&
                  (2d0 * sqrt(1d0 - pow2(b% eccentricity))) *&
                  pow2(standard_cgrav * b% m(3-s_i) / pow2(v_wind)) *&
                  pow(1d0 + pow2(v_orb/v_wind),-1.5d0)
                  
      ! limit to provided maximum
      b% wind_xfer_fraction(s_i) = min(max_xfer, b% wind_xfer_fraction(s_i))
      
   end subroutine Bondi_Hoyle_wind_transfer
   
   subroutine Tout_enhance_wind(b, s)
      type (binary_info), pointer :: b
      type (star_info), pointer :: s

      ! Tidaly enhance wind mass loss as described by
      ! Tout & Eggleton 1988,MNRAS,231,823 (eq. 2)
      real(dp) :: B_wind  ! enhancement parameter, B in eq. 2
      integer :: i, s_i
      real(dp) :: dm
      real(dp), DIMENSION(b% anomaly_steps):: rl_d, r_rl, mdot

      if (s% id == b% s1% id) then
         if (.not. b% do_enhance_wind_1) return
         B_wind = b% tout_B_wind_1
         s_i = 1
      else
         if (.not. b% do_enhance_wind_2) return
         B_wind = b% tout_B_wind_2
         s_i = 2
      end if
      
      do i = 1,b% anomaly_steps !limit radius / roche lobe
         ! phase dependent roche lobe radius
         rl_d(i) = (1d0-pow2(b%eccentricity)) / (1+b%eccentricity*cos(b% theta_co(i))) * b% rl(s_i)
         r_rl(i) = min(pow6(b% r(s_i) / rl_d(i)), pow6(0.5d0))
      end do
      
      ! actual enhancement
      mdot = s% mstar_dot * (1 + B_wind * r_rl)
      
      dm = 0d0
      do i = 2,b% anomaly_steps ! trapezoidal integration
         dm = dm + 0.5d0 * (mdot(i-1) + mdot(i)) * (b% time_co(i) - b% time_co(i-1)) 
      end do
      
      ! remember mass-loss is negative!
      !b% mdot_wind_theta = b% mdot_wind_theta + mdot ! store theta dependance for edot
      s% mstar_dot = dm ! return enhanced wind mass loss
      
   end subroutine Tout_enhance_wind

      subroutine eval_accreted_material_j(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b
         real(dp) :: qratio, min_r
         logical, parameter :: dbg = .false.
         include 'formats.inc'

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         qratio = b% m(b% a_i) / b% m(b% d_i)
         qratio = min(max(qratio,0.0667d0),15d0)
         min_r = 0.0425d0*b% separation*pow(qratio+qratio*qratio, 0.25d0)

         !TODO: MUST USE EQUATORIAL RADIUS
         if (dbg) write(*,*) "radius, impact_radius, separation: ", &
             b% r(b% a_i), min_r/rsun, b% separation/rsun
         if (b% r(b% a_i) < min_r) then
            b% accretion_mode = 2
            b% s_accretor% accreted_material_j = &
               sqrt(standard_cgrav * b% m(b% a_i) * b% r(b% a_i)) 
         else
            b% accretion_mode = 1
            b% s_accretor% accreted_material_j = &
               sqrt(standard_cgrav * b% m(b% a_i) * 1.7d0*min_r)
         end if
         b% acc_am_div_kep_am = b% s_accretor% accreted_material_j / &
             sqrt(standard_cgrav * b% m(b% a_i) * b% r(b% a_i))

          !TODO: when using wind mass transfer donor star can end up
          ! with positive mdot, need to properly set jdot in that case

      end subroutine eval_accreted_material_j
      
      subroutine eval_mdot_edd(binary_id, mdot_edd, mdot_edd_eta, ierr)
         use utils_lib, only: is_bad

         integer, intent(in) :: binary_id
         real(dp), intent(out) :: mdot_edd ! eddington accretion rate
         real(dp), intent(out) :: mdot_edd_eta ! fraction of rest mass energy released as radiation
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b
         include 'formats.inc'

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         if (b% point_mass_i == 0) then
            if (b% limit_retention_by_mdot_edd) then
               write(*,*) "Default mdot_edd calculation cannot be used when evolving both stars"
               write(*,*) "Maybe you want to set limit_retention_by_mdot_edd=.false. in binary_controls?"
               write(*,*) "Setting mdot_edd to zero"
            end if
            mdot_edd = 0
            mdot_edd_eta = 0
            return
         end if

         if (b% use_this_for_mdot_edd_eta > 0) then
            mdot_edd_eta = b% use_this_for_mdot_edd_eta
         else
            ! eg., eq. (6) of Podsiadlowski, Rappaport & Han 2003, MNRAS, 341, 385
            mdot_edd_eta = 1d0 &
                 - sqrt(1d0 - pow2(min(b% m(b% a_i),sqrt(6d0)*b% eq_initial_bh_mass)/(3d0*b% eq_initial_bh_mass)))
         end if

         if (b% use_this_for_mdot_edd > 0) then
            mdot_edd = b% use_this_for_mdot_edd*(Msun/secyer)
         else
            ! eg., eq. (9) of Podsiadlowski, Rappaport & Han 2003, MNRAS, 341, 385
            if (.not. b% use_es_opacity_for_mdot_edd) then
               mdot_edd = 4d0*pi*standard_cgrav*b% m(b% a_i) &
                  /(clight*b% s_donor% opacity(1)*mdot_edd_eta)
            else
               mdot_edd = 4d0*pi*standard_cgrav*b% m(b% a_i)&
                  /(clight*0.2d0*(1d0+b% s_donor% surface_h1)*mdot_edd_eta)
            end if
         end if

         if (is_bad(mdot_edd_eta) .or. mdot_edd_eta<0) then
            write(*,*) "ERROR while computing mdot_edd_eta"
            ierr = -1
            write(*,*) "mdot_edd_eta, b% m(b% a_i), b% eq_initial_bh_mass", &
               mdot_edd_eta, b% m(b% a_i), b% eq_initial_bh_mass
            stop
         end if

         if (is_bad(mdot_edd) .or. mdot_edd<0) then
            write(*,*) "ERROR while computing mdot_edd"
            ierr = -1
            write(*,*) "mdot_edd, b% m(b% a_i), b% s_donor% opacity(1), b% s_donor% surface_h1", &
               mdot_edd, b% m(b% a_i), b% s_donor% opacity(1), b% s_donor% surface_h1
            stop
         end if

      end subroutine eval_mdot_edd
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      subroutine my_adjust_mdots(binary_id, ierr)
         use const_def, only: dp
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: fixed_xfer_fraction, actual_mtransfer_rate, mdot_b
         actual_mtransfer_rate = 0d0
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         !calculation of transfer efficiency
         b% fixed_xfer_fraction = 1 - b% mass_transfer_alpha - b% mass_transfer_beta - &
            b% mass_transfer_delta
         ! calculation of eddington's rate of transfer'
         if (.not. b% use_other_mdot_edd) then
            call eval_mdot_edd(b% binary_id, b% mdot_edd, b% mdot_edd_eta, ierr)
         else
            call b% other_mdot_edd(b% binary_id, b% mdot_edd, b% mdot_edd_eta, ierr)
         end if

         ! Add tidal enhancement of wind
         call Tout_enhance_wind(b, b% s_donor)
         if (b% point_mass_i == 0) then
            ! do not repeat if using the implicit wind
            if (.not. (b% num_tries >0 .and. b% s_accretor% was_in_implicit_wind_limit)) &
               call Tout_enhance_wind(b, b% s_accretor)
         end if

         ! solve wind mass transfer
         ! b% mdot_wind_transfer(b% d_i) is a negative number that gives the
         ! amount of mass transferred by unit time from the donor to the
         ! accretor.
         call eval_wind_xfer_fractions(b% binary_id, ierr)
         if (ierr/=0) then
            write(*,*) "Error in eval_wind_xfer_fractions"
            return
         end if
         b% mdot_wind_transfer(b% d_i) = b% s_donor% mstar_dot * &
            b% wind_xfer_fraction(b% d_i)
         if (b% point_mass_i == 0) then
            b% mdot_wind_transfer(b% a_i) = b% s_accretor% mstar_dot * &
               b% wind_xfer_fraction(b% a_i)
         else
            b% mdot_wind_transfer(b% a_i) = 0d0
         end if

         ! Set mdot for the donor
         b% s_donor% mstar_dot = b% s_donor% mstar_dot + b% mtransfer_rate - &
            b% mdot_wind_transfer(b% a_i)

         ! Set mdot for the accretor
         if (b% point_mass_i == 0 .and. .not. b% CE_flag) then
            ! do not repeat if using the implicit wind
            if (.not. (b% num_tries >0 .and. b% s_accretor% was_in_implicit_wind_limit)) then
               b% accretion_mode = 0
               b% acc_am_div_kep_am = 0.0d0
               b% s_accretor% mstar_dot = b% s_accretor% mstar_dot - &
                  b% mtransfer_rate*b% fixed_xfer_fraction - b% mdot_wind_transfer(b% d_i)

               !set angular momentum accretion as described in A.3.3 of de Mink et al. 2013
               if (b% do_j_accretion) then
                  if (.not. b% use_other_accreted_material_j) then
                     call eval_accreted_material_j(b% binary_id, ierr)
                  else
                     call b% other_accreted_material_j(b% binary_id, ierr)
                  end if
                  if (ierr /= 0) then
                     write(*,*) 'error in accreted_material_j'
                     return
                  end if
               end if
            end if

            b% accretion_luminosity = 0d0 !only set for point mass

         else if (.not. b% CE_flag) then
            ! accretor is a point mass
            if (.not. b% model_twins_flag) then
               !combine wind and RLOF mass transfer
               actual_mtransfer_rate = b% mtransfer_rate*b% fixed_xfer_fraction+b% mdot_wind_transfer(b% d_i) !defined negative
               b% component_mdot(b% a_i) = -actual_mtransfer_rate
               ! restrict accretion to the Eddington limit
               if (b% limit_retention_by_mdot_edd .and. b% component_mdot(b% a_i) > b% mdot_edd) then
                  b% component_mdot(b% a_i) = b% mdot_edd ! remove all accretion above the edd limit
               end if
               b% accretion_luminosity = &
                  b% mdot_edd_eta*b% component_mdot(b% a_i)*clight*clight
               ! remove rest mass radiated away
               if (b% use_radiation_corrected_transfer_rate) then
                  b% component_mdot(b% a_i) = (1 - b% mdot_edd_eta) * b% component_mdot(b% a_i)
               end if
               ! the following variable stores baryonic rate of accretion in Msun/year
               b% xtra(2) = b% component_mdot(b% a_i)*(3.154e+7/1.989e+33)
               write(*,*) "mdot :::" , ".....", b% xtra(2) , "............" , b% component_mdot(b% a_i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!modification in the accretor mdot induced
! As written in the Equation 19 9a) and (b) of the paper Cipolleta et al.            
               b% component_mdot(b% a_i) = b% component_mdot(b% a_i) - 0.1*((b% m(b% a_i)/(1.989d33)) + 0.065*(b% m(b% a_i)/(1.989d33))*(b% m(b% a_i)/(1.989d33)))*b% component_mdot(b% a_i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              

            end if
         else
            !doing CE, just be sure to set mdot for a point mass to zero
            b % accretion_luminosity = 0d0
            if (b% point_mass_i /= 0) then
               b% component_mdot(b% a_i) = 0d0
            end if
         end if

         ! mdot_system_transfer is mass lost from the vicinity of each star
         ! due to inefficient rlof mass transfer, mdot_system_cct is mass lost
         ! from a circumbinary coplanar toroid.
         if (b% mtransfer_rate+b% mdot_wind_transfer(b% d_i) >= 0 .or. b% CE_flag) then
            b% mdot_system_transfer(b% d_i) = 0d0
            b% mdot_system_transfer(b% a_i) = 0d0
            b% mdot_system_cct = 0d0
         else 
            b% mdot_system_transfer(b% d_i) = b% mtransfer_rate * b% mass_transfer_alpha
            b% mdot_system_cct = b% mtransfer_rate * b% mass_transfer_delta
            if (b% point_mass_i == 0 .or. b% model_twins_flag) then
               b% mdot_system_transfer(b% a_i) = b% mtransfer_rate * b% mass_transfer_beta
            else
               ! do not compute mass lost from the vicinity using just mass_transfer_beta, as
               ! mass transfer can be stopped also by going past the Eddington limit
               b% mdot_system_transfer(b% a_i) = (actual_mtransfer_rate + b% component_mdot(b% a_i)) &
                  + b% mtransfer_rate * b% mass_transfer_beta
            end if
         end if
         
      end subroutine my_adjust_mdots

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         
!          b% s1% job% warn_run_star_extras = .false.
          extras_binary_startup = keep_going
      end function  extras_binary_startup
      
      integer function extras_binary_start_step(binary_id,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr

         extras_binary_start_step = keep_going
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
      
      end function  extras_binary_start_step
      
      !Return either keep_going, retry or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_check_model = keep_going
        
      end function extras_binary_check_model
      
      
      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_finish_step = keep_going
         
      end function extras_binary_finish_step
      
      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if 
      end subroutine extras_binary_after_evolve   
      
      
      subroutine rcl_jdot_ml(binary_id, ierr)
         !rlo_start_age: mdot_act > 1.0e-50
         !rlo_age : age - rlo_start_age
         !mdot_act : rlo + accreted wind
         !radio_start_age : when accreted mass > 0.05 (Antoniadis 2012) and rlo > mdot_critical (Dubus 1999)
         !radio_age : age - radio_start_age
         !
         INTEGER*4 getcwd, status
         character(LEN=100) dirname,dirname_postfix,dirname_prefix
         character n_dirname
         integer n_string
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         integer model_number_last,i,store_read
         integer recyelled
         real :: L1_to_acc, L1_to_cm, w_Kopal


         real :: mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
                 cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
                 ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
                 rlo_start_age,rlo_age,radio_start_age,radio_age,vesc_sur_d,eva_wind,&
                 mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,hoc_accretor,hoc_donor
         !lp_chen: Chen,H-L et al.,2013 ApJ
         common /dirname/ dirname
         common /binary/   mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
                           cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
                           ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
                           rlo_start_age,rlo_age,radio_start_age,radio_age,vesc_sur_d,eva_wind,&
                           mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,&
                           hoc_accretor,hoc_donor,model_number_last,recyelled

         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         if (b% s_donor% model_number .eq. 1)   then
		rlo_start_age = 0.0
		rlo_age = 0.0
		radio_start_age = 0.0
		radio_age = 0.0
		maccretor_i = (b% m(2)/Msun)
		eva_wind = 0.0
            mtransfer_rate_dot = 0.0
            mtransfer_rate_dot_div_mtransfer = 0.0
            recyelled = 0
         end if

         accreted_mass = (b% m(2)/Msun) - maccretor_i
         if (accreted_mass .ge. 0.1) recyelled = 1
	   mdot_critical_dubus = 3.2D-9*(((b% m(2)/Msun)/1.4)**0.5)*&
			(((b% m(1)/Msun)/1.0)**(-0.2))*((b% period/86400)**1.4)    !disk instability (Dubus et al. 1999)
         mdot_act = abs((b% mtransfer_rate) + (b% mdot_wind_transfer(b% d_i)))*secyer/Msun              
         ! rlo + wind_acc
         rlo = abs(b% mtransfer_rate)*secyer/Msun
         mdonor_wind = (b% mdot_system_transfer(b% d_i) + b% mdot_system_wind(b% d_i))*secyer/Msun
         maccretor_wind = (b% mdot_system_transfer(b% a_i) + b% mdot_system_wind(b% a_i))*secyer/Msun
         cbd_wind = b% mdot_system_cct*secyer/Msun
         mdot_edd_here = 3.6d-8*(b% m(2)/(1.4*Msun))*(1.0d5*(clight**2.0)/(6.67428d-8*b% m(2)))*&
                   (1.7/(1.0+b% s1% X(1)))
         vesc_sur_d = (2*6.672D-8*b% m(1)/b% r(1))**0.5

         if (rlo .ge. mdot_critical_dubus) radio_start_age = 0.0
         if ((rlo .gt. 1.0e-50) .and. (rlo_start_age .eq. 0.0)) rlo_start_age = b% s_donor% time/secyer
         if (rlo_start_age .gt. 0.0) rlo_age = b% s_donor% time/secyer - rlo_start_age
         if (radio_start_age .gt. 0.0) then
            radio_age = ( b% s_donor% time/secyer - radio_start_age)
            eva_wind = (b% s_donor% x_ctrl(1)*(1.0/(2.0*(vesc_sur_d**2)))*lp_chen*(b% r(1)/b% separation)**2)*secyer/Msun
         else
            radio_age = 0.0
            eva_wind = 0.0
         end if


         if (((radio_start_age .eq. 0) .and. (recyelled .eq. 1)   &
         .and. (rlo .le. mdot_critical_dubus)) ) then
         !(mtransfer_rate_dot_div_mtransfer .LE. 1E-7 .or. rlo .eq. 0)
            radio_start_age = b% s_donor% time/secyer

         end if



         !!!!!!!!!!!!!!!!!!!!!!!!
         ps_chen = 2.0*pi/(8.11156D11/(1.5D17+(radio_age)*3.15E7)**0.5)    !Chen et. al 2013 APJ 775:27 P0=3ms B~10^8G
         psdot_chen = (6.9813D-15/(2.0*3.14))*ps_chen**2.0
         lp_chen = 4.0*(pi**2)*(1.0D45)*psdot_chen/(ps_chen**3.0)
         rcl_chen = 2.99792458D10*ps_chen/(2.0*pi)

         ! wind in units of Msun/year  (value is >= 0)

         if (radio_age .eq. 0 .or. b% s_donor% x_ctrl(1) .eq. 0) then !no evaporation : disk stable + accreted_mass > 0.05

         	!mass lost from vicinity of donor
         	b% jdot_ml = (b% mdot_system_transfer(b% d_i) + b% mdot_system_wind(b% d_i))*&
             		(b% m(b% a_i)/(b% m(b% a_i)+b% m(b% d_i))*b% separation)**2*2*pi/b% period
         	!mass lost from vicinity of accretor
         	b% jdot_ml = b% jdot_ml + (b% mdot_system_transfer(b% a_i) + b% mdot_system_wind(b% a_i))*&
             		(b% m(b% d_i)/(b% m(b% a_i)+b% m(b% d_i))*b% separation)**2*2*pi/b% period
         	!mass lost from circumbinary coplanar toroid
         	b% jdot_ml = b% jdot_ml + b% mdot_system_cct * b% mass_transfer_gamma * &
             		sqrt(b% s_donor% cgrav(1) * (b% m(1) + b% m(2)) * b% separation)

         else   !evaporation start : no accretion b% m(2) = conste, take CL SAM
                b% m(2) = b% m_old(2)
         	                                       !mass lost from vicinity of donor
         	b% jdot_ml = (b% mdot_system_transfer(b% d_i) + b% mdot_system_wind(b% d_i))*&
             		(b% m(b% a_i)/(b% m(b% a_i)+b% m(b% d_i))*b% separation)**2*2*pi/b% period

            if (b% m(b% a_i)/b% m(b% d_i) .le. 10.0) then
            L1_to_acc = 0.5 + 0.227*log10(b% m(b% a_i)/b% m(b% d_i))               !distance from L1 to accretor in unit of separation
            else
            w_Kopal = (1.0/(3.0*(1+(b% m(b% a_i)/b% m(b% d_i)))))**(1.0/3.0)       !Beer 2007;Kopal1959
            L1_to_acc = 1.0 - w_Kopal + ((w_Kopal**2.0)/3.0) + ((w_Kopal**3.0)/9.0)           !distance from L1 to accretor in unit of separation
            endif
            L1_to_cm = abs(b% m(b% d_i)/(b% m(b% a_i)+b% m(b% d_i)) - L1_to_acc)   !distance from L1 to center of mass in unit of separation

         	!mass lost from L1 point
         	!b% jdot_ml = (b% mdot_system_transfer(b% d_i) + b% mdot_system_wind(b% d_i))*&
            !		(L1_to_cm*b% separation)**2*2*pi/b% period


         	!evaporation wind takes no AML
         	!b% jdot_ml = 0.0
         	!mass lost from vicinity of CL of accretor
         	b% jdot_ml = b% jdot_ml + (b% mdot_system_transfer(b% a_i) + b% mdot_system_wind(b% a_i))*&
             		((b% m(b% d_i)/(b% m(b% a_i)+b% m(b% d_i))*b% separation) - rcl_chen)**2*2*pi/b% period
         	!mass lost from circumbinary coplanar toroid
         	b% jdot_ml = b% jdot_ml + b% mdot_system_cct * b% mass_transfer_gamma * &
             		sqrt(b% s_donor% cgrav(1) * (b% m(1) + b% m(2)) * b% separation)
         end if

         mdonor_dot = ((b% m(1) - b% m_old(1))/(b% s_donor% time_step))/Msun
         maccretor_dot = ((b% m(2) - b% m_old(2))/(b% s_donor% time_step))/Msun
         !maccretor_dot = maccretor_dot*100
         mtransfer_rate_dot = abs(b% mtransfer_rate - b% mtransfer_rate_old)/b% s_donor% time_step

         if (b% mtransfer_rate .ne. 0) then
         	mtransfer_rate_dot_div_mtransfer = mtransfer_rate_dot/abs(b% mtransfer_rate)
         end if

         !if (b% s_donor% model_number .eq. 1) then
         !open (unit=10, file='fort.10')	
         !write (10,"(A7,100(A20))") 'No.','dt(yr)','age(yr)','Md(Msun)','Porb(d)','Sep','omega',&
         ! 	'Ma(Msun)','Md_dot','Ma_dot','Mdot_edd','Md_wind','Ma_wind','CBD_wind','Rlo(Msun/yr)',&
         ! 	'rlo_fraction','wind_fra(d_to_a)','wind_fra(a_to_d)','R_donor','Rl_donor',&
         ! 	'Ts_donor','L_donor','logP_sur_d','logrho_sur_d','logg_sur_d','vesc_sur_d','eva_wind',&
         !	'Jdot','Jdot_mb','Jdot_gr','Jdot_ml', 'jdot_ls', 'jdot_missing_wind', 'extra_jdot','orbital_am','rlo_dot','rlo_dot/Mdot'
         !close (10)
         !end if

         !open (unit=10, file='fort.10', access='append')	
         !write (10,"(I7,100(E20.4E3))") b% s_donor% model_number,&
	   !		b% s_donor% time_step/secyer,b% s_donor% time/secyer,b% m(1)/Msun,&
		!	b% period/86400,b% separation,b% s_donor% omega(1), b% m(2)/Msun,&
		!	mdonor_dot,maccretor_dot,mdot_edd_here,mdonor_wind,maccretor_wind,&
		!	cbd_wind,rlo,b% xfer_fraction,&
		!	b% wind_xfer_fraction(b% d_i),b% wind_xfer_fraction(b% a_i),&
		!	b% r(1),b% rl(1),b% s_donor% Teff,b% s_donor% L_phot,&
		!	b% s_donor% log_surface_pressure,b% s_donor% log_surface_density,&
		!	b% s_donor% log_surface_gravity,vesc_sur_d,eva_wind,&
		!	b% jdot, b% jdot_mb, b% jdot_gr, b% jdot_ml, b% jdot_ls, &
		!	b% jdot_missing_wind, b% extra_jdot, b% angular_momentum_j,&
		!	mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer
			
         !close (10)


         if (b% s_donor% model_number .eq. 1) then
         	model_number_last = 0
            status = getcwd(dirname)
         	do  n_string=1,len(dirname)-3
                  n_dirname = dirname(len(dirname)-n_string-2:len(dirname)-n_string-1)
         		if (n_dirname .eq. '/' ) dirname=dirname(len(dirname)-1-n_string:len(dirname))
            end do
            dirname_postfix = '.59'
            dirname = trim(dirname)//trim(dirname_postfix)
         	open (unit=59, file=dirname)
         	close (59,status='delete')
         end if





         if (b% s_donor% model_number .gt. model_number_last) then

         open (unit=59, file=dirname, access='append')	
         write (59,"(I7,2I5,100(E20.4E3))") &
            b% s_donor% model_number, b% d_i, b% accretion_mode,           &
            b% s_donor% star_age,                                          &
            radio_start_age, b% m(1)/Msun, b% m(2)/Msun, mdonor_dot,       &
            maccretor_dot, mdot_edd_here, rlo, b% wind_xfer_fraction,           &
            mdonor_wind, eva_wind, maccretor_wind, mdot_act,               &
            mdot_critical_dubus, b% period/86400, b% separation/Rsun,      &
            radio_age, ps_chen, lp_chen, b% eccentricity, b% r(1)/Rsun,    &
            b% r(2)/Rsun, b% rl(1)/Rsun, b% rl(2)/Rsun,                    &
            b% angular_momentum_j, b% s1% total_angular_momentum,          &
            b% jdot, b% jdot_mb, b% jdot_gr, b% jdot_ml, b% jdot_ls,       &
            b% jdot_missing_wind, b% extra_jdot,                           &
            b% s_donor% L_phot, b% s_donor% Teff,                          &
            b% s_donor% log_surface_density, b% s_donor% he_core_mass,     &
            b% s_donor% c_core_mass, b% s_donor% o_core_mass,              &
            b% s_donor% si_core_mass, b% s_donor% fe_core_mass,            &
            b% s_donor% power_h_burn, b% s_donor% power_he_burn,           &
            b% s_donor% power_c_burn

         else

         open (unit=59, file=dirname, position='REWIND')	
         do i=1,model_number_last-1
         read (59,*) store_read
         end do

         write (59,"(I7,2I5,100(E20.4E3))") &
            b% s_donor% model_number, b% d_i, b% accretion_mode,           &
            b% s_donor% star_age,                                          &
            radio_start_age, b% m(1)/Msun, b% m(2)/Msun, mdonor_dot,       &
            maccretor_dot, mdot_edd_here, rlo, b% wind_xfer_fraction,           &
            mdonor_wind, eva_wind, maccretor_wind, mdot_act,               &
            mdot_critical_dubus, b% period/86400, b% separation/Rsun,      &
            radio_age, ps_chen, lp_chen, b% eccentricity, b% r(1)/Rsun,    &
            b% r(2)/Rsun, b% rl(1)/Rsun, b% rl(2)/Rsun,                    &
            b% angular_momentum_j, b% s1% total_angular_momentum,          &
            b% jdot, b% jdot_mb, b% jdot_gr, b% jdot_ml, b% jdot_ls,       &
            b% jdot_missing_wind, b% extra_jdot,                           &
            b% s_donor% L_phot, b% s_donor% Teff,                          &
            b% s_donor% log_surface_density, b% s_donor% he_core_mass,     &
            b% s_donor% c_core_mass, b% s_donor% o_core_mass,              &
            b% s_donor% si_core_mass, b% s_donor% fe_core_mass,            &
            b% s_donor% power_h_burn, b% s_donor% power_he_burn,           &
            b% s_donor% power_c_burn

         end if

         model_number_last = b% s_donor% model_number

         close (59)

      end subroutine rcl_jdot_ml 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer function how_many_extra_binary_history_header_items(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_header_items = 0
      end function how_many_extra_binary_history_header_items


      subroutine data_for_extra_binary_history_header_items( &
           binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id, n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         !names(1) = "Mdot_b"
         !vals(1) = b% ixtra(1)
         !names(1) = "mdot_edd_eta"
         !vals(1) = ixtra(1)
      end subroutine data_for_extra_binary_history_header_items


      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_columns = 2
      end function how_many_extra_binary_history_columns


      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         names(1) = "Mdot_b"
         vals(1) = b% xtra(2)
         names(2) = "mdot_edd_eta"
         vals(2) = b% xtra(1)

         !names(2) = "dipole_moment"
         !names(3) = "r_m"
         !names(4) = "r_{co}"
         !names(5) = "nu_(Hz)"
         !names(6) = "r_{lc}"
         !names(4) = "Period_days_kepler"
         !vals(1) = 0.028*(b% m(b% a_i)/(1.989d33))**(-0.33)
         !vals(2) = 1d8*vals(1)
         !vals(3) = 1*((vals(2)**4)/(2*(6.67d-11)*(b% component_mdot(b% a_i))*(b% m(b% a_i)/(1.989d33))**2))**0.143
         !vals(4) = ((GM)/(4*pi**2*nu))
         !vals(5) = 1
         !vals(6) = c/2*pi*vals(5)
         !vals(4) = 2*pi*sqrt((pow3(b% separation))/(standard_cgrav*(b% m(b% a_i) + b% m(b% d_i))))
         
         
      end subroutine data_for_extra_binary_history_columns

      
!       subroutine cz_jdot_mb(binary_id, ierr)      !consider convective envelop hoc factor
!         use mlt_def
!         integer model_number_last
!         real :: mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
!                 cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
!                 ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
!                 rlo_start_age,rlo_age,radio_start_age,radio_age,vesc_sur_d,eva_wind,&
!                 mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,hoc_accretor,hoc_donor
!         !lp_chen: Chen,H-L et al.,2013 ApJ
!         integer recyelled
!         common /binary/   mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
!                           cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
!                           ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
!                           rlo_start_age,rlo_age,radio_start_age,radio_age,vesc_sur_d,eva_wind,&
!                           mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,&
!                           hoc_accretor,hoc_donor,model_number_last,recyelled
!
!         integer, intent(in) :: binary_id
!         integer, intent(out) :: ierr
!         integer :: i, k, id
!         type (binary_info), pointer :: b
!         type (star_info), pointer :: s
!         real(dp) :: rsun4,two_pi_div_p3
!         ierr = 0
!         call binary_ptr(binary_id, b, ierr)
!         if (ierr /= 0) then
!            write(*,*) 'failed in binary_ptr'
!            return
!         end if
!         b% jdot_mb = 0
!         rsun4 = rsun*rsun*rsun*rsun
!         call check_radiative_core(b)
!         call check_convective_envelop(b)     !hoc_factor
!
!         two_pi_div_p3 = (2.0*pi/b% period)*(2.0*pi/b% period)*(2.0*pi/b% period)
!         ! use the formula from rappaport, verbunt, and joss.  apj, 275, 713-731. 1983.
!         if (b% have_radiative_core(b% d_i) .or. b% keep_mb_on) &
!            b% jdot_mb = -3.8d-30*b% m(b% d_i)*rsun4* &         
!                           pow_cr(min(b% r(b% d_i),b% rl(b% d_i))/rsun,b% magnetic_braking_gamma)* &
!                           two_pi_div_p3
!            b% jdot_mb = b% jdot_mb * hoc_donor
!         if (b% evolve_both_stars .and. b% include_accretor_mb .and. &
!             (b% have_radiative_core(b% a_i) .or. b% keep_mb_on)) then
!             b% jdot_mb = b% jdot_mb - &
!                           3.8d-30*b% m(b% a_i)*rsun4* &
!                           pow_cr(min(b% r(b% a_i),b% rl(b% a_i))/rsun,b% magnetic_braking_gamma)* &
!                           two_pi_div_p3
!             b% jdot_mb = b% jdot_mb * hoc_accretor
!         end if
!
!
!      end subroutine cz_jdot_mb
!
!      subroutine cz_jdot_mb(binary_id, ierr)      !consider convective envelop hoc factor
!         use mlt_def
!         integer model_number_last
!         real :: mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
!                 cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
!                 ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
!                 rlo_start_age,rlo_age,radio_start_age,radio_age,vesc_sur_d,eva_wind,&
!                 mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,hoc_accretor,hoc_donor
!         !lp_chen: Chen,H-L et al.,2013 ApJ
!         integer recyelled
!         common /binary/   mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
!                           cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
!                           ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
!                           mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,&
!                           hoc_accretor,hoc_donor,model_number_last,recyelled
!
!         integer, intent(in) :: binary_id
!         integer, intent(out) :: ierr
!         integer :: i, k, id
!         type (binary_info), pointer :: b
!         type (star_info), pointer :: s
!         real(dp) :: rsun4,two_pi_div_p3
!         ierr = 0
!         call binary_ptr(binary_id, b, ierr)
!         if (ierr /= 0) then
!            write(*,*) 'failed in binary_ptr'
!            return
!         end if
!         b% jdot_mb = 0
!         rsun4 = rsun*rsun*rsun*rsun
!         call check_radiative_core(b)
!         call check_convective_envelop(b)     !hoc_factor
!
!         two_pi_div_p3 = (2.0*pi/b% period)*(2.0*pi/b% period)*(2.0*pi/b% period)
!         ! use the formula from rappaport, verbunt, and joss.  apj, 275, 713-731. 1983.
!         if (b% have_radiative_core(b% d_i) .or. b% keep_mb_on) &
!            b% jdot_mb = -3.8d-30*b% m(b% d_i)*rsun4* &         
!                           pow_cr(min(b% r(b% d_i),b% rl(b% d_i))/rsun,b% magnetic_braking_gamma)* &
!                           two_pi_div_p3
!            b% jdot_mb = b% jdot_mb * hoc_donor
!         if (b% evolve_both_stars .and. b% include_accretor_mb .and. &
!             (b% have_radiative_core(b% a_i) .or. b% keep_mb_on)) then
!             b% jdot_mb = b% jdot_mb - &
!                           3.8d-30*b% m(b% a_i)*rsun4* &
!                           pow_cr(min(b% r(b% a_i),b% rl(b% a_i))/rsun,b% magnetic_braking_gamma)* &
!                           two_pi_div_p3
!             b% jdot_mb = b% jdot_mb * hoc_accretor
!         end if
!
!
!      end subroutine cz_jdot_mb
        
      
      end module run_binary_extras
