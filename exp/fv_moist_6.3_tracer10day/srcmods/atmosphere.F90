!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module atmosphere_mod

!-----------------------------------------------------------------------
!
!    interface for FV dynamical core with Held-Suarez forcing
!
!-----------------------------------------------------------------------


use constants_mod, only: grav, kappa, cp_air, pi, rdgas, rvgas, SECONDS_PER_DAY, &
                         cp_vapor, cp_water, latent_heat, rho_cp, rho0, r_vir, cp_vir, &
                         ps0, deltat_rad

use fms_mod,       only: file_exist, open_namelist_file,   &
                         error_mesg, FATAL,                &
                         check_nml_error, stdlog, stdout,  &
                         write_version_number,             &
                         close_file, set_domain, nullify_domain, mpp_pe, mpp_root_pe
use time_manager_mod, only: time_type, get_time, set_time, operator(+)
use mpp_domains_mod,  only: domain2d

!------------------
! FV specific codes:
!------------------
use fv_arrays_mod, only: fv_atmos_type
use fv_control_mod,only: fv_init, domain, fv_end, adiabatic, p_ref
use fv_phys_mod,   only: fv_phys, fv_nudge
use fv_diagnostics_mod, only: fv_diag_init, fv_diag, fv_time
use fv_timing_mod,   only: timing_on, timing_off
use fv_restart_mod, only: fv_restart
use fv_dynamics_mod, only: fv_dynamics
use fv_grid_tools_mod, only: grid_type, globalsum
!use lin_cld_microphys_mod, only: lin_cld_microphys_init, lin_cld_microphys_end

!FD, LBL
use diag_manager_mod,   only: register_diag_field, send_data
use radiance_mod,       only: radiance_init, radiance_update
!use qe_moist_convection_mod, only: dry_convection 
!FD, PBL
use vert_turb_driver_mod, only: vert_turb_driver_init, vert_turb_driver, vert_turb_driver_end
use vert_diff_mod,      only: vert_diff_init, gcm_vert_diff_down, gcm_vert_diff_up, vert_diff_end, surf_diff_type
use surface_flux_mod,   only: surface_flux
use mixed_layer_mod,    only: mixed_layer_init, mixed_layer, mixed_layer_end
use fv_update_phys_mod, only: update_dwinds_phys 
use fv_mp_mod,          only: gid, masterproc
!FD, moist
use lscale_cond_mod, only: lscale_cond_init, lscale_cond
use  qe_moist_convection_mod, only: qe_moist_convection_init, &
                                    qe_moist_convection,   &
                                    qe_moist_convection_end
use lastsaturation_tracer_mod, only: lastsaturation_tracer
!-----------------------------------------------------------------------

implicit none
private

public   atmosphere_init, atmosphere,  atmosphere_end, atmosphere_domain

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: atmosphere.F90,v 18.0 2010/03/02 23:26:58 fms Exp $'
character(len=128) :: tag = '$Name: riga $'
character(len=10), parameter :: mod_name='atmosphere'

!-----------------------------------------------------------------------
!---- private data ----

type        (time_type) :: Time_step_atmos
real                    :: dt_atmos
integer :: sec
integer days, seconds

logical :: cold_start      = .false.       ! read in initial condition

type(fv_atmos_type), allocatable :: Atm(:)
!-----------------------------------------------------------------------
!FD
logical :: do_bettsmiller = .false.
!FD, PBL
logical :: used, doing_edt, doing_entrain
logical :: module_is_initialized =.false.
logical :: turb = .false.
logical :: do_virtual = .false. ! whether virtual temp used in gcm_vert_diff
real :: roughness_heat = 0.05
real :: roughness_moist = 0.05
real :: roughness_mom = 0.05

namelist/atmosphere_nml/ do_bettsmiller

!===================================================================
real, allocatable, dimension(:,:)   ::                                        &
     z_surf,               &   ! surface height
     t_surf,               &   ! surface temperature
     q_surf,               &   ! surface moisture
     u_surf,               &   ! surface U wind
     v_surf,               &   ! surface V wind
     rough_mom,            &   ! momentum roughness length for surface_flux
     rough_heat,           &   ! heat roughness length for surface_flux
     rough_moist,          &   ! moisture roughness length for surface_flux
     gust,                 &   ! gustiness constant
     z_pbl,                &   ! gustiness constant
     flux_t,               &   ! surface sensible heat flux
     flux_q,               &   ! surface moisture flux
     flux_r,               &   ! surface radiation flux
     flux_u,               &   ! surface flux of zonal mom.
     flux_v,               &   ! surface flux of meridional mom.
     drag_m,               &   ! momentum drag coefficient
     drag_t,               &   ! heat drag coefficient
     drag_q,               &   ! moisture drag coefficient
     w_atm,                &   ! wind speed
     ustar,                &   ! friction velocity
     bstar,                &   ! buoyancy scale
     qstar,                &   ! moisture scale
     dhdt_surf,            &   ! d(sensible heat flux)/d(surface temp)
     dedt_surf,            &   ! d(latent heat flux)/d(surface temp)???
     dedq_surf,            &   ! d(latent heat flux)/d(surface moisture)???
     drdt_surf,            &   ! d(upward longwave)/d(surface temp)
     dhdt_atm,             &   ! d(sensible heat flux)/d(atmos.temp)
     dedq_atm,             &   ! d(latent heat flux)/d(atmospheric mixing rat.)
     dtaudv_atm,           &   ! d(stress component)/d(atmos wind)
     dtaudu_atm,           &   ! d(stress component)/d(atmos wind)
     fracland,             &   ! fraction of land in gridbox
     rough                     ! roughness for vert_turb_driver

real, allocatable, dimension(:,:,:) :: tg_tmp, qg_tmp
real, allocatable, dimension(:,:,:) :: p_full, p_half, z_full, z_half
real, allocatable, dimension(:,:,:) :: dt_ug, dt_vg, dt_tg
real, allocatable, dimension(:,:,:) :: dt_ua, dt_va
real, allocatable, dimension(:,:,:,:) :: dt_tracers
real, allocatable, dimension(:,:,:) ::                                        &
     diff_m,               &   ! momentum diffusion coeff.
     diff_t,               &   ! temperature diffusion coeff.
     tdtlw,                &   ! place holder. appears in calling arguments of vert_turb_driver but not used unless do_edt=.true. -- pjp
     diss_heat,            &   ! heat dissipated by vertical diffusion
     non_diff_dt_ug,       &   ! zonal wind tendency except from vertical diffusion
     non_diff_dt_vg,       &   ! merid. wind tendency except from vertical diffusion
     non_diff_dt_tg,       &   ! temperature tendency except from vertical diffusion
     non_diff_dt_qg,       &   ! moisture tendency except from vertical diffusion
     conv_dt_tg,           &   ! temperature tendency from convection
     conv_dt_qg,           &   ! moisture tendency from convection
     cond_dt_tg,           &   ! temperature tendency from condensation
     cond_dt_qg                ! moisture tendency from condensation

logical, allocatable, dimension(:,:) ::                                       &
     avail,                &   ! generate surf. flux (all true)
     land,                 &   ! land points (all false)
     coldT,                &   ! should precipitation be snow at this point
     convect                   ! place holder. appears in calling arguments of vert_turb_driver but not used unless do_entrain=.true. -- pjp

real, allocatable, dimension(:,:) ::                                          &
     klzbs,                &   ! stored level of zero buoyancy values
     cape,                 &   ! convectively available potential energy
     cin,                  &   ! convective inhibition (this and the above are before the adjustment)
     invtau_q_relaxation,  &   ! temperature relaxation time scale
     invtau_t_relaxation,  &   ! humidity relaxation time scale
     rain,                 &   !
     snow

real, allocatable, dimension(:,:,:) :: &
     t_ref,          &   ! relaxation temperature for bettsmiller scheme
     q_ref               ! relaxation moisture for bettsmiller scheme

integer, allocatable, dimension(:,:) :: & convflag ! indicates which qe convection subroutines are used

integer ::           &
     id_diff_dt_ug,  &   ! zonal wind tendency from vertical diffusion
     id_diff_dt_vg,  &   ! merid. wind tendency from vertical diffusion
     id_diff_dt_tg,  &   ! temperature tendency from vertical diffusion
     id_diff_dt_qg,  &   ! moisture tendency from vertical diffusion
     id_conv_rain,   &   ! rain from convection
     id_cond_rain,   &   ! rain from condensation
     id_conv_dt_tg,  &   ! temperature tendency from convection
     id_conv_dt_qg,  &   ! temperature tendency from convection
     id_cond_dt_tg,  &   ! temperature tendency from convection
     id_cond_dt_qg       ! temperature tendency from convection

integer :: id_phalf, id_pfull, id_zhalf, id_zfull
integer :: id_drag, id_ustar, id_bstar
integer :: id_fluxt, id_ediff
real,    allocatable, dimension(:,:) :: ediff

type(surf_diff_type) :: Tri_surf ! used by gcm_vert_diff
integer i,j, isc, iec, jsc, jec, num_levels

!=================================================================================================================================

contains

!#######################################################################

  subroutine atmosphere_init ( Time_init, Time, Time_step )

    type (time_type), intent(in) :: Time_step
    type (time_type), intent(in) :: Time_init
    type (time_type), intent(in) :: Time

    ! local:
    integer :: axes(5) !4)
    integer :: ss, ds
    integer :: ntiles=1
!    integer i,j, isc, iec, jsc, jec, num_levels
    integer :: unit, ierr, io
    real pp(2)


  !----- write version and namelist to log file -----

    call write_version_number ( version, tag )

  !---- compute physics/atmos time step in seconds ----

    Time_step_atmos = Time_step
    call get_time (Time_step_atmos, sec)
    dt_atmos = real(sec)

!FD, read namelist
unit = open_namelist_file ()
ierr=1
do while (ierr /= 0)
  read  (unit, nml=atmosphere_nml, iostat=io, end=10)
  ierr = check_nml_error (io, 'atmosphere_nml')
enddo
10 call close_file (unit)
if ( mpp_pe() == mpp_root_pe() )   write (stdlog(), nml=atmosphere_nml)


  !----- initialize FV dynamical core -----
!   cold_start = (.not.file_exist('INPUT/fv_rst.res.nc').and. .not.file_exist('INPUT/'//fms_tracers_file))
    cold_start = (.not.file_exist('INPUT/fv_core.res.nc'))

    allocate(Atm(ntiles))
    call fv_init(Atm(:),dt_atmos)  ! allocates Atm components

    Atm(1)%moist_phys = .false.

    ! Init model data
         call timing_on('fv_restart')
    call fv_restart(domain, Atm, dt_atmos, seconds, days, cold_start, grid_type)
         call timing_off('fv_restart')

#ifdef PERTURB_IC
    isc = Atm(1)%isc
    iec = Atm(1)%iec
    jsc = Atm(1)%jsc
    jec = Atm(1)%jec
#ifdef SW_DYNAMICS
! TEST CASE-7
!     if ( .not. cold_start ) then
         do j=jsc,jec
            do i=isc,iec
               pp(1) = Atm(1)%agrid(i,j,1)
               pp(2) = Atm(1)%agrid(i,j,2)
               Atm(1)%delp(i,j,1) = Atm(1)%delp(i,j,1) + 120.*grav*cos(pp(2)) *  &
               exp( -(3.*(pp(1)-pi))**2 ) * exp( -(15.*(pp(2)-pi/4.))**2 )
            enddo
         enddo
!     endif
#endif
#endif

! Need to implement time in fv_restart files
!   if ( .not. cold_start ) then 
! !
! ! Check consistency in Time
! !
!     fv_time = set_time (seconds, days)
!     call get_time (Time, ss,  ds)
!
!     if(seconds /= ss .or. days /= ds) then
!        unit = stdout()
!       write(unit,*) 'FMS:', ds, ss
!       write(unit,*) 'FV:', days, seconds
!       call error_mesg('FV_init:','Time inconsistent between fv_rst and INPUT/atmos_model.res',FATAL)
!     endif
!   else
      fv_time = time
!   endif



    call fv_diag_init(Atm, axes, Time, Atm(1)%npx, Atm(1)%npy, Atm(1)%npz, p_ref)
    !call lin_cld_microphys_init(axes, Time)

!   if( nlev > 1 ) call hs_forcing_init ( axes, Time )

!-----------------------------------------------------------------------
!FD, ts_init
    isc = Atm(1)%isc
    iec = Atm(1)%iec
    jsc = Atm(1)%jsc
    jec = Atm(1)%jec
    num_levels = Atm(1)%npz
    if (cold_start) then
       do j=jsc,jec
          do i=isc,iec
             Atm(1)%ts(i,j) = Atm(1)%pt(i,j,Atm(1)%npz) + 1.
          enddo
       enddo
       Atm(1)%q(:,:,:,1) = 1e-6  !has to be zero for pure CO2 rad
    else
!u_srf here is used for restart
        Atm(1)%ts(:,:) = Atm(1)%u_srf(:,:)
    endif
    Atm(1)%q(:,:,:, 5:60) = 0.

!FD, PBL
allocate(dt_ua(Atm(1)%isd:Atm(1)%ied, &
        Atm(1)%jsd:Atm(1)%jed, num_levels))
allocate(dt_va(Atm(1)%isd:Atm(1)%ied, &
        Atm(1)%jsd:Atm(1)%jed, num_levels))

allocate(tg_tmp      (isc:iec, jsc:jec, num_levels))
allocate(qg_tmp      (isc:iec, jsc:jec, num_levels))
allocate(dt_ug       (isc:iec, jsc:jec, num_levels))
allocate(dt_vg       (isc:iec, jsc:jec, num_levels))
allocate(dt_tg       (isc:iec, jsc:jec, num_levels))
allocate(dt_tracers  (isc:iec, jsc:jec, num_levels, Atm(1)%ncnst))
allocate(non_diff_dt_tg (isc:iec, jsc:jec, num_levels))
allocate(non_diff_dt_ug (isc:iec, jsc:jec, num_levels))

allocate(p_full      (isc:iec, jsc:jec, num_levels))
allocate(z_full      (isc:iec, jsc:jec, num_levels))
allocate(p_half      (isc:iec, jsc:jec, num_levels+1))
allocate(z_half      (isc:iec, jsc:jec, num_levels+1)); z_half = 0.
allocate(z_surf      (isc:iec, jsc:jec)); z_surf = 0.0
!allocate(t_surf      (isc:iec, jsc:jec))
allocate(q_surf      (isc:iec, jsc:jec)); q_surf = 0.0
allocate(u_surf      (isc:iec, jsc:jec)); u_surf = 0.0
allocate(v_surf      (isc:iec, jsc:jec)); v_surf = 0.0
allocate(rough_mom   (isc:iec, jsc:jec)); rough_mom = roughness_mom
allocate(rough_heat  (isc:iec, jsc:jec)); rough_heat = roughness_heat
allocate(rough_moist (isc:iec, jsc:jec)); rough_moist = roughness_moist
allocate(gust        (isc:iec, jsc:jec)); gust = 1.0
allocate(z_pbl       (isc:iec, jsc:jec))
allocate(flux_t      (isc:iec, jsc:jec))
allocate(flux_q      (isc:iec, jsc:jec))
allocate(flux_r      (isc:iec, jsc:jec))
allocate(flux_u      (isc:iec, jsc:jec))
allocate(flux_v      (isc:iec, jsc:jec))
allocate(drag_m      (isc:iec, jsc:jec))
allocate(drag_t      (isc:iec, jsc:jec))
allocate(drag_q      (isc:iec, jsc:jec))
allocate(w_atm       (isc:iec, jsc:jec))
allocate(ustar       (isc:iec, jsc:jec))
allocate(bstar       (isc:iec, jsc:jec))
allocate(qstar       (isc:iec, jsc:jec))
allocate(dhdt_surf   (isc:iec, jsc:jec))
allocate(dedt_surf   (isc:iec, jsc:jec))
allocate(dedq_surf   (isc:iec, jsc:jec))
allocate(drdt_surf   (isc:iec, jsc:jec))
allocate(dhdt_atm    (isc:iec, jsc:jec))
allocate(dedq_atm    (isc:iec, jsc:jec))
allocate(dtaudv_atm  (isc:iec, jsc:jec))
allocate(dtaudu_atm  (isc:iec, jsc:jec))
allocate(land        (isc:iec, jsc:jec)); land = .false.
allocate(avail       (isc:iec, jsc:jec)); avail = .true.
allocate(fracland    (isc:iec, jsc:jec)); fracland = 0.0
allocate(rough       (isc:iec, jsc:jec))
allocate(diff_t      (isc:iec, jsc:jec, num_levels))
allocate(diff_m      (isc:iec, jsc:jec, num_levels))
allocate(diss_heat   (isc:iec, jsc:jec, num_levels))

allocate(tdtlw       (isc:iec, jsc:jec, num_levels)); tdtlw = 0.0
allocate(convect     (isc:iec, jsc:jec)); convect = .false.
allocate(ediff       (isc:iec, jsc:jec)); ediff = 0.       

!moist
allocate(conv_dt_tg  (isc:iec, jsc:jec, num_levels))
allocate(conv_dt_qg  (isc:iec, jsc:jec, num_levels))
allocate(cond_dt_tg  (isc:iec, jsc:jec, num_levels))
allocate(cond_dt_qg  (isc:iec, jsc:jec, num_levels))
allocate(coldT        (isc:iec, jsc:jec)); coldT = .false.
allocate(klzbs        (isc:iec, jsc:jec))
allocate(cape         (isc:iec, jsc:jec))
allocate(cin          (isc:iec, jsc:jec))
allocate(invtau_q_relaxation  (isc:iec, jsc:jec))
allocate(invtau_t_relaxation  (isc:iec, jsc:jec))
allocate(rain         (isc:iec, jsc:jec)); rain = 0.0
allocate(snow         (isc:iec, jsc:jec)); snow = 0.0
allocate(convflag     (isc:iec, jsc:jec))
allocate(t_ref (isc:iec, jsc:jec, num_levels)); t_ref = 0.0
allocate(q_ref (isc:iec, jsc:jec, num_levels)); q_ref = 0.0

!LBL
do i=1,num_levels
!   p_half(:,:,i) = Atm(1)%pe(isc:iec, i, jsc:jec)
   p_full(:,:,i) = Atm(1)%delp(isc:iec,jsc:jec,i) / &
        (Atm(1)%peln(:,i+1,:) - Atm(1)%peln(:,i,:))
enddo

!Atm(1)%q(:,:,:,1) = 1e-6  !has to be zero for pure CO2 rad
!!for 1D comparison
!Atm(1)%ts(:,:) = 300.
!Atm(1)%pt(isc:iec, jsc:jec, :) = 300* &
!        (p_full(isc:iec, jsc:jec, :) / ps0)**kappa

call radiance_init(isc, iec, jsc, jec, &
                num_levels, axes, Time, &
                Atm(1)%agrid(isc:iec,jsc:jec,:), &
                Atm(1)%q(isc:iec,jsc:jec,:,1), &
                Atm(1)%pt(isc:iec,jsc:jec,:), & 
                Atm(1)%ps(isc:iec,jsc:jec), &
                p_full(:,:,:))
     !from 1 to npz

call mixed_layer_init(isc, iec, jsc, jec, num_levels, t_surf, axes(1:4), Time) 
call vert_diff_init (Tri_surf, iec-isc+1, jec-jsc+1, num_levels, .true., do_virtual)
call vert_turb_driver_init (Atm(1)%agrid(isc:iec,jsc:jec,1), Atm(1)%agrid(isc:iec,jsc:jec,2), & 
        iec-isc+1,jec-jsc+1, &
        num_levels,axes(1:4),Time, doing_edt, doing_entrain)
!isc:iec not boundary

id_phalf = register_diag_field(mod_name, 'ph',        &
     axes((/1,2,4/)), Time, 'half pressure','Pa')
id_pfull = register_diag_field(mod_name, 'pf',        &
     axes((/1,2,3/)), Time, 'full pressure','Pa')
id_zhalf = register_diag_field(mod_name, 'zhalf',        &
     axes((/1,2,4/)), Time, 'half height','m')
id_zfull = register_diag_field(mod_name, 'zfull',        &
     axes((/1,2,3/)), Time, 'full height','m')

id_fluxt = register_diag_field(mod_name, 'fluxt',        &
     axes((/1,2/)), Time, 'sensible heat','W/m**2')
id_ediff = register_diag_field(mod_name, 'ediff',        &
     axes((/1,2/)), Time, 'ediff','W/m**2')

id_drag = register_diag_field(mod_name, 'dragm',        &
     axes((/1,2/)), Time, 'dragm','W/m**2')
id_ustar = register_diag_field(mod_name, 'ustar',        &
     axes((/1,2/)), Time, 'ustar','W/m**2')
id_bstar = register_diag_field(mod_name, 'bstar',        &
     axes((/1,2/)), Time, 'bstar','W/m**2')

call lscale_cond_init()
id_cond_dt_qg = register_diag_field(mod_name, 'dt_qg_condensation',        &
     axes(1:3), Time, 'Moisture tendency from condensation','kg/kg/s')
id_cond_dt_tg = register_diag_field(mod_name, 'dt_tg_condensation',        &
     axes(1:3), Time, 'Temperature tendency from condensation','K/s')
id_cond_rain = register_diag_field(mod_name, 'condensation_rain',          &
     axes(1:2), Time, 'Rain from condensation','kg/m/m/s')

if (do_bettsmiller) then
   call qe_moist_convection_init()
   id_conv_dt_qg = register_diag_field(mod_name, 'dt_qg_convection',          &
        axes(1:3), Time, 'Moisture tendency from convection','kg/kg/s')
   id_conv_dt_tg = register_diag_field(mod_name, 'dt_tg_convection',          &
        axes(1:3), Time, 'Temperature tendency from convection','K/s')
   id_conv_rain = register_diag_field(mod_name, 'convection_rain',            &
        axes(1:2), Time, 'Rain from convection','kg/m/m/s')
endif


  end subroutine atmosphere_init


!#######################################################################

  subroutine atmosphere (Time)
    type(time_type), intent(in) :: Time

    real:: zvir
    real:: time_total
    real:: tau_winds, tau_press, tau_temp

    integer i,j!, isc, iec, jsc, jec
    real :: hr_top

#ifdef NUDGE_IC
    tau_winds =  1. * 3600.
    tau_press = -1.
    tau_temp  = -1.
#else
    tau_winds = -1.
    tau_press = -1.
    tau_temp  = -1.
#endif

    fv_time = Time + Time_step_atmos
    call get_time (fv_time, seconds,  days)

    time_total = days*SECONDS_PER_DAY + seconds

    if ( tau_winds>0. .or. tau_press>0. .or. tau_temp>0. )     &
    call  fv_nudge(Atm(1)%npz, Atm(1)%isc, Atm(1)%iec, Atm(1)%jsc, Atm(1)%jec, Atm(1)%ng, &
                   Atm(1)%u, Atm(1)%v, Atm(1)%delp, Atm(1)%pt, dt_atmos,    &
                   tau_winds, tau_press, tau_temp)

  !---- call fv dynamics -----
  !  if ( adiabatic .or. Atm(1)%do_Held_Suarez ) then
  !       zvir = 0.         ! no virtual effect
  !  else
         zvir = rvgas/rdgas - 1.
  !  endif

    call set_domain(Atm(1)%domain)  ! needed for diagnostic output done in fv_dynamics
    call timing_on('fv_dynamics')
    call fv_dynamics(Atm(1)%npx, Atm(1)%npy, Atm(1)%npz, Atm(1)%ncnst-Atm(1)%pnats,         &
                     Atm(1)%ng, dt_atmos, Atm(1)%consv_te,                                  &
                     Atm(1)%fill, Atm(1)%reproduce_sum, kappa, cp_air, zvir,                &
                     Atm(1)%ks, Atm(1)%ncnst, Atm(1)%n_split, Atm(1)%q_split,               &
                     Atm(1)%u, Atm(1)%v, Atm(1)%um, Atm(1)%vm,                              &
                     Atm(1)%w, Atm(1)%delz, Atm(1)%hydrostatic,                             &
                     Atm(1)%pt, Atm(1)%delp, Atm(1)%q, Atm(1)%ps,                           &
                     Atm(1)%pe, Atm(1)%pk, Atm(1)%peln, Atm(1)%pkz, Atm(1)%phis,            &
                     Atm(1)%omga, Atm(1)%ua, Atm(1)%va, Atm(1)%uc, Atm(1)%vc,               &
                     Atm(1)%ak, Atm(1)%bk, Atm(1)%mfx, Atm(1)%mfy,                          &
                     Atm(1)%cx, Atm(1)%cy, Atm(1)%ze0, Atm(1)%hybrid_z, time_total)
    call timing_off('fv_dynamics')

!DF LBL===================================================================
isc = Atm(1)%isc
iec = Atm(1)%iec
jsc = Atm(1)%jsc
jec = Atm(1)%jec
num_levels = Atm(1)%npz
!p_atm = Atm(1)%delp(isc:iec,jsc:jec,num_levels) / &
!        (Atm(1)%peln(:,num_levels+1,:) - &
!        Atm(1)%peln(:,num_levels,:) )
!z_atm = (p_atm / Atm(1)%ps(isc:iec,jsc:jec)) * &
!        rdgas * Atm(1)%pt(isc:iec,jsc:jec,num_levels) /grav
do i=1,num_levels
   p_half(:,:,i) = Atm(1)%pe(isc:iec, i, jsc:jec)
   p_full(:,:,i) = Atm(1)%delp(isc:iec,jsc:jec,i) / &
        (Atm(1)%peln(:,i+1,:) - Atm(1)%peln(:,i,:))
enddo
p_half(:,:,num_levels+1) = Atm(1)%ps(isc:iec, jsc:jec)
z_half(:,:,1:num_levels) = Atm(1)%delp(isc:iec,jsc:jec,:) / &
        p_full(:,:,:) * rdgas /grav * (1+zvir*Atm(1)%q(isc:iec,jsc:jec,:,1)) *&
        Atm(1)%pt(isc:iec,jsc:jec,:)
z_full(:,:,:) = (p_half(:,:,2:num_levels+1) - p_full(:,:,:)) / &
        p_full(:,:,:) * rdgas /grav * (1+zvir*Atm(1)%q(isc:iec,jsc:jec,:,1)) *&
        Atm(1)%pt(isc:iec,jsc:jec,:)

do i=num_levels, 1, -1
   z_half(:,:,i) = z_half(:,:,i) + z_half(:,:,i+1)
   z_full(:,:,i) = z_full(:,:,i) + z_half(:,:,i+1)
enddo

!!DF DRY CONVEC=================================================================== 
!do i=isc, iec
!   do j=jsc, jec
!      call dry_convection(Atm(1)%pt(i,j,:), p_full(i,j,:), p_half(i,j,:), &
!      tg_tmp(i,j,:))
!   enddo
!enddo
!Atm(1)%pt(isc:iec, jsc:jec, :) = tg_tmp(isc:iec, jsc:jec, :)*1.


!if (mod(seconds, 14400) .eq. int(dt_atmos)) then !time_step_size
!       call radiance_update(Atm(1)%isc, Atm(1)%iec, &
!                Atm(1)%jsc, Atm(1)%jec, &
!                Atm(1)%ng,  Atm(1)%ncnst, Atm(1)%npz, &
!                fv_time, Atm(1)%agrid, &
!                Atm(1)%q, Atm(1)%pt, & 
!                Atm(1)%ts,Atm(1)%delp, Atm(1)%peln, &
!                non_diff_dt_tg       , Atm(1)%dTs_rad)
!                Atm(1)%dTrad, Atm(1)%dTs_rad)
if (mod(seconds, deltat_rad) .eq. int(dt_atmos)) then !time_step_size
       call radiance_update(isc, iec, jsc, jec, &
                num_levels, fv_time, &
                Atm(1)%agrid(isc:iec,jsc:jec,:), &
                Atm(1)%q(isc:iec,jsc:jec,:,1), &
                Atm(1)%pt(isc:iec,jsc:jec,:), & 
                Atm(1)%ts(isc:iec,jsc:jec), &
                p_half(:,:,:), p_full(:,:,:), &
                non_diff_dt_tg(:,:,:), &
                Atm(1)%dTs_rad(isc:iec,jsc:jec))

!!stablize the top layer, without breaking the energy conservation
!hr_top = globalsum(non_diff_dt_tg(isc:iec,jsc:jec,1), &
!        Atm(1)%npx, Atm(1)%npy, isc, iec, jsc, jec)
!non_diff_dt_tg(isc:iec, jsc:jec, 1) = hr_top
!if(gid==masterproc) write(6,*) 'TOA avergaed heating rate',hr_top 

endif

!update ta,ts
!    do i=Atm(1)%isc,Atm(1)%iec
!    do j=Atm(1)%jsc,Atm(1)%jec
!       Atm(1)%pt(i,j,:) = Atm(1)%pt(i,j,:) + &
!        Atm(1)%dTrad(i,j,:) * dt_atmos
!    enddo
!    enddo
!    Atm(1)%ts(:,:) = Atm(1)%ts(:,:) &
!        + Atm(1)%dTs_rad(:,:) * dt_atmos / 5.

!DF PBL===================================================================
dt_tg = non_diff_dt_tg *1. !radiative heating rate now
dt_ug = 0.; dt_vg = 0.; dt_tracers = 0.
dt_ua = 0.; dt_va = 0.

call surface_flux(                                                          &
        Atm(1)%pt(isc:iec,jsc:jec,num_levels),                              &
       Atm(1)%q(isc:iec,jsc:jec,num_levels,1),                              &
        Atm(1)%ua(isc:iec,jsc:jec,num_levels),                              &
        Atm(1)%va(isc:iec,jsc:jec,num_levels),                              &
                       p_full(:,:,num_levels),                              &
           z_full(:,:,num_levels)-z_surf(:,:),                              &
                   Atm(1)%ps(isc:iec,jsc:jec),                              &
                   Atm(1)%ts(isc:iec,jsc:jec),                              &
                   Atm(1)%ts(isc:iec,jsc:jec),                              &
                                  q_surf(:,:),                              & ! is intent(inout)
                                  u_surf(:,:),                              &
                                  v_surf(:,:),                              &
                               rough_mom(:,:),                              &
                              rough_heat(:,:),                              &
                             rough_moist(:,:),                              &
                               rough_mom(:,:),                              & ! using rough_mom in place of rough_scale -- pjp
                                    gust(:,:),                              &
                                  flux_t(:,:),                              & ! is intent(out)
                                  flux_q(:,:),                              & ! is intent(out)
                                  flux_r(:,:),                              & ! is intent(out)
                                  flux_u(:,:),                              & ! is intent(out)
                                  flux_v(:,:),                              & ! is intent(out)
                                  drag_m(:,:),                              & ! is intent(out)
                                  drag_t(:,:),                              & ! is intent(out)
                                  drag_q(:,:),                              & ! is intent(out)
                                   w_atm(:,:),                              & ! is intent(out)
                                   ustar(:,:),                              & ! is intent(out)
                                   bstar(:,:),                              & ! is intent(out)
                                   qstar(:,:),                              & ! is intent(out)
                               dhdt_surf(:,:),                              & ! is intent(out)
                               dedt_surf(:,:),                              & ! is intent(out)
                               dedq_surf(:,:),                              & ! is intent(out)
                               drdt_surf(:,:),                              & ! is intent(out)
                                dhdt_atm(:,:),                              & ! is intent(out)
                                dedq_atm(:,:),                              & ! is intent(out)
                              dtaudu_atm(:,:),                              & ! is intent(out)
                              dtaudv_atm(:,:),                              & ! is intent(out)
                                     dt_atmos,                              &
                                    land(:,:),                              &
                               .not.land(:,:),                              &
                                   avail(:,:)  )


  call vert_turb_driver(            1,                              1, &
                            Time,                 Time+Time_step_atmos, &
                              dt_atmos, tdtlw(:,:,:),    fracland(:,:), &
                 p_half(:,:,:),          p_full(:,:,:), &
                 z_half(:,:,:),          z_full(:,:,:), &
                            ustar(:,:),                     bstar(:,:), &
                            qstar(:,:),                     rough(:,:), &
          Atm(1)%agrid(isc:iec,jsc:jec,2),                convect(:,:), &
          Atm(1)%ua(isc:iec,jsc:jec,:),   Atm(1)%va(isc:iec,jsc:jec,:), &
          Atm(1)%pt(isc:iec,jsc:jec,:),                                 &
          Atm(1)%q(isc:iec,jsc:jec,:,1), Atm(1)%q(isc:iec,jsc:jec,:,:), &
          Atm(1)%ua(isc:iec,jsc:jec,:),   Atm(1)%va(isc:iec,jsc:jec,:), &
          Atm(1)%pt(isc:iec,jsc:jec,:),                                 &
          Atm(1)%q(isc:iec,jsc:jec,:,1), Atm(1)%q(isc:iec,jsc:jec,:,:), &
                          dt_ug(:,:,:),                   dt_vg(:,:,:), &
                          dt_tg(:,:,:),            dt_tracers(:,:,:,1), &
                   dt_tracers(:,:,:,:),                  diff_t(:,:,:), &
                         diff_m(:,:,:),                      gust(:,:), &
                            z_pbl(:,:) )
!
!! Don't zero these derivatives as the surface flux depends implicitly
!! on the lowest level values
!! However it should be noted that these derivatives do not take into
!! account the change in the Monin-Obukhov coefficients, and so are not
!! very accurate.
!
!!$   dtaudv_atm = 0.0
!!$   dhdt_atm   = 0.0
!!$   dedq_atm   = 0.0

!if (.false.) then

   call gcm_vert_diff_down (1, 1,  dt_atmos, &
             Atm(1)%ua(isc:iec,jsc:jec,:),   Atm(1)%va(isc:iec,jsc:jec,:), &
             Atm(1)%pt(isc:iec,jsc:jec,:),                                 &
             Atm(1)%q(isc:iec,jsc:jec,:,1), Atm(1)%q(isc:iec,jsc:jec,:,:), &
             diff_m(:,:,:),         diff_t(:,:,:),          p_half(:,:,:), &
             p_full(:,:,:),  z_full(:,:,:),     flux_u(:,:),  flux_v(:,:), &
                            dtaudu_atm(:,:),              dtaudv_atm(:,:), &
                            dt_ug(:,:,:),                    dt_vg(:,:,:), &
                            dt_tg(:,:,:),             dt_tracers(:,:,:,1), &
                            dt_tracers(:,:,:,:),         diss_heat(:,:,:), &
                            Tri_surf)

!
! update surface temperature
!
   call mixed_layer(                                                       &
                              Time,                                        &
                           Atm(1)%ts(:,:),                                 & ! t_surf is intent(inout)
                              flux_t(:,:),                                 &
                              flux_q(:,:),                                 &
                            fracland(:,:),                                 &
                                 dt_atmos,                                 &
                      Atm(1)%dTs_rad(:,:),                                 &
                            fracland(:,:),                                 & !dTs_rad = swd + lwn
                            Tri_surf,                                      & ! Tri_surf is intent(inout)
                           dhdt_surf(:,:),                                 &
                           dedt_surf(:,:),                                 &
                           dedq_surf(:,:),                                 &
                           drdt_surf(:,:),                                 &
                            dhdt_atm(:,:),                                 &
                            dedq_atm(:,:))

   call gcm_vert_diff_up (1, 1, dt_atmos, Tri_surf, dt_tg(:,:,:), dt_tracers(:,:,:,1), dt_tracers(:,:,:,:))

!endif

!DF Check Energy===================================================================
!non_diff_dt_ug = ((dt_tg - non_diff_dt_tg) * cp_air + &
!        Atm(1)%ua(isc:iec,jsc:jec,:) * dt_ug + &
!        Atm(1)%va(isc:iec,jsc:jec,:) * dt_vg) *&
!        Atm(1)%delp(isc:iec,jsc:jec,:) /grav
!ediff = 0.
!do i=1,num_levels
!   ediff = ediff + non_diff_dt_ug(:,:,i)
!enddo
!DF Betts-Miller convection ===================================================
if (do_bettsmiller) then
   rain = 0.0; snow = 0.0
   call qe_moist_convection (dt_atmos,    Atm(1)%pt(isc:iec,jsc:jec,:),      &
        Atm(1)%q(isc:iec,jsc:jec,:,1),                     p_full     ,      &
                          p_half     ,                           coldT,      &
                                 rain,                            snow,      &
                           conv_dt_tg,                      conv_dt_qg,      &
                                q_ref,                        convflag,      &
                                klzbs,                            cape,      &
                                  cin,             invtau_q_relaxation,      &
                  invtau_t_relaxation,                           t_ref)
   tg_tmp = conv_dt_tg + Atm(1)%pt(isc:iec,jsc:jec,:) 
   qg_tmp = conv_dt_qg + Atm(1)%q(isc:iec,jsc:jec,:,1) 
!  note the delta's are returned rather than the time derivatives
   conv_dt_tg = conv_dt_tg/dt_atmos
   conv_dt_qg = conv_dt_qg/dt_atmos
   rain       = rain/dt_atmos
   dt_tg      = dt_tg + conv_dt_tg
   dt_tracers(:,:,:,1) = dt_tracers(:,:,:,1) + conv_dt_qg
   if(id_conv_dt_qg > 0) used = send_data(id_conv_dt_qg, conv_dt_qg, Time)
   if(id_conv_dt_tg > 0) used = send_data(id_conv_dt_tg, conv_dt_tg, Time)
   if(id_conv_rain > 0) used = send_data(id_conv_rain, rain, Time)
else 
   tg_tmp = Atm(1)%pt(isc:iec,jsc:jec,:) *1.
   qg_tmp = Atm(1)%q(isc:iec,jsc:jec,:,1) *1.
endif

!DF Lagescale condensation, t, q profiles after convection ===================================================
rain = 0.0
call lscale_cond (         tg_tmp,                   qg_tmp,               &
                    p_full(:,:,:),                   p_half(:,:,:),        &
                            coldT,                            rain,        &
                             snow,                      cond_dt_tg,        &
                       cond_dt_qg )                                
cond_dt_tg = cond_dt_tg/dt_atmos
cond_dt_qg = cond_dt_qg/dt_atmos
rain       = rain/dt_atmos      
dt_tg      = dt_tg + cond_dt_tg
dt_tracers(:,:,:,1) = dt_tracers(:,:,:,1) + cond_dt_qg
if(id_cond_dt_qg > 0) used = send_data(id_cond_dt_qg, cond_dt_qg, Time)
if(id_cond_dt_tg > 0) used = send_data(id_cond_dt_tg, cond_dt_tg, Time)
if(id_cond_rain  > 0) used = send_data(id_cond_rain, rain, Time)

!DF Last saturation tracers===================================================
       call lastsaturation_tracer(isc, iec, jsc, jec, &
                num_levels, &
                Atm(1)%agrid(isc:iec,jsc:jec,:), &
                p_full(:,:,:), &
                cond_dt_qg(isc:iec,jsc:jec,:), &
                Atm(1)%q(isc:iec,jsc:jec,:,5:60) )

!DF UPDATE===================================================================
Atm(1)%ua(isc:iec,jsc:jec,:) = Atm(1)%ua(isc:iec,jsc:jec,:) + &
        dt_ug(:,:,:) * dt_atmos
Atm(1)%va(isc:iec,jsc:jec,:) = Atm(1)%va(isc:iec,jsc:jec,:) + &
        dt_vg(:,:,:) * dt_atmos
Atm(1)%pt(isc:iec,jsc:jec,:) = Atm(1)%pt(isc:iec,jsc:jec,:) + &
        dt_tg(:,:,:) * dt_atmos
Atm(1)%q(isc:iec,jsc:jec,:,1) = Atm(1)%q(isc:iec,jsc:jec,:,1) + &
        dt_tracers(:,:,:,1) * dt_atmos
!not now

                                                    call timing_on(' Update_dwinds')
    dt_ua(isc:iec,jsc:jec,:) = dt_ug
    dt_va(isc:iec,jsc:jec,:) = dt_vg
    call update_dwinds_phys(isc, iec, jsc, jec, &
        Atm(1)%isd, Atm(1)%ied, Atm(1)%jsd, Atm(1)%jed, &
        dt_atmos, dt_ua, dt_va, Atm(1)%u, Atm(1)%v)

                                                    call timing_off(' Update_dwinds')

if (id_fluxt > 0) used = send_data(id_fluxt, flux_t, Time)
if (id_ediff > 0) used = send_data(id_ediff, ediff, Time)

if (id_phalf > 0) used = send_data(id_phalf, p_half, Time)
if (id_pfull > 0) used = send_data(id_pfull, p_full, Time)
if (id_zhalf > 0) used = send_data(id_zhalf, z_half, Time)
if (id_zfull > 0) used = send_data(id_zfull, z_full, Time)

if (id_drag > 0) used = send_data(id_drag, drag_m, Time)
if (id_ustar > 0) used = send_data(id_ustar, ustar, Time)
if (id_bstar > 0) used = send_data(id_bstar, bstar, Time)

!DF UPDATE===================================================================
!    if(Atm(1)%npz /=1 .and. .not. adiabatic)then

!                                                         call timing_on('FV_PHYS')
!       call fv_phys(Atm(1)%npx, Atm(1)%npy, Atm(1)%npz, Atm(1)%isc, Atm(1)%iec,          &
!                    Atm(1)%jsc, Atm(1)%jec, Atm(1)%ng, Atm(1)%ncnst,                     &
!                    Atm(1)%u, Atm(1)%v, Atm(1)%w, Atm(1)%pt, Atm(1)%q, Atm(1)%pe,        &
!                    Atm(1)%delp, Atm(1)%peln, Atm(1)%pkz, dt_atmos,                      &
!                    Atm(1)%ua, Atm(1)%va, Atm(1)%phis, Atm(1)%agrid,                     &
!                    Atm(1)%ak, Atm(1)%bk, Atm(1)%ks, Atm(1)%ps, Atm(1)%pk,               &
!                    Atm(1)%u_srf, Atm(1)%v_srf, Atm(1)%ts, Atm(1)%delz,                  &
!                    Atm(1)%hydrostatic, Atm(1)%phys_hydrostatic, Atm(1)%oro, .true.,     &
!                    .false., p_ref, Atm(1)%fv_sg_adj, (mpp_pe()==mpp_root_pe()),         &
!                    Atm(1)%do_Held_Suarez, fv_time, time_total)
!                                                        call timing_off('FV_PHYS')
!    endif

    call nullify_domain()

  !---- diagnostics for FV dynamics -----

    call timing_on('FV_DIAG')

    call fv_diag(Atm, zvir, fv_time, Atm(1)%print_freq)

    call timing_off('FV_DIAG')


! used to restart Ts
    Atm(1)%u_srf(:,:) = Atm(1)%ts(:,:)

 end subroutine atmosphere


 subroutine atmosphere_end

    call get_time (fv_time, seconds,  days)

    !if ( Atm(1)%nwat==6 )    & 
    !call lin_cld_microphys_end
    
    if (do_bettsmiller) call qe_moist_convection_end

    call fv_end(Atm)
    deallocate(Atm)

  end subroutine atmosphere_end

 subroutine atmosphere_domain ( fv_domain )
 type(domain2d), intent(out) :: fv_domain

!  returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos
        
   fv_domain = domain
        
 end subroutine atmosphere_domain

end module atmosphere_mod


