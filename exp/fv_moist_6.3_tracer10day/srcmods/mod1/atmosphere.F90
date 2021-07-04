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


use time_manager_mod, only: time_type, get_time, set_time, operator(+)

use fms_mod,          only: file_exist, open_namelist_file,   &
                            error_mesg, FATAL,                &
                            check_nml_error, stdlog,          &
                            write_version_number,             &
                            close_file, set_domain

  use hs_forcing_mod,   only: hs_forcing_init
  use constants_mod,    only: omega, cp_air, rdgas, kappa, radius, grav, rvgas, &
                           cp_vapor, cp_water, latent_heat, rho_cp, rho0, r_vir, cp_vir
  use dimensions,       only: nS 

!------------------
! FV specific codes:
!------------------
!  use fv_pack
!use            fv_pack, only: nlon, mlat, nlev, beglat, endlat, beglon, endlon, &
!                              u, v, pt, q, ua, va, delp, phis, ps, pe, peln,    &
!                              pk, pkz, omga, u_srf, v_srf, rlonb, rlatb,        &
!                              cold_start, ncnst, pnats, consv_te, ptop,         &
!                              fv_init, fv_domain, fv_end, change_time, map_dt,  &
!                              adiabatic
use fv_pack, only: nlon, mlat, nlev, beglat, endlat, beglon, endlon, &
     rlonb, rlatb, rlat, rlon,       &
     cold_start, ncnst, pnats, consv_te, ptop,         &
     fv_init, fv_domain, fv_end, change_time, map_dt,  &
     adiabatic, restart_format, &
     get_eta_level, p_var, ak, bk, ng_d, master, &
     nt_phys
use update_fv_phys_mod, only: update_fv_phys
use fv_arrays_mod, only: fv_array_sync

  use fv_diagnostics, only: fv_diag_init, fv_diag, fv_time
  use timingModule,   only: timing_on, timing_off
  use fv_restart_mod, only: fv_restart, write_fv_rst
  use fv_dynamics_mod, only: fv_dynamics
  use fv_phys_mod, only: fv_phys

  use  diag_manager_mod, only: register_diag_field, send_data
  use  radiation_mod, only: radiation_init, radiation_down, radiation_up, &
                                    radiation_end
  use  surface_flux_mod, only: surface_flux
  use  qe_moist_convection_mod, only: moist_convection, compute_k, d622, d608
  use  simple_sat_vapor_pres_mod, only: escomp
  use  lscale_cond_mod, only: lscale_cond

  use radiance_mod, only: radiance_init, radiance_update
  
!-----------------------------------------------------------------------

implicit none
private

public   atmosphere_init, atmosphere,  atmosphere_end

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: atmosphere.F90,v 13.0.2.2 2006/05/19 16:44:39 wfc Exp $'
character(len=128) :: tag = '$Name:  $'
character(len=10), parameter :: mod_name='atmosphere'

!-----------------------------------------------------------------------
!---- private data ----

type        (time_type) :: Time_step_atmos
real                    :: dt_atmos
integer :: sec
integer days, seconds

!-----------------------------------------------------------------------
!df 1401010
integer :: i, j
integer :: ij, nx, tsiz, isiz
integer :: is, ie
!real :: eff_heat_capacity = rho_cp*1.
real :: mld = 1.
integer :: id_t_surf, id_conv_rain, id_cond_rain, id_pme
integer :: id_conv_rain_profile, id_cond_rain_profile
integer :: id_flux_t, id_flux_q, id_flux_r, id_rh
logical :: used
real, allocatable, dimension(:,:)   :: t_surf
real, allocatable, dimension(:,:,:) :: p_half, p_full, rh, Tlev
real, allocatable, dimension(:,:)   :: flux_t, flux_q, flux_r
real, allocatable, dimension(:,:)   :: conv_rain, cond_rain, pme
real, allocatable, dimension(:,:)   :: net_surf_sw_down, surf_lw_down
real, allocatable, dimension(:,:,:) :: conv_rain_profile, cond_rain_profile

real, allocatable, dimension(:,:,:) :: dt_tg_rad 
real, allocatable, dimension(:,:,:) :: OLRnu
real, allocatable, dimension(:,:)   :: solar, OLR 
integer ::  id_OLRnu, id_OLR, id_solar  
!-----------------------------------------------------------------------
contains

!#######################################################################

 subroutine atmosphere_init ( Time_init, Time, Time_step )

#include <fv_arrays.h>

 type (time_type), intent(in) :: Time_step
 type (time_type), intent(in) :: Time_init
 type (time_type), intent(in) :: Time
 logical cold_start

! local:
 integer axes(5) !4)
 integer ss, ds

#include <fv_point.inc>
!----- write version and namelist to log file -----

    call write_version_number ( version, tag )

!---- compute physics/atmos time step in seconds ----

    Time_step_atmos = Time_step
    call get_time (Time_step_atmos, sec)
    dt_atmos = real(sec)

!----- initialize FV dynamical core -----

    call fv_init( sec )
    call fv_restart( days, seconds )
   
    if ( .not. cold_start ) then 
!
! Check consistency in Time
!
    fv_time = set_time (seconds, days)
    call get_time (Time, ss,  ds)

    if( seconds /= ss .or. days /= ds ) then
!       write(6,*) 'FMS:', ds, ss
!       write(6,*) 'FV:', days, seconds
        call error_mesg('FV_init:', &
     'Time inconsistent between fv_rst and INPUT/atmos_model.res', &
                         FATAL)
    endif
    else
    fv_time = time
    endif

    call fv_diag_init( axes, Time )

    if( nlev > 1 ) call hs_forcing_init ( axes, Time )

!-----------------------------------------------------------------------
!df 1401010
allocate(t_surf      (1:nlon, beglat:endlat))
!t_surf(:,:) = tsini 
!----------------------------------------------------------------------
    isiz = nlon/4
    nx = nlon/isiz
    tsiz = nx * ( endlat - beglat + 1 )         ! Total loop length
! Calling phys one latitude at a time using a big/fat OpenMP loop
! For cache performance, this could be changed to finer decomposition if needed
!$omp parallel do private(ij,i,j,k,m,is,ie,p_full,p_half)
if (.not. file_exist('INPUT/atmos_tracers.res.nc')) then
    do ij=1,tsiz
       j  = beglat + (ij-1) / nx
       is = 1 + isiz * mod(ij-1, nx)
       ie = is + isiz - 1
       t_surf(is:ie,j) = pt(is:ie,j,nlev) + 1.
    end do
else
    do ij=1,tsiz
       j  = beglat + (ij-1) / nx
       is = 1 + isiz * mod(ij-1, nx)
       ie = is + isiz - 1
       t_surf(is:ie,j) = q(is:ie,j,nlev,2) !-20.
       if (master) write (6,*) "restart Ts"
    end do 
end if

    q(:,:,nlev,2) = 0. !1e-4

allocate(solar       (1:nlon, beglat:endlat))
allocate(OLR         (1:nlon, beglat:endlat))
allocate(OLRnu       (1:nlon, beglat:endlat, nS))
allocate(dt_tg_rad   (1:nlon, beglat:endlat, nlev))
allocate(Tlev        (1:nlon, beglat:endlat, nlev+1)) !for radiative calculation

allocate(p_half      (1:nlon, beglat:endlat, nlev+1))
allocate(p_full      (1:nlon, beglat:endlat, nlev))
allocate(flux_t      (1:nlon, beglat:endlat))
allocate(flux_q      (1:nlon, beglat:endlat))
allocate(flux_r      (1:nlon, beglat:endlat))
allocate(conv_rain   (1:nlon, beglat:endlat))
allocate(cond_rain   (1:nlon, beglat:endlat))
allocate(pme         (1:nlon, beglat:endlat))
allocate(net_surf_sw_down   (1:nlon, beglat:endlat))
allocate(surf_lw_down       (1:nlon, beglat:endlat))
allocate(conv_rain_profile  (1:nlon, beglat:endlat, nlev))
allocate(cond_rain_profile  (1:nlon, beglat:endlat, nlev))
allocate(rh                 (1:nlon, beglat:endlat, nlev))

p_half = 0.; p_full =0.
!     ----- register diagnostic fields -----
id_t_surf = register_diag_field(mod_name, 't_surf',        &
     axes(1:2), Time, 'Surface temperature','Kelvin')

id_flux_t = register_diag_field(mod_name, 'flux_t',        &
     axes(1:2), Time, 'Sensible heat flux','W/m**2')

id_flux_q = register_diag_field(mod_name, 'flux_q',        &
     axes(1:2), Time, 'Evaporation flux','W/m**2')

id_flux_r = register_diag_field(mod_name, 'flux_r',        &
     axes(1:2), Time, 'Upward IR at surface','W/m**2')

id_conv_rain = register_diag_field(mod_name, 'conv_rain',        &
     axes(1:2), Time, 'Convection rain','kg/s/m**2')

id_cond_rain = register_diag_field(mod_name, 'cond_rain',        &
     axes(1:2), Time, 'Condensation rain','kg/s/m**2')

id_pme = register_diag_field(mod_name, 'pme',        &
     axes(1:2), Time, 'P-E','kg/s/m**2')

id_conv_rain_profile = register_diag_field(mod_name, 'conv_profile',  &
     axes(1:3), Time, 'Convection rain profile','kg/s/m**2')

id_cond_rain_profile = register_diag_field(mod_name, 'cond_profile',  &
     axes(1:3), Time, 'Condensation rain profile','kg/s/m**2')

id_rh = register_diag_field(mod_name, 'rh',        &
     axes(1:3), Time, 'Relative humidity','%')

id_OLRnu = register_diag_field(mod_name, 'OLRnu',        &
     axes((/1,2,5/)), Time, 'OLR spectrum','W/m**2/cm-1') 

id_OLR   = register_diag_field(mod_name, 'OLR',        &
     axes((/1,2/)), Time, 'OLR','W/m**2') 

id_solar   = register_diag_field(mod_name, 'solar',        &
     axes((/1,2/)), Time, 'solar','W/m**2') 

!-----------------------------------------------------------------------
!call radiation_init(1, nlon, beglat, endlat, nlev, axes, Time,rlat(:,:))

!if(master) write(6,*) q(1,beglat,:,1)
do i=1,nlon
   do j=beglat,endlat
      call get_eta_level(nlev, ps(i,j), p_full(i,j,:), p_half(i,j,:))
   end do
end do
call radiance_init(1, nlon, beglat, endlat, nlev, axes, Time, rlat(:,:), &
        q(1:nlon,beglat:endlat,nlev:1:-1,1), p_full(:,:,nlev:1:-1), &
        pt(1:nlon,beglat:endlat,nlev:1:-1), ps(:,:))
dt_tg_rad = 0.
!-----------------------------------------------------------------------
!allocate(OLRs                 (1:nlon, beglat:endlat, wn))
!id_OLRs = register_diag_field(mod_name, 'OLRs',        &
!     axes((/1,2,5/)), Time, 'radiance spectrum','%')

!-----------------------------------------------------------------------
 end subroutine atmosphere_init


!#######################################################################

 subroutine atmosphere (Time)
#include <fv_arrays.h>
 type(time_type), intent(in) :: Time

! local:
 real    zvir
 logical p_map                      ! Perform only partial remapping
#include <fv_point.inc>

!df 141011 used in my physics schemes---------------------------------
real, dimension(1:nlon, beglat:endlat, nlev) :: tg_tmp, qg_tmp, cp, &
                                        u_tmp, v_tmp, rh_tmp, esat, &
                                        dt_ug, dt_vg, dt_tg, dt_qg, &
                                        diff_u, diff_v
real, dimension(1:nlon, beglat:endlat) :: q1, p1, dp1, dmass, lmass, rl, rt, ct, kappa1, t11, &
                                 zh, Ee, sigma, psg_tmp, Ek, dt_psg, &
                                 delta_t_surf_rad

real, dimension(1:nlon, beglat:endlat) ::      &
     gust,                 &   !gustiness constant
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
     dtaudv_atm                ! d(stress component)/d(atmos wind)
logical, dimension(1:nlon, beglat:endlat) ::    land, avail

real, dimension(nlev) :: t_prev, q_prev, pfull_prev, t_after, q_after, &
                               pfull_after, rain_profile, lnpfull_prev, &
                               lnp_full, water, u_after, v_after, &
                               tmean, qmean
                       
real, dimension(nlev+1) :: phalf_prev, phalf_after, lnp_half

real, dimension(1:nlon, nlev) :: tm, qm, um, vm
real, dimension(1:nlon)       :: psm, tsm
real :: psmean, tsmean

integer :: k,n  
integer :: ngroup, nn
!integer :: i, j, k, n
!integer ij, nx, tsiz, isiz
!integer is, ie

real :: Ep, delta_t_surf, Emass
real :: lh, rain1
real :: k3, k4, Ep_2, water_1, water_2
real :: vcoeff, delta_t
real :: ftop
real gmean

land = .false.; avail = .true.; gust = 1.
dt_ug = 0.; dt_vg = 0.; dt_tg = 0.; dt_qg = 0.; dt_psg = 0.
diff_u = 0.; diff_v = 0.
delta_t = dt_atmos

! fix initial q profile
    isiz = nlon/4
    nx = nlon/isiz
    tsiz = nx * ( endlat - beglat + 1 )         ! Total loop length
! Calling phys one latitude at a time using a big/fat OpenMP loop
! For cache performance, this could be changed to finer decomposition if needed
!$omp parallel do private(ij,i,j,k,m,is,ie,p_full,p_half)
    do ij=1,tsiz              
       j  = beglat + (ij-1) / nx
       is = 1 + isiz * mod(ij-1, nx)
       ie = is + isiz - 1
       do i=is,ie
          qg_tmp(i,j,:) = q(i,j,:,1) 
       end do
    end do


!-----------------------------------------------------------------------

  fv_time = Time + Time_step_atmos
  call get_time (fv_time, seconds,  days)

!---- call fv dynamics -----
  zvir = 0.         ! no virtual effect if not full physics
  zvir = d608

  if ( mod(seconds, map_dt) == 0 ) then
       p_map = .false.
  else
       p_map = .true.
  endif
                                call timing_on('fv_dynamics')

#ifndef USE_LIMA
  call fv_dynamics (nlon,    mlat,   nlev,    beglat,   endlat,    &
                    ncnst,   pnats,  p_map,   consv_te,            &
                    u,       v,      delp,    pt,       q,         &
                    ps,      pe,     pk,      pkz,      phis,      &
                    omga,    peln,   ptop,    omega,    sec,       &  
                    zvir,    cp_air, rdgas,   kappa,  radius, ua, va, fv_time )
#else
  call fv_dynamics (nlon,    mlat,   nlev,    beglat,   endlat,    &
                    ncnst,   pnats,  p_map,   consv_te,            &
                    u,       v,      delp,    pt,       q,         &
                    ps,      pe,     pk,      pkz,      phis,      &
                    omga,    peln,   ptop,    omega,    sec,       &  
                    zvir,    cp_air, rdgas,   kappa,  radius, ua, va )
#endif

                                call timing_off('fv_dynamics')

      if( nlev /=1 .and. .not. adiabatic ) then
                                call timing_on('FV_PHYS')
!                call fv_phys ( fv_time, sec )
!df physics 141011
!----------------------------------------------------------------------
    isiz = nlon/4
    nx = nlon/isiz
    tsiz = nx * ( endlat - beglat + 1 )         ! Total loop length
! Calling phys one latitude at a time using a big/fat OpenMP loop
! For cache performance, this could be changed to finer decomposition if needed
!$omp parallel do private(ij,i,j,k,m,is,ie,p_full,p_half)
    do ij=1,tsiz              
       j  = beglat + (ij-1) / nx
       is = 1 + isiz * mod(ij-1, nx)
       ie = is + isiz - 1
!       do k=1,nlev
       do i=is,ie
          call get_eta_level(nlev, ps(i,j), p_full(i,j,:), p_half(i,j,:))
          u_tmp(i,j,:) = ua(i,j,:)
          v_tmp(i,j,:) = va(i,j,:)
          tg_tmp(i,j,:) = pt(i,j,:)
!          qg_tmp(i,j,:) = q(i,j,:,1)
          psg_tmp(i,j) = ps(i,j)
          !interpolate Tlev for radiative calculation
          lnp_half = log(p_half(i,j,:))
          lnp_full = log(p_full(i,j,:))
          call interp(nlev, lnp_full, tg_tmp(i,j,:), &
                nlev+1, lnp_half, Tlev(i,j,:))
          Tlev(i,j,1) = 2.*tg_tmp(i,j,1) - Tlev(i,j,2)
          Tlev(i,j,nlev+1) = 2.*tg_tmp(i,j,nlev) - Tlev(i,j,nlev)
       end do
    end do

       n = nlev
       cp = cp_air * (1 - qg_tmp) + cp_vapor * qg_tmp

!       call radiation_down(1, beglat, Time,                   &
!                       rlat(:,:),                &
!                       rlon(:,:),                &
!                       p_half(:,:,:),                  &
!                       qg_tmp(:,:,:),      &
!                       tg_tmp(:,:,:),             &
!                       net_surf_sw_down(:,:),          &
!                       surf_lw_down(:,:))

!       call radiation_up(1, beglat, Time,                   &
!                     rlat(:,:),                &
!                     p_half(:,:,:),                  &
!                     t_surf(:,:),                    &
!                     tg_tmp(:,:,:),             &
!                     dt_tg(:,:,:))

        if (mod(seconds, 14400) .eq. int(dt_atmos)) then !time_step_size
        call radiance_update(1, beglat, Time, &
                rlat(:,:), rlon(:,:), &
                qg_tmp(:,:,nlev:1:-1), &
                tg_tmp(:,:,nlev:1:-1), &
                Tlev(:,:,nlev+1:1:-1), &
                t_surf(:,:), psg_tmp(:,:), &
                p_full(:,:,nlev:1:-1), &
                p_half(:,:,nlev+1:1:-1), &
                dt_tg_rad(:,:,:), &
                net_surf_sw_down(:,:), &
                surf_lw_down(:,:), &
                solar(:,:), &
                OLR(:,:), &
                OLRnu(:,:,:))  !dt_tg is flux divergence here
        end if

       dt_tg = dt_tg_rad(:,:,nlev:1:-1) *grav / cp(:,:,:) &
        / (p_half(:,:,2:n+1) - p_half(:,:,1:n))
       tg_tmp = tg_tmp + dt_tg * delta_t

       delta_t_surf_rad = (surf_lw_down + net_surf_sw_down  &
                     - 5.67e-8*t_surf **4) &
                    * delta_t / rho_cp / mld !eff_heat_capacity
       t_surf = t_surf + delta_t_surf_rad

    isiz = nlon/4
    nx = nlon/isiz
    tsiz = nx * ( endlat - beglat + 1 )         ! Total loop length
! Calling phys one latitude at a time using a big/fat OpenMP loop
! For cache performance, this could be changed to finer decomposition if needed
!$omp parallel do private(ij,i,j,k,m,is,ie,p_full,p_half)
    do ij=1,tsiz              
       j  = beglat + (ij-1) / nx
       is = 1 + isiz * mod(ij-1, nx)
       ie = is + isiz - 1
       do i=is,ie
          pt(i,j,:) = tg_tmp(i,j,:) 
          q(i,j,:,1)= qg_tmp(i,j,:)
       end do
    end do


!if (master) write(6,*) 'sec=',seconds,'s'
!if (master) write(6,*) t_surf

!if(master) write(6,*) 'surface temperature correction (K)=', ftop
!if(master) write(6,*) 'p_map=', p_map
!if(master) write(6,*) r_vir, cp_vir, zvir
!--------------------------------------------------------------
                                call timing_off('FV_PHYS')
      endif

!---- diagnostics for FV dynamics -----

                                call timing_on('FV_DIAG')

  call fv_diag(fv_time, nlon, mlat, nlev, beglat, endlat, ncnst, zvir,   &
               dt_atmos, .true.)

                                call timing_off('FV_DIAG')

!--------------------------------------------------------
!df 141010
if(id_t_surf > 0) used = send_data(id_t_surf, t_surf, Time)
if(id_flux_t > 0) used = send_data(id_flux_t, flux_t, Time)
if(id_flux_q > 0) used = send_data(id_flux_q, flux_q, Time)
if(id_flux_r > 0) used = send_data(id_flux_r, flux_r, Time)
if(id_conv_rain > 0) used = send_data(id_conv_rain, conv_rain, Time)
if(id_cond_rain > 0) used = send_data(id_cond_rain, cond_rain, Time)
if(id_pme       > 0) used = send_data(id_pme      , pme      , Time)
if(id_conv_rain_profile > 0) used = send_data(id_conv_rain_profile, conv_rain_profile, Time)
if(id_cond_rain_profile > 0) used = send_data(id_cond_rain_profile, cond_rain_profile, Time)

if(id_rh > 0)     used = send_data(id_rh, rh, Time)

if (id_OLRnu>0) used = send_data(id_OLRnu, OLRnu, Time)
if (id_OLR>0)   used = send_data(id_OLR, OLR, Time)
if (id_solar>0) used = send_data(id_solar, solar, Time)

!--------------------------------------------------------

 end subroutine atmosphere


 subroutine atmosphere_end

#include <fv_arrays.h>
    isiz = nlon/4
    nx = nlon/isiz
    tsiz = nx * ( endlat - beglat + 1 )         ! Total loop length
! Calling phys one latitude at a time using a big/fat OpenMP loop
! For cache performance, this could be changed to finer decomposition if needed
!$omp parallel do private(ij,i,j,k,m,is,ie,p_full,p_half)
    do ij=1,tsiz
       j  = beglat + (ij-1) / nx
       is = 1 + isiz * mod(ij-1, nx)
       ie = is + isiz - 1
       q(is:ie,j,nlev,2) = t_surf(is:ie,j)
    end do

!----- initialize domains for writing global physics data -----

    call set_domain ( fv_domain )
    call get_time (fv_time, seconds,  days)
    call write_fv_rst( 'RESTART/fv_rst.res', days, seconds, grav, &
         restart_format )


    call fv_end(days, seconds)
   
    deallocate(t_surf)
    deallocate(p_half, p_full)
    deallocate(flux_t, flux_q, flux_r)
    deallocate(conv_rain, cond_rain)
    deallocate(net_surf_sw_down, surf_lw_down)
    deallocate(conv_rain_profile, cond_rain_profile, rh)

    !call radiation_end

 end subroutine atmosphere_end

end module atmosphere_mod
