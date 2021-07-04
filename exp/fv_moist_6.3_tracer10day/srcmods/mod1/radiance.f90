module radiance_mod

! ==================================================================================
! index order is different from the main program
! 1-- surface now
! nLay-- TOA
! ==================================================================================

   use fms_mod,               only: open_file, check_nml_error, &
                                    mpp_pe, close_file, &
                                    error_mesg, FATAL

   use constants_mod,         only: stefan, cp_air, grav, pstd_mks, &
                                B_nu, &
                                c, c2, atm, Rstar, Tref, mpr
   use dimensions,            only : nGas, nAng, nS, nlinesMAX, &
                                nTem, datadir 
   use line_profs,            only : lorentz, doppler, voigt


   use    diag_manager_mod,   only: register_diag_field, send_data

   use    time_manager_mod,   only: time_type, &
                                    operator(+), operator(-), operator(/=)
   use      topography_mod,   only: get_ocean_mask
!   use      transforms_mod,   only: get_grid_boundaries
   use fv_pack, only : master
 
!==================================================================================
implicit none
private
!==================================================================================

! version information 
 
character(len=128) :: version='$Id: radiance.f90 $'
character(len=128) :: tag='homemade'

!==================================================================================

! public interfaces

public :: radiance_init , radiance_update ! 
!==================================================================================


! module variables
  logical :: initialized =.false.
  integer :: vgas = 1 
  real    :: solar_constant  = 1000. !1360.0
  real    :: del_sol         = 1.4
! modif omp: winter/summer hemisphere
  real    :: del_sw          = 0.
  real    :: kappa_gray = 0.01 !5.0d-3 

  integer  :: nlines_beguier = 7786 ! number of lines in Beguier+ JQRST (2015) CH4 data
  logical  :: use_beguier_2015 = .true. ! include Beguier+ JQRST (2015) CH4 data in near-IR?

  integer :: iGas_H2O = 1
  integer :: iGas_CO2 = 2
  integer :: iGas_O3  = 3
  integer :: iGas_CH4 = 6

    !------- inputs --------
  logical :: HITEMP = .False.                  ! use HITEMP line lists?
!  real :: deltanu_trunc = 25.0d0      ! truncation wavenumber [cm^-1]
  real :: dTlay = 50.0 !25.0
  real :: nu_lw1 = 1.0 !1.0
  real :: nu_lw2 = 2000.0 !3000.0
  logical, parameter :: gray_debug   = .false. ! uniform opacity vs. wavenumber

    ! constants once loaded
  integer, dimension(nGas) :: nlines, deltanu_trunc            ! actual number of lines
  real, dimension(nGas)    :: mu_i 
  integer, dimension(nGas, nlinesMAX) ::  mol,  iso    ! molecule number
  real, dimension(nGas, nlinesMAX) ::  nu0, Sref, einstein, gam_air, gam_self,&
          Egnd, ncoeff, delta , &
          gam_lore, gam_dopp, Strue

!  integer iLay, iLev, iAng, iS                    ! do loop variables
!  integer iRev, iGas, ierr

  real, dimension(nAng) :: cosa, ang_wt
  real, dimension(nS) ::  nu_lw, nu_sw, sigma_lw, &
          sigma_CIA  
  real :: dnu_lw

  real, allocatable, dimension(:,:) :: ss, cos_lat, p2, solar
  real, allocatable, dimension(:,:) :: T_k_grid !nlev, nTem
  integer, allocatable, dimension(:,:,:) :: iT1, iT2, iT1_out
  real, allocatable, dimension(:,:,:,:) ::  T_lin_weight !(nS,nLay) 
  real, allocatable, dimension(:,:,:,:) :: f_i !nGas
  real, allocatable, dimension(:,:,:)   :: dp, mu_avg
  real, allocatable, dimension(:,:,:,:) :: sigma_lw_ar 
  real, allocatable, dimension(:,:,:)   :: sigma_total_dry
  real, allocatable, dimension(:,:,:,:) :: sigma_t_i, &
          dtau_lw, dtau_lw_a, dTran, &
          Bnu_lay, Cterm, sigma_d_i, sigma_v_i
  real, allocatable, dimension(:,:,:,:) :: tau_lw, Bnu_lev, &
          I_lev_up, I_lev_dn
  real, allocatable, dimension(:,:,:,:,:) :: log_sig_d, log_sig_v 

real, allocatable, dimension(:,:,:) :: OLRnu, OLRnu_temp, OSRnu, OSRnu_temp, GSRnu, GSRnu_temp, Bnu_s, tau_lw_inf, &
Beff
real, allocatable, dimension(:,:) :: OLR_temp, OSR_temp, GSR_temp, OLR, OSR, GSR
real, allocatable, dimension(:,:,:) :: Ilev_up_lw, Ilev_dn_lw, Flev_up_lw, Flev_dn_lw, &
       Ilev_up_sw, Ilev_dn_sw, Flev_up_sw, Flev_dn_sw 
real, allocatable, dimension(:,:,:) :: dF_lw, dF_sw, dTdt_lw, dTdt_sw

real, save :: pi, deg_to_rad , rad_to_deg


!namelist/radiance_nml/ os08, constant_albedo, solar_constant, del_sol, 
!==================================================================================
!-------------------- diagnostics fields -------------------------------

integer :: id_flux_lw, id_tau_lw_inf, &
        id_sigma_lw, id_sigma_di, id_sigma_vi
logical :: used
character(len=10), parameter :: mod_name = 'lbl'

real :: missing_value = -999.

contains
! ==================================================================================
! ==================================================================================

subroutine radiance_init(is, ie, js, je, num_levels, axes, Time, lat, qlay, play,Tlay,psurf)
!-------------------------------------------------------------------------------------
integer, intent(in), dimension(4) :: axes
type(time_type), intent(in)       :: Time
integer, intent(in)               :: is, ie, js, je, num_levels
real, intent(in) , dimension(:,:)   :: lat, psurf
real, intent(in),  dimension(:,:,:) :: qlay, play, Tlay

!-------------------------------------------------------------------------------------
integer, dimension(3) :: half = (/1,2,4/)
integer :: ierr, io, unit, i
logical :: water_file_exists
!-----------------------------------------------------------------------------------------
! read namelist and copy to logfile

!unit = open_file ('input.nml', action='read')
!ierr=1
!do while (ierr /= 0)
!   read  (unit, nml=radiance_nml, iostat=io, end=10)
!   ierr = check_nml_error (io, 'radiance_nml')
!enddo
!10 call close_file (unit)

!unit = open_file ('logfile.out', action='append')
!if ( mpp_pe() == 0 ) then
!  write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
!  write (unit, nml=radiance_nml)
!endif
!call close_file (unit)

pi    = 4.0*atan(1.)
deg_to_rad = 2.*pi/360.
rad_to_deg = 360.0/2./pi

initialized = .true.
deltanu_trunc(:) = 25.0d0 

allocate(T_k_grid       (num_levels, nTem))
allocate(iT1		(ie-is+1, je-js+1, num_levels))
allocate(iT1_out        (ie-is+1, je-js+1, num_levels))
allocate(iT2            (ie-is+1, je-js+1, num_levels))
allocate(T_lin_weight   (ie-is+1, je-js+1, num_levels, nS))

allocate (sigma_lw_ar           (num_levels, nS, nGas, nTem))
allocate (sigma_total_dry       (num_levels, nS, nTem))

allocate (log_sig_d     (ie-is+1, je-js+1, num_levels, nS, 2))
allocate (log_sig_v     (ie-is+1, je-js+1, num_levels, nS, 2))

allocate (sigma_t_i     (ie-is+1, je-js+1, num_levels, nS))
allocate (sigma_d_i     (ie-is+1, je-js+1, num_levels, nS))
allocate (sigma_v_i     (ie-is+1, je-js+1, num_levels, nS))
allocate (dtau_lw       (ie-is+1, je-js+1, num_levels, nS))
allocate (dtau_lw_a     (ie-is+1, je-js+1, num_levels, nS))
allocate (dTran         (ie-is+1, je-js+1, num_levels, nS))
allocate (Bnu_lay       (ie-is+1, je-js+1, num_levels, nS))
allocate (Cterm         (ie-is+1, je-js+1, num_levels, nS))

allocate (tau_lw        (ie-is+1, je-js+1, num_levels+1, nS))
allocate (Bnu_lev       (ie-is+1, je-js+1, num_levels+1, nS))
allocate (I_lev_up      (ie-is+1, je-js+1, num_levels+1, nS))
allocate (I_lev_dn      (ie-is+1, je-js+1, num_levels+1, nS))

allocate (OLRnu         (ie-is+1, je-js+1, nS))
allocate (OLRnu_temp    (ie-is+1, je-js+1, nS))
allocate (OSRnu         (ie-is+1, je-js+1, nS))
allocate (OSRnu_temp    (ie-is+1, je-js+1, nS))
allocate (GSRnu         (ie-is+1, je-js+1, nS))
allocate (GSRnu_temp    (ie-is+1, je-js+1, nS))
allocate (Bnu_s         (ie-is+1, je-js+1, nS))
allocate (Beff          (ie-is+1, je-js+1, nS))
allocate (tau_lw_inf    (ie-is+1, je-js+1, nS))

allocate (OLR         (ie-is+1, je-js+1))
allocate (OLR_temp    (ie-is+1, je-js+1))
allocate (OSR         (ie-is+1, je-js+1))
allocate (OSR_temp    (ie-is+1, je-js+1))
allocate (GSR         (ie-is+1, je-js+1))
allocate (GSR_temp    (ie-is+1, je-js+1))

allocate (Ilev_up_lw    (ie-is+1, je-js+1, num_levels+1))
allocate (Ilev_dn_lw    (ie-is+1, je-js+1, num_levels+1))
allocate (Flev_up_lw    (ie-is+1, je-js+1, num_levels+1))
allocate (Flev_dn_lw    (ie-is+1, je-js+1, num_levels+1))
allocate (Ilev_up_sw    (ie-is+1, je-js+1, num_levels+1))
allocate (Ilev_dn_sw    (ie-is+1, je-js+1, num_levels+1))
allocate (Flev_up_sw    (ie-is+1, je-js+1, num_levels+1))
allocate (Flev_dn_sw    (ie-is+1, je-js+1, num_levels+1))

allocate (dF_lw    (ie-is+1, je-js+1, num_levels))
allocate (dF_sw    (ie-is+1, je-js+1, num_levels))
allocate (dTdt_lw  (ie-is+1, je-js+1, num_levels))
allocate (dTdt_sw  (ie-is+1, je-js+1, num_levels))
!allocate (Tcold    (ie-is+1, je-js+1, num_levels))
allocate (dp       (ie-is+1, je-js+1, num_levels))
allocate (mu_avg   (ie-is+1, je-js+1, num_levels))
allocate (f_i      (ie-is+1, je-js+1, num_levels, nGas))

allocate (ss               (ie-is+1, je-js+1))
allocate (cos_lat          (ie-is+1, je-js+1))
allocate (solar            (ie-is+1, je-js+1))
allocate (p2               (ie-is+1, je-js+1))

f_i(:,:,:,1) = qlay  !should be molar concentration, q is mass concentration now
f_i(:,:,:,2) = 400.e-6
mu_i(1) = 18.01528     ! molar mass of H2O [g/mol]
mu_i(2) = 44.0095      ! molar mass of CO2 [g/mol]
!-----------------------------------------------------------------------
    !-------- read line data --------
    ! line data is only read if calc_sigma(iGas) = .true.
    call read_line_data()

    !-------- set up temperature grid on which cross-sections are defined --------
    call setup_T_k_grid(Tlay(1, 1, :),T_k_grid)

!    if (master) write(6,*) T_k_grid(1,:), T_k_grid(30,:)
!    write(6,*) T_k_grid(num_levels,:) !play(:,:,1), Tlay(:,:,1)

    call define_cosa_quad(cosa,ang_wt)


    !-------- set up radiative transfer calculations --------
!    call setup_shortwave(Tlay,dTlay,play,ps,f_i,mu_i,grav,calc_sigma,nu_sw)

!    if (master) write(6,*) play(:,:,1), Tlay(:,:,1)
    call setup_longwave(T_k_grid,play,f_i,mu_i,nu_lw)


    !-------- define lower Tlay for sigma interpolation --------
!    Tcold = Tlay - dTlay


!    if (master) write(6,*) nu0

!------------ initialize diagnostic fields ---------------


id_tau_lw_inf  = register_diag_field(mod_name, 'tau_lw_inf',        &
     axes((/1,2,3/)), Time, 'longwave optical thickness','none') 

id_flux_lw   = register_diag_field(mod_name, 'flux_lw',        &
     axes((/1,2,4/)), Time, 'longwave flux(positive upward)','W/m**2') 

id_sigma_lw   = register_diag_field(mod_name, 'sigma_lw',        &
     axes((/3,5/)), Time, 'longwave abs','') 

id_sigma_di   = register_diag_field(mod_name, 'sigma_di',        &
     axes((/1,2,5/)), Time, 'longwave abs','') 

id_sigma_vi   = register_diag_field(mod_name, 'sigma_vi',        &
     axes((/1,2,5/)), Time, 'longwave abs','') 

!if (id_sigma_lw>0) used = send_data(id_sigma_lw, sigma_lw_ar(:,:,1,1) , Time)

return
end subroutine radiance_init


! ==================================================================================

subroutine radiance_update (is, js, Time_diag, lat, lon, qlay, Tlay, Tlev, Ts, ps, play,  &
                           plev, dTdt,net_surf_sw_down, surf_lw_down, &
                           sol, OL, OLnu)

integer, intent(in)                 :: is, js
type(time_type), intent(in)         :: Time_diag
real, intent(in) , dimension(:,:)   :: lat, lon, Ts, ps
real, intent(in) , dimension(:,:,:) :: Tlay, play, qlay
real, intent(in) , dimension(:,:,:) :: plev, Tlev
real, intent(out) , dimension(:,:,:) :: dTdt
real, intent(out) , dimension(:,:)  :: net_surf_sw_down, surf_lw_down, sol, OL
real, intent(out) , dimension(:,:,:) :: OLnu

integer :: iGas, iAng, iLev, iLay, nLay

nLay = size(Tlay,3)

ss  = sin(lat)
p2 = (1. - 3.*ss*ss)/4.
solar(:,:) = 0.25*solar_constant*(1.0 + del_sol*p2 + del_sw * ss)
!cos_lat  = cos(lat(1,:))
!   do i = 1, size(t,1)
!      solar(i,:) = 0.25*solar_constant*(1.0 + del_sol*p2 + del_sw * ss)
!      solar(i,:) = solar_constant*cos_lat/pi
!   enddo
!endif

dp(:,:,:)       = plev(:,:,1:nLay) - plev(:,:,2:nLay+1) !! the order has been reversed
f_i(:,:,:,1)    = qlay  !should be molar concentration, q is mass concentration now
!f_i(:,:,:,2)   = 4e-6 !should vary with global mean mass
mu_avg(:,:,:)   = 28.
do iGas=1,nGas
   mu_avg = mu_avg + (mu_i(iGas) - 28.) *f_i(:,:,:,iGas)
end do

    !-------- temperature grid interpolation calculation --------
    T_lin_weight = -1.0e8
    if(nTem>1) call calculate_T_k_grid(Tlay(:,:,:),iT1,T_lin_weight)

    !-------- visible radiative transfer calculation --------
    
    OSRnu(:,:,:)        = 0.0d0
    OSR                 = 0.0d0
    GSR                 = 0.0d0
    Flev_up_sw(:,:,:)   = 0.0d0
    Flev_dn_sw(:,:,:)   = 0.0d0

    !if(verbose)then
    !   write(*,*) ' Calculating shortwave direct beam attenuation.'
    !endif
    !sw_angle_loop: do iAng = 1, nAng ! integration over propagation angle
!
!       if(verbose)then
!          write(*,'(a44,f8.2,a9)') ' Calculating shortwave with emission angle = ', &
!               acos(cosa(iAng))*180.0d0/pi, ' degrees.'
!       endif
!       call calculate_shortwave(iT1,cosa(iAng),T_lin_weight,ps,dp,f_i,mu_avg,grav,kappa_gray, &
!            OSRnu_temp,OSR_temp,GSRnu_temp,GSR_temp,Ilev_up_sw,Flev_dn_sw)

       ! c.f. e.g. Stamnes et al. (2000), DISORT tech. report, eqn. (9a)
!       OSR        = OSR        + 2*pi*OSR_temp  *cosa(iAng)*ang_wt(iAng)
!       OSRnu      = OSRnu      + 2*pi*OSRnu_temp*cosa(iAng)*ang_wt(iAng)
!       Flev_up_sw = Flev_up_sw + 2*pi*Ilev_up_sw*cosa(iAng)*ang_wt(iAng)

       ! at the moment GSR does not depend on propagation angle as we are in the
       ! no atmospheric scattering regime. So no need to integrate it over cosa.
       !GSR        = GSR        + 2*pi*GSR_temp  *cosa(iAng)*ang_wt(iAng)
       !GSRnu      = GSRnu      + 2*pi*GSRnu_temp*cosa(iAng)*ang_wt(iAng)

!    end do sw_angle_loop
    ! these two are fluxes in units of W/m2 and W/m2/cm^-1 already
!    GSR        = GSR_temp
!    GSRnu      = GSRnu_temp

    ! absorbed stellar radiation = ISR - OSR [W/m2]
!    ASR = ISR - OSR

    !-------- IR radiative transfer calculation --------
    
    OLRnu(:,:,:)        = 0.0d0
    OLR                 = 0.0d0
    Flev_up_lw(:,:,:)   = 0.0d0
    Flev_dn_lw(:,:,:)   = 0.0d0

    lw_angle_loop: do iAng = 1, nAng ! integration over propagation angle

       !if(verbose)then
       !   write(*,'(a44,f8.2,a9)') ' Calculating longwave with emission angle = ', &
       !        acos(cosa(iAng))*180.0d0/pi, ' degrees.'
       !endif

       call calculate_longwave(iT1,cosa(iAng),Tlay,Tlev,Ts,T_lin_weight,play,dp,f_i,mu_avg, &
            kappa_gray,OLRnu_temp,OLR_temp,Ilev_up_lw,Ilev_dn_lw)

       OLR        = OLR        + 2*pi*OLR_temp  *cosa(iAng)*ang_wt(iAng)
       OLRnu      = OLRnu      + 2*pi*OLRnu_temp*cosa(iAng)*ang_wt(iAng)
       Flev_up_lw = Flev_up_lw + 2*pi*Ilev_up_lw*cosa(iAng)*ang_wt(iAng)
       Flev_dn_lw = Flev_dn_lw + 2*pi*Ilev_dn_lw*cosa(iAng)*ang_wt(iAng)

    end do lw_angle_loop

    !-------- heating rate calculation --------
    do ilay=1,nLay
       iLev = iLay + 1
       dF_lw(:,:,iLay) = (Flev_up_lw(:,:,iLev) - Flev_up_lw(:,:,iLev-1)) - &
       (Flev_dn_lw(:,:,iLev) - Flev_dn_lw(:,:,iLev-1))
       dF_sw(:,:,iLay) = (Flev_up_sw(:,:,iLev) - Flev_up_sw(:,:,iLev-1)) - &
       (Flev_dn_sw(:,:,iLev) - Flev_dn_sw(:,:,iLev-1))
    end do
    !dTdt_lw(:,:,:) = -(grav/cp_heat)*dF_lw(:,:,:)/dp(:,:,:)
    !dTdt_sw(:,:,:) = -(grav/cp_heat)*dF_sw(:,:,:)/dp(:,:,:)
    !dTdt(:)    = dTdt_lw(:) + dTdt_sw(:)
    dTdt(:,:,:)    = -dF_lw(:,:,:)

    net_surf_sw_down(:,:)       = solar
    surf_lw_down(:,:)           = Flev_dn_lw(:,:,1)
    sol         = solar
    OL          = OLR
    OLnu        = OLRnu

    !-------- report results --------
    !if(report_results)then
    !   write(*,*) '---------------------------------'
    !   write(*,*) 'Results:'
    !   write(*,'(a9,f14.7,a5)') '  ISR  = ',ISR,' W/m2'
    !   write(*,'(a9,f14.7,a5)') '  OLR  = ',OLR,' W/m2'
    !   write(*,'(a9,f14.7,a5)') '  ASR  = ',ASR,' W/m2'
    !   write(*,'(a9,f14.7)')    '  OLRe = ',OLR/(stefan*Ts**4)
    !   write(*,'(a9,f14.7)')    '  Apla = ',OSR/ISR
    !   write(*,'(a9,f14.7)')    '  Areq = ',1.0d0 - OLR/ISR
       ! Areq is the albedo required for thermal equilibrium
    !endif

!------- downward lw flux surface -------
!      if ( id_lwdn_sfc > 0 ) then
!          used = send_data ( id_lwdn_sfc, surf_lw_down, Time_diag)
!      endif

!OLR = Tlev(:,:,1)

!if (id_tau_lw_inf>0)    used = send_data(id_tau_lw_inf, &
!        f_i(:,:,:,1), Time_diag)
if (id_flux_lw>0)       used = send_data(id_flux_lw, &
        Flev_up_lw, Time_diag)

!if (id_tau_lw_inf>0)    used = send_data(id_tau_lw_inf, &
!        iT1, Time_diag)

!if (master) write(6,*) Tlay(1,1,:), Tlev(1,1,:)
if (id_sigma_di>0) used = send_data(id_sigma_di, &
        sigma_d_i(:,:,1,:) , Time_diag)
if (id_sigma_vi>0) used = send_data(id_sigma_vi, &
        sigma_v_i(:,:,1,:) , Time_diag)


return
end subroutine radiance_update

! ==================================================================================

                                                                      
subroutine radiance_end
                                                                                                      
!deallocate (b, tdt_rad, tdt_sw, entrop_rad) 


end subroutine radiance_end

! ==================================================================================


 subroutine read_line_data()

    ! read all line data
    integer :: kk, il, ierr, iGas

    integer :: mol_temp, iso_temp      ! molecule number
    real ::  nu0_temp, Sref_temp, einstein_temp,gamair_temp, gamself_temp,Egnd_temp,ncoeff_temp,  delta_temp

    ! read lines
    if (master) write(6,*) "creating line data..."

    molec_loop : do iGas = 1, nGas
!    if(calc_sigma(iGas))then
       if(iGas==iGas_H2O)then
          if(HITEMP)then
             open(unit=111,file=trim(datadir)//'01_HITEMP2010_red2.par')
          else
             open(unit=111,file=trim(datadir)//'01_hit12.par')
          endif
       elseif(iGas==iGas_CO2)then
          if(HITEMP)then
             open(unit=111,file=trim(datadir)//'02_HITEMP2010_red2.par')
          else
             open(unit=111,file=trim(datadir)//'02_hit12.par')
          endif
!       elseif(iGas==iGas_O3)then
!          open(unit=111,file=trim(datadir)//'03_hit12.par')
!       elseif(iGas==iGas_CH4)then
!          open(unit=111,file=trim(datadir)//'06_hit12.par')
!       elseif(iGas==iGas_SO2)then
!          open(unit=111,file=trim(datadir)//'09_hit12.par')
!       elseif(iGas==iGas_H2S)then
!          open(unit=111,file=trim(datadir)//'31_hit12.par')
       else
          write(6,*) 'Note: spectral dataset not found for iGas = ',iGas,' .'
       end if

       il = 1
       read_loop : do 
          if(HITEMP)then
             read(111,'(i2, i1, f12.6, 2e10.3, 2f5.4, f10.4, f4.2)',iostat=kk) mol_temp, iso_temp, &
                  nu0_temp, Sref_temp, einstein_temp, gamair_temp, gamself_temp, Egnd_temp, ncoeff_temp!, delta_temp
          else
             read(111,'(i2, i1, f12.6, 2e10.3, 2f5.4, f10.4, f4.2, f8.2)',iostat=kk) mol_temp, iso_temp, &
                  nu0_temp, Sref_temp, einstein_temp, gamair_temp, gamself_temp, Egnd_temp, ncoeff_temp, delta_temp
          endif

          if(Sref_temp>1.0d-25)then ! don't allow weak lines
          !if(Sref_temp>1.0d-25)then ! don't allow weak lines
          !if(Sref_temp>1.0d-30)then ! don't allow weak lines
          !if(Sref_temp>1.0d-35)then ! don't allow weak lines
             mol(iGas,il)      = mol_temp
             iso(iGas,il)      = iso_temp
             nu0(iGas,il)      = nu0_temp
             Sref(iGas,il)     = Sref_temp
             einstein(iGas,il) = einstein_temp
             gam_air(iGas,il)  = gamair_temp
             gam_self(iGas,il) = gamself_temp
             Egnd(iGas,il)     = Egnd_temp
             ncoeff(iGas,il)   = ncoeff_temp
             delta(iGas,il)    = delta_temp
             il = il + 1
          endif

          
          if(kk == -1)then
             exit read_loop
          endif
          if(il>nlinesMAX)then
             exit read_loop
          endif

       end do read_loop
       close(111)
       
       nlines(iGas) = il - 1

       ! special extra to add unassigned near-IR methane lines from
       ! Beguier+, JQRST (2015)
!       if(use_beguier_2015 .and. iGas==iGas_CH4)then
!          write(*,*) 'Adding Beguier+ (2015) unassigned CH4 lines in near-IR...'
!          open(unit=112,file=trim(datadir)//'beguier_2015.par')
!          do il=nlines(iGas)+1,nlines(iGas)+nlines_beguier
!             read(112,*) nu0_temp, Sref_temp
!             mol(iGas,il)      = 6
!             iso(iGas,il)      = 1
!             nu0(iGas,il)      = nu0_temp
!             Sref(iGas,il)     = Sref_temp
!             einstein(iGas,il) = -1.0e10
!             gam_air(iGas,il)  = 1.8d-2 ! from Beguier+ (2015) Section 3.1.
!             gam_self(iGas,il) = 0.0d0
!             Egnd(iGas,il)     = 1.0d0
!             ncoeff(iGas,il)   = 0.0d0 ! TO BE DETERMINED
!             delta(iGas,il)    = 0.0d0
!          end do
!          close(112)

!          nlines(iGas) = il - 1
!       end if


!    end if
    end do molec_loop

  end subroutine read_line_data

  ! ==================================================================================

  subroutine line_strength_width(mu_i,p,p_i,T,iGas)

    ! calculate line strengths, widths and frequency pressure shifts
    ! for a given species at a given temperature

    real, intent(in) :: T    ! temperature [K]
    real, intent(in) :: p    ! total pressure [Pa]
    real, intent(in) :: p_i  ! partial pressure [Pa]
    real, intent(in) :: mu_i ! molar mass [g/mol]
    integer, intent(in) :: iGas

    real :: Tfact, Tfact1       ! temperature-dependent factors []
    real :: QrefQ               ! total internal partition sum
    !real(8) gi                  ! state independent degeneracy factor
    
    integer :: il

    logical, parameter :: do_air_shift = .false. ! calculate pressure-shifted line center
    logical, parameter :: verbose      = .true. ! print more info for diagnostic

    if(T<20.0d0)then
       call error_mesg('radiance:','In line_strength_width T = exiting',FATAL)
    end if

!    if(iGas==iGas_H2 .or. iGas==iGas_N2)then ! todo: generalize this
!       if(verbose) write(*,*) 'Treating ', gas_name(iGas), ' as vib-rot inactive.'
!    else

!       if(verbose) write(*,*) "calculating line strengths for ",gas_name(iGas),"..."

       
       read_loop : do il = 1, nlines(iGas)
          
          ! calculate (air) pressure-shifted nu
!          if(do_air_shift)then
             !nu_shift(iGas,il) = nu0(iGas,il) + delta(il)*p_i(iGas)/atm
!          end if
          
          ! calculate Lorentz and Doppler linewidths [cm^-1]
          ! Lorentz linewidth is the Lorentz HWHM
          ! Doppler linewidth is the Doppler HWHM / sqrt(ln(2))
          gam_lore(iGas,il) = (Tref/T)**ncoeff(iGas,il) &
                               *(gam_air(iGas,il)*(p - p_i)/atm + gam_self(iGas,il)*p_i/atm)
          gam_dopp(iGas,il) = nu0(iGas,il)*sqrt(2*Rstar*T/(mu_i/1.0d3))/c
          ! note mu_i not mu_avg here because it is the mean speed of the molecule doing the
          ! absorbing that counts.
          
          ! get Total Internal Partition Sum (TIPS) ratio Qref / Q
          ! simple approach based on BD_TIPS data
          ! robust up to about 1000 K for H2O and CO2
          ! up to 500-600 K for CH4
          ! up to about 350 K for O3
          if(iGas==iGas_H2O .or. iGas==iGas_CO2 .or. iGas==iGas_O3 .or. iGas==iGas_CH4)then
             QrefQ = (Tref/T)**1.5d0
          else
             call error_mesg('radiance: ','Molecule Q(T) needs to be assessed still!',FATAL)
          end if

          
          ! get actual line intensity
          Tfact1          = (1.0d0 - exp(-c2*nu0(iGas,il)/T)) &
                              / (1.0d0 - exp(-c2*nu0(iGas,il)/Tref))
          Tfact           = exp( -c2*Egnd(iGas,il)*(1.0d0/T - 1.0d0/Tref) )
          Strue(iGas,il)  = Sref(iGas,il)*Tfact*Tfact1*QrefQ

          ! no temperature scaling for Beguier+ (2015) lines
          if(use_beguier_2015 .and. iGas.eq.iGas_CH4 .and. il>(nlines(iGas)-nlines_beguier)) Strue(iGas,il) = Sref(iGas,il)
          
       end do read_loop
       
       nlines(iGas) = il - 1

!    end if
    
  end subroutine line_strength_width

 
  ! ==================================================================================


  subroutine get_line_abs_cross_section(iGas,nu,sigma,T)
    
    ! calculate line absorption cross-section vs. nu for a given species
    ! in a given wavenumber range

    integer, intent(in)  :: iGas
    real, intent(in)  :: nu(nS)     ! wavenumber [cm^-1]
    real, intent(in)  :: T          ! temperature [K]: for CO2 chi-factor calculation
    real, intent(out) :: sigma(nS)  ! total absorption cross-section [cm2/molec]

    logical :: mask(nS)
    integer :: iS, il
    real :: f_temp, gam_temp
    
    real :: nu1  ! start wavenumber [cm^-1]
    real :: nu2  ! finish wavenumber [cm^-1]
    
    logical, parameter :: use_voigt   = .false.  ! use voigt function for line profiles
    logical, parameter :: debug       = .false. ! print additional text for diagnostic

!    if(debug) write(*,*) "calculating absorption cross sections for ",gas_name(iGas),"..."

    sigma(:) = 0.0d0

    nu1 = nu(1)
    nu2 = nu(nS)

!    if(.not.iGas==iGas_H2)then ! todo: transfer update from generate_kmatrix

       ! add the lines, one by one
       ! for now we use all isotopes (with Earth-like abundances!)
       line_loop : do il = 1, nlines(iGas)

          ! this is the very slow part
          trunc : if(nu0(iGas,il)>nu1-deltanu_trunc(iGas) .and. nu0(iGas,il)<nu2+deltanu_trunc(iGas))then

             mask = abs(nu - nu0(iGas,il)) < deltanu_trunc(iGas)

             write_sig_loop : do iS = 1, nS

                ! update sigma only when mask = T
                if(mask(iS))then

                   if(iGas==iGas_CO2)then
                      gam_temp = gam_lore(iGas,il) *chi_factor(nu(iS),nu0(iGas,il),T)
                      ! check but I believe this is correct
                   else
                      gam_temp = gam_lore(iGas,il) ! no sublorentzian lineshift
                   end if

                   if(use_voigt)then
                      call voigt(gam_temp, gam_dopp(iGas,il), nu(iS)-nu0(iGas,il), f_temp)
                   else
                      call lorentz(nu(iS)-nu0(iGas,il),gam_temp,f_temp)
                      !call doppler(nu(iS)-nu0(iGas,il),gam_dopp(iGas,il),f_temp)
                   endif

                   if(debug)then
                      call voigt(0.0d0,1.0d0,0.0d0,f_temp)
                      ! f_V(0) when Lorentz HWHM = 0: should yield f_D(0) = 1/sqrt(pi)
                      print*,f_temp
                      call voigt(1.0d0,1.0d-8,0.0d0,f_temp)
                      ! f_V(0) when Doppler HWHM -> 0: should yield f_L(0) = 1/pi
                      print*,f_temp
                      stop
                   endif
                   
                   sigma(iS) = sigma(iS) + Strue(iGas,il)*f_temp
                endif
             end do write_sig_loop
          end if trunc
       end do line_loop

!    end if
    
  end subroutine get_line_abs_cross_section

  ! ==================================================================================

  real(8) function chi_factor(nu_temp,nu_L,T)

    ! calculate sublorentzian line profile chi factor using Perrin & Hartman (1989) assumptions

    ! to be double-checked

    implicit none

    real(8), intent(in) :: T       ! temperature [K]
    real(8), intent(in) :: nu_L    ! line center wavenumber [cm^-1]
    real(8), intent(in) :: nu_temp ! wavenumber [cm^-1]

    ! parameters for the empirical scheme
    real(8), parameter, dimension(3) :: alpha = [0.0888d0,  0.0d0,     0.0232d0]
    real(8), parameter, dimension(3) :: beta  = [-0.160d0,  0.0526d0,  0.0d0]
    real(8), parameter, dimension(3) :: epsi  = [0.00410d0, 0.00152d0, 0.0d0]
    real(8), parameter, dimension(3) :: sigPH = [3.0,       30.0,      120.0] ! cm^-1

    real(8), dimension(3) :: B
    real(8) deltanu ! wavenumber separation [cm^-1]
    
    B          = alpha + beta*exp(-epsi*T)
    chi_factor = 1.0d0;

    deltanu   = abs(nu_temp - nu_L)
    if(deltanu<sigPH(1))then
       chi_factor = 1.0d0
    elseif(deltanu>sigPH(1) .and. deltanu<sigPH(2))then
       chi_factor = exp(-B(1)*(deltanu - sigPH(1)))
    elseif(deltanu>sigPH(2) .and. deltanu<sigPH(3))then
       chi_factor = exp(-B(1)*(sigPH(2) - sigPH(1))-B(2)*(deltanu - sigPH(2)))
    else
       chi_factor = exp(-B(1)*(sigPH(2) - sigPH(1)) - B(2)*(sigPH(3) - sigPH(2)) &
            - B(3)*(deltanu - sigPH(3)))
    end if

    return
  end function chi_factor

  ! ==================================================================================
  

  subroutine setup_longwave(T_k_grid,play,f_i,mu_i,nu_lw_out)

    real, intent(in)  :: T_k_grid(:,:)        ! temperature in layer  [K]
    real, intent(in)  :: play(:,:,:)        ! pressure layers [Pa]
    real, intent(in)  :: f_i(:,:,:,:)    ! species molar concentration [mol/mol]
    real, intent(in)  :: mu_i(nGas)        ! molar mass of species [g/mol]
    real, intent(out) :: nu_lw_out(nS)     ! longwave wavenumber [cm^-1]

    real :: nu_lw_check                       ! longwave wavenumber [cm^-1]
    integer ::  iGas, iLay, iS, iTem
    !------- subroutine options --------------  
    logical, parameter   :: verbose = .true.

 !   ! read namelist 
 !   open(10,file='input.nml',status='old',form='formatted',iostat=ierr)
 !   if (ierr/=0) then
 !      print*, 'Cannot find required input.nml file, aborting.'
 !      call abort
 !   else
 !      read(10,longwave_nml)
 !      close(10)
 !   endif
    
    ! create fine spectral grid and initialize abs array
    nu_lw(1)  = nu_lw1
    dnu_lw = (nu_lw2-nu_lw1)/dble(nS-1)  ! Wavenumber interval [cm^-1]
    do iS = 2, nS
       nu_lw(iS) = nu_lw(iS-1) + dnu_lw
    end do
    nu_lw_out = nu_lw

    ! calculate absorption cross section at each level

    sigma_lw(:)                = 0.0d0
    tau_lw_inf(:,:,:)          = 0.0d0
    sigma_lw_ar(:,:,:,:)       = 0.0d0
    sigma_total_dry(:,:,:)     = 0.0d0
    sigma_t_i(:,:,:,:)         = 0.0d0
    dtau_lw(:,:,:,:)           = +1.0d-8 ! initialize non-zero to avoid NaN when optical depth is low

    gray_debug_if : if(gray_debug)then

       ! don't calculate anything in this case

    else

       big_species_loop : do iGas = 1, nGas
          
          ! read gas absorption cross section data from file if requested
!          if(.not.calc_sigma(iGas))then

             ! no need to check the number of temperatures for each layer match
             ! this was done in setup_shortwave
             
!             ! check the spectral arrays match
!             open(unit=2,file='saved_sigma_data/nu_lw.dat')
!             do iS = 1, nS
!                read(2,*) nu_lw_check
!                if(abs(nu_lw_check-nu_lw(iS))>1.0d-8)then
!                   write(*,*) 'Error: nu_lw data in saved_sigma_data does not match!'
!                   print*,nu_lw_check,' vs. ',nu_lw(iS)
!                   stop
!                end if
!             end do
!             close(2)

             ! no need to check the pressure layers match
             ! this was done in setup_shortwave
             
!             open(unit=3,file='saved_sigma_data/sigma_'//gas_name(iGas)//'_lw.dat')
!             do iS = 1, nS
!                do iLay = 1, nLay
!                   do iTem = 1,nTem
!                      read(3,'(e12.4)') sigma_lw_ar(iS,iLay,iGas,iTem)
!                   end do
!                end do
!             end do
!             close(3)
             
!          else

             ! calculate cross-sections from scratch
             !   write(6,*) T_k_grid(1,:), T_k_grid(30,:)
             do iLay = 1, size(T_k_grid,1) !nLay
!                if(verbose)then
!                   print*,'For longwave at layer ',iLay,': '
!                endif
                
!                if(nTem==2)then
!                   call line_strength_width(mu_i(iGas),play(iLay),f_i(iLay,iGas)*play(iLay),Tlay(iLay)-dTlay,iGas)
!                   call get_line_abs_cross_section(iGas,nu_lw,sigma_lw,Tlay(iLay)+dTlay)
!                   sigma_lw_ar(:,iLay,iGas,1)    = sigma_lw
                
!                   call line_strength_width(mu_i(iGas),play(iLay),f_i(iLay,iGas)*play(iLay),Tlay(iLay)+dTlay,iGas)
!                   call get_line_abs_cross_section(iGas,nu_lw,sigma_lw,Tlay(iLay)-dTlay)
!                   sigma_lw_ar(:,iLay,iGas,nTem) = sigma_lw
!                else               
                do iTem=1,nTem

                   call line_strength_width(mu_i(iGas),play(1,1,iLay),f_i(1,1,iLay,iGas)*play(1,1,iLay),T_k_grid(iLay,iTem),iGas)
                   call get_line_abs_cross_section(iGas,nu_lw,sigma_lw,T_k_grid(iLay,iTem))
                   sigma_lw_ar(iLay,:,iGas,iTem)    = sigma_lw ! cross-section per molecule of species
!                end if


          ! multiply the species cross-section by molar concentration to get the
          ! contribution to the total cross-section, in cm2/molecule of air.
          ! at this point we also add the CIA and compute the total 'dry' cross-section.
          ! do not add variable gas cross-section to the total yet.
          not_vgas : if( iGas .ne. vgas) then
!             if ((f_i(i,j,iLay,iGas)*play(i,j,iLay)>1.0d1)) then
! for CIA calculation
!                   sigma_CIA(:) = 0.0d0

!                   if(nTem==2)then
!                      call get_CIA_cross_section(iGas,nu_lw,sigma_CIA,Tlay(iLay)-dTlay,play(iLay),f_i(iLay,:))
!                      sigma_total_dry(:,iLay,1)    = sigma_total_dry(:,iLay,1)    + &
!                           sigma_lw_ar(:,iLay,iGas,1)*f_i(iLay,iGas)    + sigma_CIA                     
!                      call get_CIA_cross_section(iGas,nu_lw,sigma_CIA,Tlay(iLay)+dTlay,play(iLay),f_i(iLay,:))
!                      sigma_total_dry(:,iLay,nTem) = sigma_total_dry(:,iLay,nTem) + &
!                           sigma_lw_ar(:,iLay,iGas,nTem)*f_i(iLay,iGas) + sigma_CIA
!                   else
!                      call get_CIA_cross_section(iGas,nu_lw,sigma_CIA,Tlay(iLay),play(iLay),f_i(iLay,:))
                      ! note that f_i is sent as a nGas array to get_CIA_cross_section even when we are inside
                      ! an iGas loop. this is because we need to know the abundance of the partner gas.
                      sigma_total_dry(iLay,:,iTem) = sigma_total_dry(iLay,:,iTem) + &
                      sigma_lw_ar(iLay,:,iGas,iTem)*f_i(1,1,iLay,iGas) !+ sigma_CIA
!                   end if
                
!             else
!                   sigma_total_dry(:,iLay,1) = sigma_total_dry(:,iLay,1) + sigma_lw_ar(:,iLay,iGas,1)*f_i(iLay,iGas)
!                   if(nTem==2)then
!                      sigma_total_dry(:,iLay,nTem) = sigma_total_dry(:,iLay,nTem) + sigma_lw_ar(:,iLay,iGas,nTem)*f_i(iLay,iGas)
!                   end if
!             end if
          end if not_vgas
                
                   
          end do
          end do

             ! save sigma data for future use
             ! save nu_lw to allow consistency checks
!             open(unit=1,file='saved_sigma_data/sigma_'//gas_name(iGas)//'_lw.dat')
!             open(unit=2,file='saved_sigma_data/nu_lw.dat')
!             do iS = 1, nS
!                do iLay = 1, nLay
!                   do iTem = 1,nTem
!                      write(1,'(e12.4)') sigma_lw_ar(iS,iLay,iGas,iTem)
!                   end do
!                end do
!                write(2,*) nu_lw(iS)
!             end do
!             close(1)
!             close(2)

!          end if
         
       end do big_species_loop
    end if gray_debug_if

  end subroutine setup_longwave

  ! ==================================================================================

  ! ==================================================================================

  subroutine define_cosa_quad(cosa,ang_wt)

    ! calculate the gaussian nodes and weights for
    ! the mean emission angle cosine quadrature

    ! c.f. https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss.E2.80.93Legendre_quadrature

    implicit none

    real(8), intent(inout) :: cosa(nAng)
    real(8), intent(inout) :: ang_wt(nAng)

    ! gaussian quadrature nodes and weights
    real(8) xi(nAng), ci(nAng)

    ! gauss interval (cosa = 0 to 1 here)
    real(8), parameter :: a = 0.0d0
    real(8), parameter :: b = 1.0d0

    ! move to dynamic allocation eventually

    if(nAng == 1)then
       xi(1) = 0.0d0
       ci(1) = 2.0d0
    elseif(nAng == 2)then
       xi(1) = +sqrt(1.0d0/3.0d0)
       xi(2) = -sqrt(1.0d0/3.0d0)
       ci(:) = 1.0d0
    elseif(nAng == 4)then
       xi(1) = -sqrt(3./7. + (2./7.)*sqrt(6./5.))
       xi(2) = -sqrt(3./7. - (2./7.)*sqrt(6./5.))
       xi(3) = +sqrt(3./7. - (2./7.)*sqrt(6./5.))
       xi(4) = +sqrt(3./7. + (2./7.)*sqrt(6./5.))
       ci(1) = (18.0d0-sqrt(30.0d0))/36.0d0
       ci(2) = (18.0d0+sqrt(30.0d0))/36.0d0
       ci(3) = (18.0d0+sqrt(30.0d0))/36.0d0
       ci(4) = (18.0d0-sqrt(30.0d0))/36.0d0

!!$ elseif(nAng == 8)then
!!$       xi(:) = (/-0.1834346424956498, &
!!$            0.1834346424956498,  &
!!$            -0.5255324099163290, &
!!$            0.5255324099163290,  & 
!!$            -0.7966664774136267, &
!!$            0.7966664774136267,  &
!!$            -0.9602898564975363, &
!!$            0.9602898564975363 /)
!!$       ci(:) =(/ 0.3626837833783620, & 
!!$            0.3626837833783620, &
!!$            0.3137066458778873, &
!!$            0.3137066458778873, &
!!$            0.2223810344533745, &
!!$            0.2223810344533745, &
!!$            0.1012285362903763, &
!!$            0.1012285362903763 /)
    else
       write(*,*) 'Error, emission angle weighting not defined for nAng = ', nAng
    endif

    cosa   = ((b-a)*xi + (b+a))/2
    ang_wt = ci*(b-a)/2
    
    return
  end subroutine define_cosa_quad

  ! ==================================================================================

  subroutine setup_T_k_grid(Tlay,T_k_grid_out)

    real(8), intent(in)  :: Tlay(:)              ! temperature in layer  [K]
    real(8), intent(out) :: T_k_grid_out(:,:) ! temperature array for cross-sections [K]

    integer ierr, iLay, iTem, nLay
    
    nLay = size(Tlay)

    ! read namelist 
!    open(10,file='input.nml',status='old',form='formatted',iostat=ierr)
!    if (ierr/=0) then
!       print*, 'Cannot find required input.nml file, aborting.'
!       call abort
!    else
!       read(10,temperature_k_grid_nml)
!       close(10)
!    endif
    
    ! calculate temperature grid on which cross-sections are defined
    do iLay = 1,nLay
       do iTem = 1,nTem
          T_k_grid_out(iLay,iTem) = Tlay(iLay) + dble(iTem-1-(nTem-1.)/2)*dTlay
       end do
    end do

    ! default values for iT1 and iT2
    iT1(:,:,:) = 1
    iT2(:,:,:) = iT1(:,:,:) + 1

    !T_k_grid_out = T_k_grid

  end subroutine setup_T_k_grid

  ! ==================================================================================

  subroutine calculate_T_k_grid(Tlay,iT1,T_lin_weight)

    real(8), intent(in)  :: Tlay(:,:,:)            ! temperature in layer  [K]
    integer, intent(inout) :: iT1(:,:,:)         ! T-grid array points for linear interpolation [] 
    real(8), intent(out) :: T_lin_weight(:,:,:,:) ! temperature weigting for linear interpolation [] 

    integer :: i, j, iLay, nLay
    
    nLay = size(Tlay,3)

    iT1_out(:,:,:) = -1
    
    ! either just select 1st values in array, or
    ! interpolate over T grid to get cross-section
    ! at actual atmospheric temperature

    ! find lower / upper T values at each layer
do i=1,size(Tlay,1)
   do j=1,size(Tlay,2)

    do iLay = 1,nLay

       find_iT : do 

          if(Tlay(i,j,iLay)<T_k_grid(iLay,iT1(i,j,iLay)))then
             ! move to lower T_k_grid values
             iT1(i,j,iLay) = iT1(i,j,iLay) - 1
             ! halt program if Tlay is below lowest T_k_grid value
             if(iT1(i,j,iLay)==0)then
                write(*,*) 'Temperature at iLay = ',iLay,' too cold.'
                write(*,*) 'Tlay(iLay)          = ',Tlay(i,j,iLay),' K.'
                write(*,*) 'min T_k_grid        = ',T_k_grid(iLay,1),' K.'
                call error_mesg('radiance: ','too cold for interp',FATAL)
             end if
          elseif(Tlay(i,j,iLay)>T_k_grid(iLay,iT1(i,j,iLay)+1))then
             ! move to higher T_k_grid values
             iT1(i,j,iLay) = iT1(i,j,iLay) + 1
             ! halt program if Tlay is above highest T_k_grid value
             if(iT1(i,j,iLay)==nTem)then
                !write(*,*) 'Temperature at iLay = ',iLay,' too hot.'
                !write(*,*) 'Tlay(iLay)          = ',Tlay(iLay),' K.'
                !write(*,*) 'max T_k_grid        = ',T_k_grid(iLay,nTem),' K.'
                call error_mesg('radiance: ','too hot for interp',FATAL)
             end if
          else
             ! T_k_grid values are correct, exit
             exit
          end if

       end do find_iT

       T_lin_weight(i,j,iLay,:) = (Tlay(i,j,iLay) - T_k_grid(iLay,iT1(i,j,iLay)))/(T_k_grid(iLay,iT1(i,j,iLay)+1) - T_k_grid(iLay,iT1(i,j,iLay))) 

    end do
   end do
end do
    iT2(:,:,:) = iT1(:,:,:) + 1

!    do iLay = 1,nLay
!       T_lin_weight(:,iLay) = (Tlay(iLay) - T_k_grid(iLay,iT1(iLay)))/(T_k_grid(iLay,iT2(iLay)) - T_k_grid(iLay,iT1(iLay))) 
!    end do
    !if(verbose)then
    !   do iLay = 1,nLay
    !      write(*,'(a15,i3)')    'iT1          = ', iT1(iLay)
    !   end do
    !   do iLay = 1,nLay
    !      write(*,'(a15,f8.3)')  'Tlay         = ', Tlay(iLay)
    !   end do
    !   do iLay = 1,nLay
    !      write(*,'(a15,f8.3)')  'T_lin_weight = ', T_lin_weight(1,iLay)
    !   end do
    !end if

    iT1_out(:,:,:) = iT1

  end subroutine calculate_T_k_grid

  ! ==================================================================================

  subroutine calculate_longwave(iT1,cosa_i,Tlay,Tlev,Ts,T_lin_weight,play,dp,f_i,mu_avg, &
       kappa_gray,OLRnu,OLR,I_lev_up_int,I_lev_dn_int)

    !use cross_section,  only : get_CIA_cross_section

    implicit none

    integer, intent(in)  :: iT1(:,:,:)             ! T-grid array points for linear interpolation [] 
    real(8), intent(in)  :: cosa_i                ! emission angle cosine []
    real(8), intent(in)  :: Tlay(:,:,:)            ! temperature in layer  [K]
    real(8), intent(in)  :: Tlev(:,:,:)            ! temperature at level  [K]
    real(8), intent(in)  :: Ts(:,:)                    ! surface temperature [K]
    real(8), intent(in)  :: T_lin_weight(:,:,:,:) ! temperature weigting for linear interpolation [] 
    real(8), intent(in)  :: play(:,:,:)            ! pressure layers [Pa]
    real(8), intent(in)  :: dp(:,:,:)              ! pressure difference [Pa]
    real(8), intent(in)  :: f_i(:,:,:,:)        ! species molar concentration [mol/mol]
    real(8), intent(in)  :: mu_avg(:,:,:)          ! average molar mass in layer [g/mol]
    real(8), intent(in)  :: kappa_gray            ! gray gas mass absorption cross-section for debug [m2/kg]

    real(8), intent(out) :: OLRnu(:,:,:)             ! outgoing longwave spectral irradiance [W/m2/cm^-1/sr] 
    real(8), intent(out) :: OLR(:,:)                   ! total outgoing longwave irradiance [W/m2/sr] 
    real(8), intent(out) :: I_lev_up_int(:,:,:)    ! upward irradiance in each layer [W/m2/sr]
    real(8), intent(out) :: I_lev_dn_int(:,:,:)    ! downward irradiance in each layer [W/m2/sr]

    integer :: i, j, iLay, nLay, iLev, iS, iRev

    nLay = size(Tlay,3)

    sigma_d_i(:,:,:,:) = 0.0d0
    sigma_v_i(:,:,:,:) = 0.0d0
    dtau_lw(:,:,:,:)   = +1.0d-8 ! initialize non-zero to avoid NaN when optical depth is low
    
    ! ----- beginning of part that requires nTem --------------

    gray_debug_loop : if(gray_debug)then

       ! nTem = 1 for this case verified in setup_shortwave
       do i=1,size(Tlay,1)
       do j=1,size(Tlay,2)
       do iLay = 1,nLay
          dtau_lw(i,j,iLay,:) = dtau_lw(i,j,iLay,:)+kappa_gray*dp(i,j,iLay)/grav * &
          f_i(i,j,iLay,1) !H2O is gray
       end do
       end do
       end do

    else

       ! either just select 1st values in array, or
       ! interpolate over T grid to get cross-section
       ! at actual atmospheric temperature
       if(nTem==1)then   
          !sigma_d_i = sigma_total_dry(:,:,1)
          !if(vgas.ne.0) sigma_v_i = sigma_lw_ar(:,:,vgas,1)
       else
          iT2(:,:,:) = iT1(:,:,:) + 1

          ! do log-space linear interpolation
          ! calculate logarithm of cross-section arrays
          ! y = y1 + (y2-y1)*(x-x1)/(x2-x1)
          do i=1,size(Tlay,1)
          do j=1,size(Tlay,2)
          do iLay = 1,nLay
             log_sig_d(i,j,iLay,:,:) = log10(sigma_total_dry(iLay,:,iT1(i,j,iLay):iT2(i,j,iLay)) + 1.0d-50)
             if(vgas.ne.0) log_sig_v(i,j,iLay,:,:) = &
             log10(sigma_lw_ar(iLay,:,vgas,iT1(i,j,iLay):iT2(i,j,iLay)) + 1.0d-50)
          end do
          end do
          end do
          sigma_d_i(:,:,:,:) = 10.0d0**(log_sig_d(:,:,:,:,1) + &
          (log_sig_d(:,:,:,:,2) - log_sig_d(:,:,:,:,1))*T_lin_weight(:,:,:,:))
          if(vgas.ne.0) sigma_v_i(:,:,:,:) = 10.0d0**(log_sig_v(:,:,:,:,1) + &
          (log_sig_v(:,:,:,:,2) - log_sig_v(:,:,:,:,1))*T_lin_weight(:,:,:,:))
       end if
       
    
       ! get sigma_CIA for variable gas (continuum absorption if H2O)
       ! and add variable gas cross-section (single-molecule + CIA) to the total
       ! for speed, we assume that _only_ the mixing ratio of the volatile gas
       ! 'iGas_var' may change at each timestep.
       if(vgas==0)then
          sigma_t_i = sigma_d_i
       else
           do i=1,size(Tlay,1)
           do j=1,size(Tlay,2)
           do iLay = 1, nLay

             !sigma_CIA(:) = 0.0d0

             ! calculate variable gas continuum at exact instantaneous atmospheric temperature, so no T interpolation
             !call get_CIA_cross_section(vgas,nu_lw,sigma_CIA,Tlay(iLay),play(iLay),f_i(iLay,:))
             sigma_t_i(i,j,iLay,:) = sigma_d_i(i,j,iLay,:) + sigma_v_i(i,j,iLay,:)*f_i(i,j,iLay,vgas)! + sigma_CIA
          end do
          end do
          end do
       endif

       ! compute total layer optical depth
       ! note array flip [nS,nLay] ==> [nLay,nS]
       do i=1,size(Tlay,1)
       do j=1,size(Tlay,2)
       do iLay = 1, nLay
          dtau_lw(i,j,iLay,:) = dtau_lw(i,j,iLay,:) + &
          1.0d-4*sigma_t_i(i,j,iLay,:)*dp(i,j,iLay)/(grav*mu_avg(i,j,iLay)*mpr)
       end do
       end do
       end do

    end if gray_debug_loop
    
    ! calculate Planck function vs. height and wavenumber
    do iS = 1, nS
       do i=1,size(Tlay,1)
       do j=1,size(Tlay,2)
       do iLev = 1, nLay+1
          call B_nu(nu_lw(iS),Tlev(i,j,iLev),Bnu_lev(i,j,iLev,iS))
       end do
       do iLay = 1, nLay
          call B_nu(nu_lw(iS),Tlay(i,j,iLay),Bnu_lay(i,j,iLay,iS))
       end do
       call B_nu(nu_lw(iS),Ts(i,j),Bnu_s(i,j,iS))
       end do
       end do
    end do
    
    ! calculate irradiance boundary conditions
    I_lev_up(:,:,1,:)    = Bnu_s(:,:,:)
    I_lev_dn(:,:,nLay+1,:) = 0.0d0

    ! scale dtau_lw by emission angle cosine
    dtau_lw_a(:,:,:,:) = dtau_lw(:,:,:,:)/cosa_i

    ! ----- end of part that requires nTem --------------

    ! vertical path optical depth at TOA
    ! somewhat inefficient to calculate it for every propagation angle
    tau_lw_inf(:,:,:) = sum(dtau_lw(:,:,:,:),3)
    
    ! calculate dTran
    ! the transmission through each layer,
    ! individually, for a given propagation angle.
    dTran = exp(-dtau_lw_a)

    ! make sure transmission does not cause floating point exception
     do i=1,size(Tlay,1)
     do j=1,size(Tlay,2)
     do iLay = 1,nLay
       do iS = 1,nS
          if(dTran(i,j,iLay,iS) < 1.0d-16)then           
             dTran(i,j,iLay,iS) = 1.0d-16
          end if
       end do
    end do
    end do
    end do

    ! calculate Cterm (not dependent on flux direction)
    Cterm(:,:,:,:) = 1.0d0/dtau_lw_a(:,:,:,:) - dTran(:,:,:,:)/(1.0d0 - dTran(:,:,:,:))

    ! Calculate upwards spectral irradiance using Clough et al. (1992) method
    ! starting from the BOA up
    do iLev = 2, nLay+1
       iLay             = iLev - 1
       Beff(:,:,:)      = Bnu_lev(:,:,iLev,:) + &
       2.0d0*(Bnu_lay(:,:,iLay,:) - Bnu_lev(:,:,iLev,:))*Cterm(:,:,iLay,:)
       I_lev_up(:,:,iLev,:) = I_lev_up(:,:,iLev-1,:)*dTran(:,:,iLay,:) + &
       Beff*(1.0d0 - dTran(:,:,iLay,:))
    end do
    ! Calculate downwards spectral irradiance using Clough et al. (1992) method
    ! starting from the TOA down
    do iRev = 1,nLay
       iLev             = nLay+1 - iRev ! start at nLev-1, finish at 1
       iLay             = iLev
       Beff             = Bnu_lev(:,:,iLev,:) + &
       2.0d0*(Bnu_lay(:,:,iLay,:) - Bnu_lev(:,:,iLev,:))*Cterm(:,:,iLay,:)
       I_lev_dn(:,:,iLev,:) = I_lev_dn(:,:,iLev+1,:)*dTran(:,:,iLay,:) + &
       Beff*(1.0d0 - dTran(:,:,iLay,:))
    end do

    ! calculate OLR and irradiances for output
    ! trapezoidal quadrature for the irradiances
    OLRnu(:,:,:)        = I_lev_up(:,:,nLay+1,:) 
    I_lev_up_int(:,:,:) = ((nu_lw(2)-nu_lw(1))/2.0d0) * &
    (I_lev_up(:,:,:,1) + I_lev_up(:,:,:,nS) + 2.0d0*sum(I_lev_up(:,:,:,2:nS-1),4) )
    I_lev_dn_int(:,:,:) = ((nu_lw(2)-nu_lw(1))/2.0d0) * &
    (I_lev_dn(:,:,:,1) + I_lev_dn(:,:,:,nS) + 2.0d0*sum(I_lev_dn(:,:,:,2:nS-1),4) )

!    I_lev_up_int(:,:,:) = ((nu_lw(2)-nu_lw(1))/2.0d0) * &
!    (I_lev_up(:,:,:,2) + I_lev_up(:,:,:,nS-1) + 2.0d0*sum(I_lev_up(:,:,:,2:nS-1),4) )
!    I_lev_dn_int(:,:,:) = ((nu_lw(2)-nu_lw(1))/2.0d0) * &
!    (I_lev_dn(:,:,:,2) + I_lev_dn(:,:,:,nS-1) + 2.0d0*sum(I_lev_dn(:,:,:,2:nS-1),4) )


    OLR(:,:) = I_lev_up_int(:,:,nLay+1)
    !OLR = 0.0d0
    !do iS = 1, nS
    !   OLR = OLR + OLRnu(iS)
    !end do
    !OLR = OLR*dnu_lw

    !-------- end longwave radiative transfer calculation --------

  end subroutine calculate_longwave

end module radiance_mod




