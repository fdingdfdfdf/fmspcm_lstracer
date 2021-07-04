module radiation_mod

! ==================================================================================
! ==================================================================================

   use fms_mod,               only: open_file, check_nml_error, &
                                    mpp_pe, close_file

   use constants_mod,         only: stefan, cp_air, grav, pstd_mks

   use    diag_manager_mod,   only: register_diag_field, send_data

   use    time_manager_mod,   only: time_type, &
                                    operator(+), operator(-), operator(/=)
   use      topography_mod,   only: get_ocean_mask
!   use      transforms_mod,   only: get_grid_boundaries
 
!==================================================================================
implicit none
private
!==================================================================================

! version information 
 
character(len=128) :: version='$Id: radiation.f90 $'
character(len=128) :: tag='homemade'

!==================================================================================

! public interfaces

public :: radiation_init, radiation_down, radiation_up, radiation_end
!==================================================================================


! module variables
logical :: initialized =.false.

logical :: os08 = .true.
logical :: moist_tau = .false.
logical :: constant_albedo = .true.

real    :: solar_constant  = 1360.0
real    :: del_sol         = 1.4
! modif omp: winter/summer hemisphere
real    :: del_sw          = 0.0
real    :: ir_tau_eq       = 6.0
real    :: ir_tau_pole     = 1.5
real    :: atm_abs         = 0.0
real    :: sw_diff         = 0.0
real    :: linear_tau      = 0.1
real    :: albedo_value    = 0.06
real    :: albedo_eq       = 0.2
real    :: albedo_pole     = 0.6
real    :: albedo_land     = 0.1
real    :: window          = 0.0 ! spectral window transparent to LW
real    :: wv_exponent     = 4.0 
real    :: solar_exponent  = 4.0 
real    :: wv_tau          = 1.0
real    :: reference_slp   = 1.e5 
logical :: tidally_locked  =.false.
logical :: do_column_wv    = .true.
real    :: noon_longitude  = 270.0
real    :: kappa_v         = 0.01 !absorption coefficient water vapor (m^2 kg^-1)
real    :: kappa_d         = 0.0  !absorption coefficient dry air (m^2 kg^-1)


real, allocatable, dimension(:) :: ss, cos_lat, solar_tau_0, p2
real, allocatable, dimension(:,:) :: solar, lw_tau_0
real, allocatable, dimension(:,:)              :: b_surf
real, allocatable, dimension(:,:,:)   :: b, tdt_rad, entrop_rad, tdt_sw
real, allocatable, dimension(:,:,:) :: up, down, net, solar_down, flux_rad, flux_sw
real, allocatable, dimension(:,:,:) :: dtrans, lw_tau
real, allocatable, dimension(:,:)   :: solar_tau
real, allocatable, dimension(:,:)   :: olr, swin
real, allocatable, dimension(:,:)   :: albedo, column_water_vapor
logical, allocatable, dimension(:,:)   :: ocean_mask
real, allocatable, dimension(:,:,:) :: dp_half
real, allocatable, dimension(:)   :: blon, blat

real, save :: pi, deg_to_rad , rad_to_deg

namelist/radiation_nml/ os08, constant_albedo, solar_constant, del_sol, &
           ir_tau_eq, ir_tau_pole, atm_abs, sw_diff, &
           linear_tau, wv_tau, do_column_wv, del_sw, albedo_value, &
           window, wv_exponent, solar_exponent, reference_slp, &
           tidally_locked, noon_longitude, &
           albedo_eq, albedo_pole, albedo_land, &
           moist_tau, kappa_v, kappa_d


!==================================================================================
!-------------------- diagnostics fields -------------------------------

integer :: id_olr, id_swdn_sfc, id_swdn_toa, id_lwdn_sfc, id_lwup_sfc, &
           id_tdt_rad, id_flux_rad, id_flux_lw, id_flux_sw, id_entrop_rad, id_tdt_sw

character(len=10), parameter :: mod_name = 'two_stream'

real :: missing_value = -999.


contains



! ==================================================================================
! ==================================================================================


subroutine radiation_init(is, ie, js, je, num_levels, axes, Time, lat)

!-------------------------------------------------------------------------------------
integer, intent(in), dimension(4) :: axes
type(time_type), intent(in)       :: Time
integer, intent(in)               :: is, ie, js, je, num_levels
real, intent(in) , dimension(:,:)   :: lat
!-------------------------------------------------------------------------------------
integer, dimension(3) :: half = (/1,2,4/)
integer :: ierr, io, unit, i

logical :: water_file_exists

!-----------------------------------------------------------------------------------------
! read namelist and copy to logfile

unit = open_file ('input.nml', action='read')
ierr=1
do while (ierr /= 0)
   read  (unit, nml=radiation_nml, iostat=io, end=10)
   ierr = check_nml_error (io, 'radiation_nml')
enddo
10 call close_file (unit)

unit = open_file ('logfile.out', action='append')
if ( mpp_pe() == 0 ) then
  write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
  write (unit, nml=radiation_nml)
endif
call close_file (unit)

pi    = 4.0*atan(1.)
deg_to_rad = 2.*pi/360.
rad_to_deg = 360.0/2./pi

initialized = .true.

allocate (b                (ie-is+1, je-js+1, num_levels))
allocate (tdt_rad          (ie-is+1, je-js+1, num_levels))
allocate (tdt_sw           (ie-is+1, je-js+1, num_levels))
allocate (entrop_rad       (ie-is+1, je-js+1, num_levels))

allocate (dtrans           (ie-is+1, je-js+1, num_levels))
allocate (dp_half          (ie-is+1, je-js+1, num_levels))

allocate (up               (ie-is+1, je-js+1, num_levels+1))
allocate (down             (ie-is+1, je-js+1, num_levels+1))
allocate (net              (ie-is+1, je-js+1, num_levels+1))
allocate (solar_down       (ie-is+1, je-js+1, num_levels+1))
allocate (flux_rad         (ie-is+1, je-js+1, num_levels+1))
allocate (flux_sw          (ie-is+1, je-js+1, num_levels+1))

allocate (lw_tau              (ie-is+1, je-js+1, num_levels+1))

allocate (b_surf           (ie-is+1, je-js+1))
allocate (olr              (ie-is+1, je-js+1))
allocate (swin             (ie-is+1, je-js+1))
allocate (albedo           (ie-is+1, je-js+1))
allocate (ocean_mask       (ie-is+1, je-js+1))

allocate (lw_tau_0            (ie-is+1, je-js+1))
allocate (column_water_vapor(ie-is+1, je-js+1))

allocate (solar_tau        (je-js+1, num_levels+1))

allocate (ss               (je-js+1))
allocate (cos_lat          (je-js+1))
allocate (solar            (ie-is+1, je-js+1))
allocate (solar_tau_0      (je-js+1))
allocate (p2               (je-js+1))

! variables for getting ocean_mask from topography mod -Tim M.
allocate(blon(is:ie+1), blat(js:je+1))
!call get_grid_boundaries(blon, blat)
! -Tim M.

! set albedo
if ( constant_albedo ) then
   albedo(:,:) = albedo_value
else
   do i = is, ie
      albedo(i,:) = albedo_pole + (albedo_eq - albedo_pole)*cos(lat(1,:))
   enddo
   water_file_exists = get_ocean_mask(blon, blat, ocean_mask)
   where ( ocean_mask == .false. )
      albedo = albedo + albedo_land*cos(lat(:,:))
   endwhere
endif


!-----------------------------------------------------------------------
!------------ initialize diagnostic fields ---------------

    id_olr = &
    register_diag_field ( mod_name, 'olr', axes(1:2), Time, &
               'outgoing longwave radiation', &
               'watts/m2', missing_value=missing_value               )
    id_swdn_sfc = &
    register_diag_field ( mod_name, 'swdn_sfc', axes(1:2), Time, &
               'SW flux down at surface', &
               'watts/m2', missing_value=missing_value               )
    id_swdn_toa = &
    register_diag_field ( mod_name, 'swdn_toa', axes(1:2), Time, &
               'SW flux down at TOA', &
               'watts/m2', missing_value=missing_value               )
    id_lwup_sfc = &
    register_diag_field ( mod_name, 'lwup_sfc', axes(1:2), Time, &
               'LW flux up at surface', &
               'watts/m2', missing_value=missing_value               )

    id_lwdn_sfc = &
    register_diag_field ( mod_name, 'lwdn_sfc', axes(1:2), Time, &
               'LW flux down at surface', &
               'watts/m2', missing_value=missing_value               )

    id_tdt_rad = &
        register_diag_field ( mod_name, 'tdt_rad', axes(1:3), Time, &
               'Temperature tendency due to radiation', &
               'K/s', missing_value=missing_value               )

    id_tdt_sw  = &
        register_diag_field ( mod_name, 'tdt_sw', axes(1:3), Time, &
               'Temperature tendency due to SW radiation', &
               'K/s', missing_value=missing_value               )

    id_flux_rad = &
        register_diag_field ( mod_name, 'flux_rad', axes(half), Time, &
               'Total radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_lw = &
        register_diag_field ( mod_name, 'flux_lw', axes(half), Time, &
               'Net longwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_sw = &
        register_diag_field ( mod_name, 'flux_sw', axes(half), Time, &
               'Net shortwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_entrop_rad = &
            register_diag_field ( mod_name, 'entrop_rad', axes(1:3), Time, &
               'Entropy production by radiation', &
               '1/s', missing_value=missing_value               )


return
end subroutine radiation_init


! ==================================================================================

subroutine radiation_down (is, js, Time_diag, lat, lon, p_half, q, t,         &
                           net_surf_sw_down, surf_lw_down)

! Begin the radiation calculation by computing downward fluxes.
! This part of the calculation does not depend on the surface temperature.

integer, intent(in)                 :: is, js
type(time_type), intent(in)         :: Time_diag
real, intent(in) , dimension(:,:)   :: lat, lon
real, intent(out) , dimension(:,:)   :: net_surf_sw_down
real, intent(out) , dimension(:,:)   :: surf_lw_down
real, intent(in) , dimension(:,:,:) :: t, p_half, q


integer :: i, j, k, n

logical :: used


n = size(t,3)

ss  = sin(lat(1,:))
p2 = (1. - 3.*ss*ss)/4.

cos_lat  = cos(lat(1,:))

if (tidally_locked) then
   do i = 1, size(t,1)
      do j = 1, size(t,2)
         solar(i,j) = solar_constant*cos(lat(i,j))*cos(lon(i,j) - noon_longitude*deg_to_rad)
         if ( solar(i,j) .LT. 0.0 ) then
            solar(i,j) = 0.0
         endif
      enddo
   enddo
else
   do i = 1, size(t,1)
      solar(i,:) = 0.25*solar_constant*(1.0 + del_sol*p2 + del_sw * ss)
!      solar(i,:) = solar_constant*cos_lat/pi
   enddo
endif

if ( os08 ) then
   do i = 1, size(t,1)
      lw_tau_0(i,:) = ir_tau_eq + (ir_tau_pole - ir_tau_eq)*ss*ss
   enddo
elseif ( moist_tau ) then
   ! (dsa, 3/9/10)
    do k=1, n
       dp_half(:,:,k) = p_half(:,:,k+1) - p_half(:,:,k)
    end do
else
   ! Tim M. mods
   if (do_column_wv) then
      ! use column water vapor
      do k=1, n
         dp_half(:,:,k) = p_half(:,:,k+1) - p_half(:,:,k)
      end do
      
      column_water_vapor = 0.0
      do i=1, size(t,1)
         do j=1, size(t,2)
            do k=1, n
               column_water_vapor(i,j) = column_water_vapor(i,j) + &
                    q(i,j,k) * dp_half(i,j,k)
                    !q(i,j,k)/(1.0 - q(i,j,k)) * dp_half(i,j,k)
            enddo
            column_water_vapor(i,j) = column_water_vapor(i,j)/grav
            lw_tau_0(i,j) = wv_tau * column_water_vapor(i,j)/10.0
         enddo
      enddo
   endif

endif

solar_tau_0 = (1.0 - sw_diff*ss*ss)*atm_abs

!longwave optical depth
lw_tau = 0.
do k = 1, n

   lw_tau(:,:,k+1) = lw_tau(:,:,k) + q(:,:,k) * dp_half(:,:,k)
   solar_tau(:,k) = solar_tau_0(:)*(p_half(1,:,k)/reference_slp)**solar_exponent
end do
lw_tau = lw_tau * wv_tau / grav / 10.
do i = 1, size(t,1)
   do j = 1, size(t,2)
      lw_tau(i,j,:) = lw_tau(i,j,:) + &
                      linear_tau * p_half(i,j,:) / p_half(i,j,n+1)
   end do
end do

! no radiation from spectral window
b = (1.0-window)*stefan*t*t*t*t

do k = 1, n
   do i = 1, size(t,1)
      dtrans(i,:,k) = exp(-(lw_tau(i,:,k+1)-lw_tau(i,:,k)))
   enddo
end do

down(:,:,1) = 0.0
do k = 1,n
  do j =1, size(t,2)
    down(:,j,k+1) = down(:,j,k)*dtrans(:,j,k) + b(:,j,k)*(1.0 - dtrans(:,j,k))
  end do
end do

do i = 1, size(t,1)
   do j = 1, size(t,2)
      do k = 1,n+1
         solar_down(i,j,k) = solar(i,j)*exp(-solar_tau(j,k))
      end do
   end do
end do

surf_lw_down     = down(:,:,n+1)
net_surf_sw_down = solar_down(:,:,n+1)*(1. - albedo(:,:))
swin = solar_down(:,:,1)
    

!------- downward sw flux surface -------
      if ( id_swdn_sfc > 0 ) then
          used = send_data ( id_swdn_sfc, net_surf_sw_down, Time_diag)
      endif
!------- incoming sw flux toa -------
      if ( id_swdn_toa > 0 ) then
          used = send_data ( id_swdn_toa, swin, Time_diag)
      endif
!------- downward lw flux surface -------
      if ( id_lwdn_sfc > 0 ) then
          used = send_data ( id_lwdn_sfc, surf_lw_down, Time_diag)
      endif

return
end subroutine radiation_down

! ==================================================================================

subroutine radiation_up (is, js, Time_diag, lat, p_half, t_surf, t, tdt)

! Now complete the radiation calculation by computing the upward and net fluxes.

integer, intent(in)                 :: is, js
type(time_type), intent(in)         :: Time_diag
real, intent(in) , dimension(:,:)   :: lat
real, intent(in) , dimension(:,:)   :: t_surf
real, intent(in) , dimension(:,:,:) :: t, p_half
real, intent(inout), dimension(:,:,:) :: tdt


integer :: i, j, k, n

logical :: used

n = size(t,3)

! total flux from surface
b_surf = stefan*t_surf*t_surf*t_surf*t_surf

! first deal with non-window upward flux
up(:,:,n+1) = b_surf*(1.0-window)
do k = n,1,-1
  do j = 1, size(t,2)
    up(:,j,k) = up(:,j,k+1)*dtrans(:,j,k) + b(:,j,k)*(1.0 - dtrans(:,j,k))
  end do
end do

! add upward flux in spectral window
do k = 1,n+1
 up(:,:,k) = up(:,:,k) + b_surf(:,:)*window
end do

do k = 1,n+1
  net(:,:,k) = up(:,:,k)-down(:,:,k)
  flux_sw(:,:,k) = albedo(:,:)*solar_down(:,:,n+1) - solar_down(:,:,k)
  flux_rad(:,:,k) = net(:,:,k) + flux_sw(:,:,k)
end do

do k = 1,n
  tdt_rad(:,:,k) = (net(:,:,k+1) - net(:,:,k) - solar_down(:,:,k+1) + solar_down(:,:,k))  
             !*grav/(cp_air*(p_half(:,:,k+1)-p_half(:,:,k)))
  tdt_sw(:,:,k) = (- solar_down(:,:,k+1) + solar_down(:,:,k))  &
             *grav/(cp_air*(p_half(:,:,k+1)-p_half(:,:,k)))
  tdt(:,:,k) = tdt(:,:,k) + tdt_rad(:,:,k)
end do


olr = up(:,:,1)
    

!------- outgoing lw flux toa (olr) -------
      if ( id_olr > 0 ) then
          used = send_data ( id_olr, olr, Time_diag)
      endif
!------- upward lw flux surface -------
      if ( id_lwup_sfc > 0 ) then
          used = send_data ( id_lwup_sfc, b_surf, Time_diag)
      endif
!------- temperature tendency due to radiation ------------
      if ( id_tdt_rad > 0 ) then
         used = send_data ( id_tdt_rad, tdt_rad, Time_diag)
      endif
      if ( id_tdt_sw > 0 ) then
         used = send_data ( id_tdt_sw, tdt_sw, Time_diag)
      endif
!------- total radiative flux (at half levels) -----------
      if ( id_flux_rad > 0 ) then
         used = send_data ( id_flux_rad, flux_rad, Time_diag)
      endif
!------- longwave radiative flux (at half levels) --------
      if ( id_flux_lw > 0 ) then 
         used = send_data ( id_flux_lw, net, Time_diag)
      endif
      if ( id_flux_sw > 0 ) then
         used = send_data ( id_flux_sw, flux_sw, Time_diag)
      endif
      if ( id_entrop_rad > 0 ) then
         do k=1,n 
            entrop_rad(:,:,k) =tdt_rad(:,:,k)/t(:,:,k)*p_half(:,:,n+1)/1.e5
         end do
         used = send_data ( id_entrop_rad, entrop_rad, Time_diag)
      endif

return
end subroutine radiation_up

! ==================================================================================

                                                                      
subroutine radiation_end
                                                                                                      
deallocate (b, tdt_rad, tdt_sw, entrop_rad) 
deallocate (up, down, net, solar_down, flux_rad, flux_sw)
deallocate (b_surf, olr, swin, albedo)
deallocate (dtrans)   
deallocate (lw_tau, solar_tau)
deallocate (ss, solar, lw_tau_0, solar_tau_0, p2) 

end subroutine radiation_end

! ==================================================================================

end module radiation_mod




