
module constants_mod

! <CONTACT EMAIL="Bruce.Wyman@noaa.gov">
!   Bruce Wyman
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!    Defines useful constants for Earth.
! </OVERVIEW>

! <DESCRIPTION>
!   Constants are defined as real parameters.
!   Constants are accessed through the "use" statement.
! </DESCRIPTION>

implicit none
private

character(len=128) :: version='$Id: constants.f90,v 13.0 2006/03/28 21:37:37 fms Exp $'
character(len=128) :: tagname='$Name: latest $'
!dummy variable to use in HUGE initializations
real :: realnumber

!------------ physical constants ---------------
! <DATA NAME="RADIUS" UNITS="m" TYPE="real" DEFAULT="6371.e3">
!   radius of the earth
! </DATA>
! <DATA NAME="OMEGA" UNITS="1/s" TYPE="real" DEFAULT="7.292e-5">
!   rotation rate of the planet (earth)
! </DATA>
! <DATA NAME="GRAV" UNITS="m/s^2" TYPE="real" DEFAULT="9.80">
!   acceleration due to gravity
! </DATA>
! <DATA NAME="RDGAS" UNITS="J/kg/deg" TYPE="real" DEFAULT="287.04">
!   gas constant for dry air
! </DATA>
! <DATA NAME="KAPPA" TYPE="real" DEFAULT="2./7.">
!   RDGAS / CP_AIR
! </DATA>
! <DATA NAME="CP_AIR" UNITS="J/kg/deg" TYPE="real" DEFAULT="RDGAS/KAPPA">
!   specific heat capacity of dry air at constant pressure
! </DATA>
! <DATA NAME="CP_OCEAN" UNITS="J/kg/deg" TYPE="real" DEFAULT="3989.24495292815">
!   specific heat capacity taken from McDougall (2002) "Potential Enthalpy ..."
! </DATA>
! <DATA NAME="RHO0" UNITS="kg/m^3" TYPE="real" DEFAULT="1.035e3">
!   average density of sea water
! </DATA>
! <DATA NAME="RHO0R" UNITS="m^3/kg" TYPE="real" DEFAULT="1.0/RHO0">
!   reciprocal of average density of sea water
! </DATA>
! <DATA NAME="RHO_CP" UNITS="J/m^3/deg" TYPE="real" DEFAULT="RHO0*CP_OCEAN">
!   (kg/m^3)*(cal/kg/deg C)(joules/cal) = (joules/m^3/deg C)
! </DATA>

real, public, parameter :: RADIUS = 6371.0e3   
real, public, parameter :: OMEGA  = 7.292e-5 !/50. !/365.0
real, public, parameter :: GRAV   = 9.80    
real, public, parameter :: RDGAS  = 287.00421 
real, public, parameter :: KAPPA  = 0.2858607676 
real, public, parameter :: CP_AIR = 1004.0 !RDGAS/KAPPA 
real, public, parameter :: CP_OCEAN = 3989.24495292815
real, public, parameter :: RHO0    = 1.035e3
real, public, parameter :: RHO0R   = 1.0/RHO0
real, public, parameter :: RHO_CP  = RHO0*4187.  !CP_OCEAN

!------------ water vapor constants ---------------
! <DATA NAME="RVGAS" UNITS="J/kg/deg" TYPE="real" DEFAULT="461.50">
!   gas constant for water vapor
! </DATA>
! <DATA NAME="CP_VAPOR" UNITS="J/kg/deg" TYPE="real" DEFAULT="4.0*RVGAS">
!   specific heat capacity of water vapor at constant pressure
! </DATA>
! <DATA NAME="DENS_H2O" UNITS="kg/m^3" TYPE="real" DEFAULT="1000.">
!   density of liquid water
! </DATA>
! <DATA NAME="HLV" UNITS="J/kg" TYPE="real" DEFAULT="2.500e6">
!   latent heat of evaporation
! </DATA>
! <DATA NAME="HLF" UNITS="J/kg" TYPE="real" DEFAULT="3.34e5">
!   latent heat of fusion
! </DATA>
! <DATA NAME="HLS" UNITS="J/kg" TYPE="real" DEFAULT="2.834e6">
!   latent heat of sublimation
! </DATA>
! <DATA NAME="TFREEZE" UNITS="degK" TYPE="real" DEFAULT="273.16">
!   temp where fresh water freezes
! </DATA>

real, public, parameter :: RVGAS = 461.91733246 !461.50 
real, public, parameter :: CP_VAPOR = 1847.0 !4.0*RVGAS
real, public, parameter :: DENS_H2O = 1000. 
real, public, parameter :: HLV = 2.493e6 !2.500e6   
real, public, parameter :: HLF = 3.34e5   
real, public, parameter :: HLS = HLV + HLF
real, public, parameter :: TFREEZE = 273.16    

!-------------- radiation constants -----------------

! <DATA NAME="WTMAIR" UNITS="AMU" TYPE="real" DEFAULT="2.896440E+01">
!  molecular weight of air 
! </DATA>
! <DATA NAME="WTMH2O" UNITS="AMU" TYPE="real" DEFAULT="1.801534E+01">
!  molecular weight of water
! </DATA>
! <DATA NAME="WTMO3" UNITS="AMU" TYPE="real" DEFAULT="47.99820E+01">
!   molecular weight of ozone
! </DATA>
! <DATA NAME="DIFFAC" TYPE="real" DEFAULT="1.660000E+00">
! diffusivity factor
! </DATA>
! <DATA NAME="SECONDS_PER_DAY" UNITS="seconds" TYPE="real" DEFAULT="8.640000E+04">
! seconds in a day
! </DATA>
! <DATA NAME="AVOGNO" UNITS="atoms/mole" TYPE="real" DEFAULT="6.023000E+23">
!  Avogadro's number 
! </DATA>
! <DATA NAME="PSTD" UNITS="dynes/cm^2" TYPE="real" DEFAULT="1.013250E+06">
!  mean sea level pressure
! </DATA>
! <DATA NAME="PSTD_MKS" UNITS="Newtons/m^2" TYPE="real" DEFAULT="101325.0">
!  mean sea level pressure
! </DATA>

real, public, parameter :: WTMAIR = 2.896440E+01
real, public, parameter :: WTMH2O = WTMAIR*(RDGAS/RVGAS) !pjp OK to change value because not used yet.
real, public, parameter :: WTMO3  = 47.99820E+01
real, public, parameter :: DIFFAC = 1.660000E+00
real, public, parameter :: SECONDS_PER_DAY  = 8.640000E+04
real, public, parameter :: AVOGNO = 6.023000E+23
real, public, parameter :: PSTD   = 1.013250E+06
real, public, parameter :: PSTD_MKS    = 101325.0
!real, public, parameter :: REARTH  = 6.356766E+08 !pjp Not used anywhere. 

! <DATA NAME="RADCON" UNITS="deg sec/(cm day)" TYPE="real" DEFAULT="((1.0E+02*GRAV)/(1.0E+04*CP_AIR))*SECONDS_PER_DAY">
!  factor used to convert flux divergence to heating rate in degrees per day
! </DATA>
! <DATA NAME="RADCON_MKS" UNITS="deg sec/(m day)" TYPE="real" DEFAULT="(GRAV/CP_AIR)*SECONDS_PER_DAY">
!  factor used to convert flux divergence to heating rate in degrees per day
! </DATA>
! <DATA NAME="O2MIXRAT" TYPE="real" DEFAULT="2.0953E-01">
! mixing ratio of molecular oxygen in air
! </DATA>
! <DATA NAME="RHOAIR" UNITS="kg/m^3" TYPE="real" DEFAULT="1.292269">
!  reference atmospheric density
! </DATA>
! <DATA NAME="ALOGMIN" TYPE="real" DEFAULT="-50.0">
!  minimum value allowed as argument to log function
! </DATA>

real, public, parameter :: RADCON = ((1.0E+02*GRAV)/(1.0E+04*CP_AIR))*SECONDS_PER_DAY
real, public, parameter :: RADCON_MKS  = (GRAV/CP_AIR)*SECONDS_PER_DAY
real, public, parameter :: O2MIXRAT    = 2.0953E-01
real, public, parameter :: RHOAIR      = 1.292269
real, public, parameter :: ALOGMIN     = -50.0

!------------ miscellaneous constants ---------------
! <DATA NAME="STEFAN" UNITS="W/m^2/deg^4" TYPE="real" DEFAULT="5.6734e-8">
!   Stefan-Boltzmann constant
! </DATA>
! <DATA NAME="VONKARM"  TYPE="real" DEFAULT="0.40">
!   Von Karman constant
! </DATA>
! <DATA NAME="PI" TYPE="real" DEFAULT="3.14159265358979323846">
!    ratio of circle circumference to diameter
! </DATA>
! <DATA NAME="RAD_TO_DEG"  TYPE="real" DEFAULT="180.0/PI">
!   degrees per radian
! </DATA>
! <DATA NAME="DEG_TO_RAD"  TYPE="real" DEFAULT="PI/180.0">
!   radians per degree
! </DATA>
! <DATA NAME="RADIAN"  TYPE="real" DEFAULT="180.0/PI">
!   equal to RAD_TO_DEG. Named RADIAN for backward compatability.
! </DATA>
! <DATA NAME="C2DBARS" UNITS="dbars" TYPE="real" DEFAULT="1.e-4">
!   converts rho*g*z (in mks) to dbars: 1dbar = 10^4 (kg/m^3)(m/s^2)m
! </DATA>
! <DATA NAME="KELVIN" TYPE="real" DEFAULT="273.15">
!   degrees Kelvin at zero Celsius
! </DATA>
! <DATA NAME="EPSLN" TYPE="real" DEFAULT="1.0e-40">
!   a small number to prevent divide by zero exceptions
! </DATA>

real, public, parameter :: STEFAN  = 5.6734e-8 
real, public, parameter :: VONKARM = 0.40     
real, public, parameter :: PI      = 3.14159265358979323846
real, public, parameter :: RAD_TO_DEG=180./PI
real, public, parameter :: DEG_TO_RAD=PI/180.
real, public, parameter :: RADIAN  = RAD_TO_DEG
real, public, parameter :: C2DBARS = 1.e-4
real, public, parameter :: KELVIN  = 273.15
real, public, parameter :: EPSLN   = 1.0e-40

real, public, parameter :: cp_water = 4187.
real, public, parameter :: r_vir = RVGAS/RDGAS-1.
real, public, parameter :: cp_vir = CP_VAPOR/CP_AIR-1.

  real, public, parameter :: c      = 2.99792458d8    ! speed of light [m/s]
  real, public, parameter :: h      = 6.6260693d-34   ! Planck constant [J s]
  real, public, parameter :: kB     = 1.3806488d-23   ! Boltzmann constant [m2/kg/s2/K]
  real, public, parameter :: c2     = (h*c/kB)*1.0d2  ! a constant [cm/K]
  real, public, parameter :: atm    = 101325.0d0      ! Pa to atm conversion
  real, public, parameter :: Rstar  = 8.3144621       ! ideal gas constant [J/mol/K]
  real, public, parameter :: Tref   = 296.d0          ! hitran standard temperature [K]
  real, public, parameter :: mpr    = 1.672621e-27    ! mass of proton [kg]

  complex, public, parameter :: i   = dcmplx(0.0d0,+1.0d0)
 
!-----------------------------------------------------------------------
! version and tagname published
! so that write_version_number can be called for constants_mod by fms_init
public :: version, tagname
!-----------------------------------------------------------------------
public :: constants_init, latent_heat, B_nu

contains

subroutine constants_init

! dummy routine.

end subroutine constants_init


subroutine latent_heat(T, L)

real, intent(in) :: T
real, intent(out) :: L

L = hlv + (cp_vapor - cp_water)*(T - 273.15) ! phys.water.TriplePointT)

end subroutine latent_heat


  subroutine B_nu(nu_cm,T,Bnu)

    real, intent(in) :: nu_cm ! wavenumber [cm^-1]
    real, intent(in) :: T     ! temperature [K]
    real, intent(out) :: Bnu  ! Planck spectral irradiance [W/m2/cm^-1/sr]
    !real(8), parameter :: c1 = 2.0d8*h*c**2 ! W/m2/cm^-4/sr
    real :: nu, Btemp

    nu = 1.0d2*c*nu_cm ! cm^-1 --> Hz

    !Bnu = c1*nu_cm**3/(exp(c2*nu_cm/T) - 1.0d0)
    Btemp = (2*h*nu**3/c**2)/(exp(h*nu/(kB*T))  -  1.0d0) ! W/m2/sr/Hz
    Bnu   = 1.0d2*c*Btemp  ! W/m2/sr/cm^-1
    
  end subroutine B_nu


end module constants_mod

! <INFO>

!   <FUTURE>               
!   1.  Renaming of constants.
!   </FUTURE>               
!   <FUTURE>               
!   2.  Additional constants.
!   </FUTURE>
!   <NOTE>
!    Constants have been declared as type REAL, PARAMETER.
!
!    The value a constant can not be changed in a users program.
!    New constants can be defined in terms of values from the
!    constants module using a parameter statement.<br><br>
!
!    The name given to a particular constant may be changed.<br><br>
!
!    Constants can be used on the right side on an assignment statement
!    (their value can not be reassigned). 
!
!
!<TESTPROGRAM NAME="EXAMPLE">
!<PRE>
!    use constants_mod, only:  TFREEZE, grav_new => GRAV
!    real, parameter :: grav_inv = 1.0 / grav_new
!    tempc(:,:,:) = tempk(:,:,:) - TFREEZE
!    geopotential(:,:) = height(:,:) * grav_new
!</PRE>
!</TESTPROGRAM>
!   </NOTE>

! </INFO>

module dimensions
  
  implicit none

  integer, parameter :: nGas      = 2           ! number of species
  integer, parameter :: nAng      = 1 !4            ! number of angles to integrate over
  integer, parameter :: nlinesMAX = 3000000      ! maximum number of lines to read
  integer, parameter :: nS        = 2000 !2000        ! spectral resolution (same in sw and lw for now)
  integer, parameter :: nTem      = 5            ! number of cross-section temperatures at each layer (1 or 2)
  
  ! directory containing data needed by the model
  character(len=50), save :: datadir='~/lbldata/'
  !character(len=40), save :: datadir='/Users/robin/programs/PCM_LBL/data/'
  !character(len=40), save :: datadir='/Volumes/robin/programs/PCM_LBL/data/'
  
end module dimensions
