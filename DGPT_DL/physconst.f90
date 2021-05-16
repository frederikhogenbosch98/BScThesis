module physconst

!---------------------------------------------------------------------
! This module contains some physical constants
!---------------------------------------------------------------------

  use f90_kind
  implicit none

  real(dp), parameter :: pi       = 3.141592653589793_dp
  real(dp), parameter :: twopi    = 2.d0*pi
  real(dp), parameter :: J_per_eV = 1.602176620898E-19_dp
  real(dp), parameter :: m_p      = 938.2720813_dp ! Proton mass in  MeV / c^2
  real(dp), parameter :: hbc      = 0.19732697E-10_dp ! Hbar * c in  MeV * cm
  real(dp), parameter :: m_e      = 0.5109989461_dp ! Electron mass in MeV / c^2
  real(dp), parameter :: alpha    = 0.0072973525664_dp ! Fine structure constant dimensionless
  real(dp), parameter :: eps0     = 1.418284572502546E-26_dp ! vacuum permitivity in C^2 Mev^-1 cm^-1
  real(dp), parameter :: ee       = 1.6021766208E-19_dp ! elementary charge in C
  real(dp), parameter :: M_water  = 18.01528 ! Water molar mass in g / mol
  real(dp), parameter :: N_A      = 6.022140857E23_dp ! Avogadro's number in particles/mol
  real(dp), parameter :: rho_H    = 0.07 ! g/cm^3
  real(dp), parameter :: rho_O    = 1.141 ! g/cm^3
  real(dp), parameter :: M_H      = 1.008 ! g/mol
  real(dp), parameter :: M_O      = 15.999 ! g/mol
  integer, parameter, dimension(2) :: Z_t = (/ 1,8 /) ! atomic numbers of H and O
  real(dp),parameter, dimension(2) :: IPot = (/ 19.0_dp, 11.2_dp+11.7_dp * 8 /) * 1.0E-6_dp ! Ionization potentials of H and O in MeV
  
  
end module physconst
