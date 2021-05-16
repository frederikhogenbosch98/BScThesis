module heterogeneous
! This module calculates the modified stopping power S* = S(E) + 1/2 dT(E)/dE
! and the straggling coefficient T(E) for different
! materials based on eq 5.33 from Turner [1] and 4.16 from Kelsey [2]
! References:
! [1] -
! [2] - 

use f90_kind
use physconst
implicit none

contains

subroutine composition(material,I,ln_I,n,Z,N_i)
! Based on the material computes several quantities of interest
! Parameters:
! _____________________________________________________________________________
! Input/Output:
! 
! material       - "S"kin,"B"one,"W"ater,"T"umor
! I (eV)         - atomic excitation energies 
! ln_I           - log(material excitation energy)
! n (e/m^3)      - material electronic density
! Z              - atomic numbers array
!
! Other:
! 
! M(g/mol)       - Molar mass = standard atomic weight x molar mass const. (1 g/mol)
! N_i(atoms/m^3) - atomic densities of elements in the material
! percentages    - % fraction by weight of element in material 
! rho(g/cm^3)    - density of material
! _____________________________________________________________________________
! The convention for the Z,M,I arrays is: 
! _____________________________________________________________________________
! index    : 1    ,2     ,3     ,4     ,5     ,6     ,7     ,8     ,9     ,10
! element  : H    ,C     ,N     ,O     ,Na    ,P     ,S     ,Ca    ,Cl    ,K
! Z        : 1    ,6     ,7     ,8     ,11    ,15    ,16    ,20    ,17    ,19 
! M(g/mol) : 1.008,12.011,14.007,15.999,22.989,30.973,32.060,40.078,35.450,39.098
! _____________________________________________________________________________
    
character(len=1),intent(in)        :: material
real(dp),intent(out)               :: ln_I,n
real(dp),dimension(10),intent(out) :: I,N_i
integer,dimension(10),intent(out)  :: Z

real(dp),dimension(10):: M,ln_Ii,percentages
real(dp)              :: rho

Z = (/ 1, 6, 7, 8, 11, 15, 16, 20, 17, 19 /)
M = (/1.008_dp,12.011_dp,14.007_dp,15.999_dp,22.989_dp,30.973_dp,32.060_dp,40.078_dp,35.450_dp,39.098_dp /)
    
I(1)       = 19.0_dp
ln_Ii(1)   = log(I(1))
I(2:5)     = 11.2_dp+11.7_dp*Z(2:5)
ln_Ii(2:5) = log(I(2:5))
I(6:10)    = 52.8_dp+8.71_dp*Z(6:10)
ln_Ii(6:10)= log(I(6:10))

if (material=="S") then 
  percentages = (/10.0_dp,25.0_dp,4.6_dp,59.4_dp,0.2_dp,0.1_dp,0.3_dp,0.0_dp,0.3_dp,0.1_dp/)/100.0_dp
  rho = 1.09_dp
elseif (material=="B") then
  percentages = (/3.4_dp,11.5_dp,4.2_dp,43.5_dp,0.1_dp,0.0_dp,0.3_dp,22.5_dp,0.0_dp,0.0_dp/)/100.0_dp
  rho = 1.85_dp
elseif (material=="W") then
  percentages = (/11.19_dp,0.0_dp,0.0_dp,88.81_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/)/100.0_dp
  rho = 1.0_dp
elseif (material=="T") then
  percentages = (/0.6_dp,19.4_dp,2.5_dp,66.1_dp,0.2_dp,0.3_dp,0.2_dp,0.0_dp,0.3_dp,0.3_dp/)/100.0_dp
  rho = 1.2_dp
endif

N_i = 1.0E6_dp*rho*percentages*N_A/M
n = sum(N_i*Z)
ln_I = sum(N_i*Z*ln_Ii)/n

end subroutine composition

real(dp) function S(E,material)
! Gives the *modified* stopping power (Mev/cm) based on energy and
! material according to equation 5.33 from Turner (2.20 in Burlacu report)
! if onT then S = S + 1/2 dT/dE
! Parameters:
! ____________________________________________________________________
! E  (MeV) - Energy of proton
! material - "S"kin,"B"one,"W"ater,"T"umor
! 
! ____________________________________________________________________
    
character(len=1),intent(in):: material
real(dp),intent(in)        :: E
    
real(dp)                   :: ln_I, n, v_p, F, corr,dT_dE
logical                    :: onS,onT
integer,dimension(10)      :: Z
real(dp),dimension(10)     :: I,N_i

S = 0.0_dp

corr = 1.0_dp/(4.0_dp * pi * eps0)**2

call composition(material,I,ln_I,n,Z,N_i)

v_p = sqrt(2.0_dp*E/m_p) ! v_p is actualy Beta as m_p is in units of MeV/c^2
F = log( 1.02E6_dp * v_p**2 / (1.0_dp - v_p**2) ) - v_p**2

S = 5.08E-31_dp * n / v_p**2 * (F - ln_I)

I   = I   * 1.0E-6_dp
N_i = N_i * 1.0E-6_dp

dT_dE = sum( corr * N_i * &
    4.0_dp * pi * ee**4 * Z * &
    (-2.0_dp*I*m_p/(3.0_dp*m_e) * log(4.0_dp*m_e*E/(I*m_p))/E**2 + &
    2.0_dp*I*m_p/(3.0_dp*m_e*E)))

S = S + 1.0_dp/2.0_dp * dT_dE

end function S

real(dp) function T(E,material)
! Gives the straggling coefficient based on the energy
! and material type.
! T based on 2.23 in report Burlacu (Williams equation)
! Parameters:
! corr_fact [MeV^2cm^2/C^4] - correction factor used for T 

real(dp),intent(in)        :: E
character(len=1),intent(in):: material
    
real(dp)                   :: corr,ln_I,n,v_p
logical                    :: onT
integer,dimension(10)      :: Z
real(dp),dimension(10)     :: I,N_i

T = 0.0_dp

corr = 1.0_dp/(4.0_dp * pi * eps0)**2

call composition(material,I,ln_I,n,Z,N_i)

I   = I   * 1.0E-6_dp
N_i = N_i * 1.0E-6_dp

v_p = sqrt(2.0_dp*E/m_p)
T =  sum( corr * N_i * 4.0_dp * pi * ee**4 * Z * ( 1 + 4.0_dp*I/(3.0_dp*m_e*v_p**2) * log(2.0_dp*m_e*v_p**2/I) ) )

end function T
  
end module heterogeneous
