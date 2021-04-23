program test
  ! Description: This program solves d/dx phi(x,E) = d/dE S*phi(x,E) + 1/2 (d/dE)^2 T phi(x,E)
  !              For the T part SIPG is used. For the S part DG is used.
  !
  ! Parameters : no_grps (int) - Number of groups in which the E domain is divided
  !              no_steps(int) - Number of elements in which the x domain is divided
  !              dx       (cm) - Width of one element in the x domain

use f90_kind
use quadrature
use physconst
use heterogeneous
use fem
use CN
  
implicit none

real(dp)            :: x_max=8.5_dp
integer, parameter  :: no_grps=500
integer,parameter   :: no_steps=2000

integer,parameter   :: kl_G=3
integer,parameter   :: ku_G=3
real(dp),allocatable,dimension(:) :: point,weight
real(dp), dimension(:,:), allocatable :: M_band
real(dp), dimension(:,:), allocatable :: G_band

integer, parameter :: no_dof=2*no_grps
real(dp), dimension(no_dof) :: phi,phi_old
real(dp) :: dE(no_grps)
real(dp) :: E_bounds(no_grps+1)
real(dp) :: E_min,E_max,eval_point,Sum,k,A,mu,sigma,dx,E_low,E_high,phi_low,phi_high
integer  :: gr,start_row,end_row,start_col,end_col,row,col,step,pos,i,j,idx,kl_M,ku_M

! Set the energy domain and discretization (uniform)

E_max = 201.0_dp
E_min = 1.0_dp
dE = (E_max-E_min) / real(no_grps,dp)
E_bounds(1) = E_max
do gr=1,no_grps
  E_bounds(gr+1) = E_bounds(gr) - dE(gr)
enddo

!print *,E_bounds

dx = x_max / no_steps

! Init

phi     = 0.0_dp
phi_old = 0.0_dp

! Gaussian boundary condition

A     = 1.0_dp
mu    = 100.0_dp
sigma = 1.0_dp / sqrt(2.0_dp)

call project_gaussian_on_dg(A,mu,sigma,E_bounds,dE,phi)
phi_old = phi
print *,phi
! Plot flux

!do gr=no_grps,1,-1
!  E_low  = E_bounds(gr+1)
!  E_high = E_bounds(gr)
!  phi_low  = phi(2*(gr-1)+1) - phi(2*(gr-1)+2)
!  phi_high = phi(2*(gr-1)+1) + phi(2*(gr-1)+2)
!  print *,E_low,  phi_low
!  print *,E_high, phi_high
!enddo

! Calc mass matrix

call calc_M(M_band,kl_M,ku_M,dE)

! Stepping

do step=1,no_steps
  ! Construct G matrix with CSD and straggling

  call build_G_band(no_dof,E_bounds,dE,G_band,kl_G,ku_G)

  ! Do single CN step

  call CN_1step(no_dof,M_band,kl_M,ku_M,G_band,kl_G,ku_G,phi_old,phi,dx)

  ! Time copy

  phi_old = phi
enddo

! Plot flux

do gr=no_grps,1,-1
  E_low  = E_bounds(gr+1)
  E_high = E_bounds(gr)
  phi_low  = phi(2*(gr-1)+1) - phi(2*(gr-1)+2)
  phi_high = phi(2*(gr-1)+1) + phi(2*(gr-1)+2)
!  print *,E_low,  phi_low
!  print *,E_high, phi_high
enddo

end program test
