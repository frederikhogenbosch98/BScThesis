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
integer, parameter :: int_size=40
integer,parameter   :: kl_G=3
integer,parameter   :: ku_G=3
real(dp),allocatable,dimension(:) :: point,weight
real(dp), dimension(:,:), allocatable :: M_band
real(dp), dimension(:,:), allocatable :: G_band

integer, parameter :: no_dof=2*no_grps
real(dp), dimension(no_dof) :: phi,phi_old
real(dp) :: dE(no_grps)
real(dp) :: E_bounds(no_grps+1)
real(dp) :: E_bounded(int_size)
real(dp) :: E_bounded_old(int_size)
real(dp) :: phi_bounded(int_size*2)
real(dp) :: phi_bounded_old(int_size*2)
real(dp) :: E_min,E_max,eval_point,Sum,k,A,mu,sigma,dx,E_low,E_high,phi_low,phi_high
integer  :: gr,start_row,end_row,start_col,end_col,row,col,step,pos,i,j,idx,kl_M,ku_M,n_max

! timing variables
integer count_0, count_1
integer count_rate, count_max
double precision time_init, time_final, elapsed_time


! Set the energy domain and discretization (uniform)

E_max = 201.0_dp
E_min = 1.0_dp
dE = (E_max-E_min) / real(no_grps,dp)
E_bounds(1) = E_max
do gr=1,no_grps
  E_bounds(gr+1) = E_bounds(gr) - dE(gr)
enddo

!print *,E_bounds
!print *,dE
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
!print *,phi
! Plot flux
open(unit=12, file='data.txt')
!do i=1,no_grps
!    print *,phi(i)
!    write(12,*) phi(i)
!enddo

!close(12)

!do gr=no_grps,1,-1
!  E_low  = E_bounds(gr+1)
!  E_high = E_bounds(gr)
!  phi_low  = phi(2*(gr-1)+1) - phi(2*(gr-1)+2)
!  phi_high = phi(2*(gr-1)+1) + phi(2*(gr-1)+2)
!  print *,E_low,  phi_low
!  print *,E_high, phi_high
!enddo
open(unit=14, file='dataphi.txt')
! Calc mass matrix

call calc_M(M_band,kl_M,ku_M,dE)
!print *,"size(M_band) = ", size(M_band)
!print *,"shape(M_band) = ",shape(M_band)
n_max = 230
call det_bounds(E_bounds, dE, phi_old, E_bounded, phi_bounded, n_max)

!write(14,*) phi_bounded
!write(12,*) E_bounded

!Timing
call system_clock(count_0, count_rate, count_max)
time_init=count_0*1.0/count_rate

!print *, size(phi)
!print *, size(E_bounds)
E_bounded_old = E_bounded
! Stepping
do step=1,2000
  ! Construct G matrix with CSD and straggling
  if (mod(step,200)==0) then    
  call update_bounds(E_bounds, phi_old, E_bounded_old, dE, phi_bounded_old, E_bounded, phi_bounded, step, n_max, int_size, no_steps)
  write(14,*) phi_bounded
  write(12,*) E_bounded
  !print *, step
  endif

  !call det_E_bounds(E_bounds, dE, E_max, E_min, step)

  call build_G_band(no_dof,E_bounds,dE,G_band,kl_G,ku_G,step)
  ! Do single CN step

  call CN_1step(no_dof,M_band,kl_M,ku_M,G_band,kl_G,ku_G,phi_old,phi,dx)

  ! Time copy
  phi_old = phi
  !print *, phi
enddo

close(12)
close(14)

call system_clock(count_1, count_rate, count_max)
time_final = count_1*1.0/count_rate
elapsed_time = time_final-time_init

print *,"elasped time: ",elapsed_time

!print *,"size of phi = ", size(phi)
!print *,"shape of phi = ", shape(phi)
!print *,"size(G_band) = ", shape(G_band)
! Plot flux
do gr=no_grps,1,-1
    E_low  = E_bounds(gr+1)
    E_high = E_bounds(gr)
    phi_low  = phi(2*(gr-1)+1) - phi(2*(gr-1)+2)
    phi_high = phi(2*(gr-1)+1) + phi(2*(gr-1)+2)
!    print *,E_low,  phi_low
!    print *,E_high, phi_high
    !write(14,*) phi_high
    !write(12,*) E_high
enddo

!close(14)
!close(12)


end program test
