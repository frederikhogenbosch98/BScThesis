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
integer, parameter :: int_size=75
integer, parameter :: n_max=241
!integer, parameter :: n_max=241-12
integer,parameter   :: kl_G=3
integer,parameter   :: ku_G=3
real(dp),allocatable,dimension(:) :: point,weight
real(dp), dimension(:,:), allocatable :: M_band
real(dp), dimension(:,:), allocatable :: G_band

integer, parameter :: no_dof=3*no_grps
real(dp), dimension(no_dof) :: phi,phi_old,phi_un
real(dp) :: dE(no_grps)
real(dp) :: dE_bounded(int_size)
real(dp) :: E_bounds(no_grps+1)
real(dp) :: E_bounded(int_size+1)
real(dp) :: E_bounded_old(int_size+1)
real(dp) :: phi_bounded(int_size*3)
real(dp) :: phi_bounded_old(int_size*3)
real(dp) :: phi_bounded_un(int_size*3)
real(dp) :: phi_non_CN(int_size*3)
real(dp), dimension(size(phi_old)) :: phi_proj
real(dp), dimension(int_size) :: E_lowb, E_highb, E_avg
real(dp), dimension(int_size) :: phi_lowb, phi_highb, phi_avg
real(dp) :: deltaE, E_min,E_max,E_plot,eval_point,Sum,k,A,mu,sigma,dx,E_low,E_high,phi_low,phi_high, E_mid, phi_mid, phi_xE
integer  :: gr,dE_plot,start_row,end_row,start_col,end_col,row,col,step,pos,i,j,idx,kl_M,ku_M
integer  :: updated
! timing variables
integer count_0, count_1
integer count_rate, count_max
double precision time_init, time_final, elapsed_time
real(dp), dimension(5000) :: phi_quad_plot

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
phi_un = phi


open(unit=12, file='data.txt')
open(unit=14, file='dataphi.txt')


! Calc mass matrix
call calc_M(M_band,kl_M,ku_M,dE)




!Timing
call system_clock(count_0, count_rate, count_max)
time_init=count_0*1.0/count_rate



! Stepping
do step=1,1

  call build_G_band(no_dof,E_bounds,dE,G_band,kl_G,ku_G)
  !Do single CN step

  call CN_1step(no_dof,M_band,kl_M,ku_M,G_band,kl_G,ku_G,phi_old,phi,dx)

  ! Time copy
  phi_old = phi
enddo

call system_clock(count_1, count_rate, count_max)
time_final = count_1*1.0/count_rate
elapsed_time = time_final-time_init

print *, "number of steps: ", steps
print *,"elasped time: ",elapsed_time



! Plot flux

call project_phi(no_grps-int_size, int_size, phi_bounded_old, phi_proj)

do gr=no_grps,1,-1
    E_low  = E_bounds(gr+1)
    E_high = E_bounds(gr)
    phi_low  = phi_un(3*(gr-1)+1) - phi_un(3*(gr-1)+2) - phi_un(3*(gr-1)+3)
    phi_high = phi_un(3*(gr-1)+1) + phi_un(3*(gr-1)+2) + phi_un(3*(gr-1)+3)
enddo


do gr=no_grps,1,-1
    E_low = E_bounds(gr+1)
    E_high = E_bounds(gr)

    do dE_plot=1, 9, 1
        E_plot = E_low + 0.04_dp*dE_plot
    !    print *, E_plot
        call phi_at_E(gr, E_plot, E_high, E_low, phi, phi_xE)
        write(12,*) phi_xE
    enddo

enddo    

close(12)
close(14)


end program test
