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
integer, parameter :: n_max=241-12
!integer, parameter :: n_max=24-12
integer,parameter   :: kl_G=3
integer,parameter   :: ku_G=3
real(dp),allocatable,dimension(:) :: point,weight
real(dp), dimension(:,:), allocatable :: M_band
real(dp), dimension(:,:), allocatable :: G_band

integer, parameter :: no_dof=2*no_grps
real(dp), dimension(no_dof) :: phi,phi_old,phi_un
real(dp) :: dE(no_grps)
real(dp) :: phi_plot(no_grps)
real(dp) :: dE_bounded(int_size)
real(dp) :: E_bounds(no_grps+1)
real(dp) :: E_bounded(int_size+1)
real(dp) :: E_bounded_old(int_size+1)
real(dp) :: phi_bounded(int_size*2)
real(dp) :: phi_bounded_old(int_size*2)
real(dp) :: phi_bounded_un(int_size*2)
real(dp) :: phi_non_CN(int_size*2)
real(dp), dimension(size(phi_old)) :: phi_proj
real(dp), dimension(int_size) :: E_lowb, E_highb, E_avg
real(dp), dimension(int_size) :: phi_lowb, phi_highb, phi_avg
real(dp) :: E_min,E_max,eval_point,Sum,k,A,mu,sigma,dx,E_low,E_high,phi_low,phi_high
integer  :: gr,grdos,start_row,end_row,start_col,end_col,row,col,step,pos,i,j,idx,kl_M,ku_M
integer  :: updated, interval, new_interval
real(dp) :: phi_boundary
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

dx = x_max / no_steps


! Init phi
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
dE_bounded = dE(1:int_size)
call calc_M(M_band,kl_M,ku_M,dE_bounded)


! Init E and phi bounds
call det_bounds(E_bounds, dE, phi, E_bounded, phi_bounded, n_max, int_size)
phi_bounded_old = phi_bounded
phi_bounded_un = phi_bounded
updated = 0
E_bounded_old = E_bounded



!Timing
call system_clock(count_0, count_rate, count_max)
time_init=count_0*1.0/count_rate

! Init interval
interval = n_max


! Stepping
do step=1,2000

  ! Define update boundary
  if (step<1500) then
      phi_boundary = 0.001_dp
  else
      phi_boundary = 0.0023_dp
  endif

  ! Update E and phi bounds
  if (phi_bounded(int_size*2)> phi_boundary) then
    print *, 'updated at step: ', step
    call update_bounds(E_bounds, phi_old, E_bounded_old, dE, phi_bounded_old, E_bounded, phi_bounded, interval, step, n_max, int_size,no_steps, updated, phi_proj, new_interval)    
    updated = step
    phi_non_cn = phi_bounded
    interval = new_interval
  endif

  ! Plot phi at step
    if(step==2000) then
      do gr=int_size,1,-1
        phi_high = phi_bounded(2*(gr-1)+1) + phi_bounded(2*(gr-1)+2)
      enddo
  endif

  ! Plot phi every 100 steps
  if (mod(step,100)==0) then
      !print *, step
      call project_phi(interval, int_size, phi_bounded_old, phi_proj)
      do gr=no_grps,1,-1
        phi_plot(gr)=phi_proj(2*(gr-1)+1)+phi_proj(2*(gr-1)+2)
      enddo
        write(14,*) phi_plot
  endif

  phi_bounded_old = phi_bounded
  
  ! Construct G matrix
  call build_G_band(size(phi_bounded),E_bounded,dE_bounded,G_band,kl_G,ku_G)

  ! Do single CN step
  call CN_1step(size(phi_bounded),M_band,kl_M,ku_M,G_band,kl_G,ku_G,phi_bounded_old,phi_bounded,dx)

  ! Copy
  phi_old = phi
  phi_bounded_old = phi_bounded
enddo

! End system clock
call system_clock(count_1, count_rate, count_max)
time_final = count_1*1.0/count_rate
elapsed_time = time_final-time_init

print *,"elasped time: ",elapsed_time

! Project phi
call project_phi(no_grps-int_size, int_size, phi_bounded_old, phi_proj)

! Plot flux
do gr=no_grps,1,-1
    E_low  = E_bounds(gr+1)
    E_high = E_bounds(gr)
    phi_low  = phi_proj(2*(gr-1)+1) - phi_proj(2*(gr-1)+2)
    phi_high = phi_proj(2*(gr-1)+1) + phi_proj(2*(gr-1)+2)
enddo


do gr=int_size,1,-1
    phi_low = phi_bounded_un(2*(gr-1)+1) - phi_bounded_un(2*(gr-1)+2)
    phi_high = phi_bounded_un(2*(gr-1)+1) + phi_bounded_un(2*(gr-1)+2)
enddo


close(14)
close(12)


end program test
