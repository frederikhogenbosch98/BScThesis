module fem

contains

subroutine calc_M(M_band,kl_M,ku_M,dE)
use f90_kind
implicit none

real(dp), dimension(:,:), allocatable, intent(inout) :: M_band
integer, intent(out) :: kl_M
integer, intent(out) :: ku_M
real(dp), dimension(:), intent(in) :: dE

integer :: gr,no_grps,row,col,row_band,col_band,N
real(dp) :: M11,M22

kl_M = 0
ku_M = 0

no_grps = size(dE)
N = 2 * no_grps

if (allocated(M_band)) deallocate(M_band)
allocate(M_band(2*kl_M+ku_M+1,N))

M_band = 0.0_dp

do gr=1,no_grps
  M11 = dE(gr)
  M22 = dE(gr) / 3.0_dp
!print *, M11
  row = (gr-1)*2 + 1
  col = (gr-1)*2 + 1

  row_band = kl_M + ku_M + 1 + row - col
  col_band = col

  M_band(row_band,col_band) = M11

  row = (gr-1)*2 + 2
  col = (gr-1)*2 + 2

  row_band = kl_M + ku_M + 1 + row - col
  col_band = col

  M_band(row_band,col_band) = M22
enddo
!print *, size(M_band)
end subroutine calc_M

subroutine project_gaussian_on_dg(A,mu,sigma,E_bounds,dE,phi)
use f90_kind
use quadrature
use functions
use shape_fun_E
implicit none

real(dp), intent(in) :: A
real(dp), intent(in) :: mu
real(dp), intent(in) :: sigma
real(dp), dimension(:), intent(in) :: E_bounds
real(dp), dimension(:), intent(in) :: dE
real(dp), dimension(:), intent(out) :: phi

integer :: no_grps,qp,gr
real(dp) :: rhs(2),E_low,E_high,Jac,M11,M22,E_qp
real(dp),allocatable, dimension(:) :: fun_E
integer, parameter :: nqp=7
real(dp), dimension(nqp) :: points,weights

call legendre_set(points,weights)

phi = 0.0_dp
!print *,dE
no_grps = size(dE)
do gr=1,no_grps
  E_low  = E_bounds(gr+1)
  E_high = E_bounds(gr)
  M11 = dE(gr)
  M22 = dE(gr) / 3.0_dp

  rhs = 0.0_dp
  Jac = dE(gr) / 2.0_dp

  do qp=1,nqp
    E_qp = E_low + ((points(qp) + 1.0_dp) / 2.0_dp) * dE(gr)

    call calc_shape_fun_E(E_qp,E_low,E_high,2,.false.,fun_E)

    rhs = rhs + Jac * weights(qp) * gaussian(E_qp,A,mu,sigma) * fun_E
  enddo

  phi(2*(gr-1)+1) = rhs(1) / M11
  phi(2*(gr-1)+2) = rhs(2) / M22
enddo
!print *, phi
end subroutine project_gaussian_on_dg


subroutine det_bounds(E_bounds, dE, old_phi, E_bounded, phi_bounded, n_max, int_size)

use quadrature
use functions
use f90_kind
implicit none

real(dp), dimension(:), intent(in) :: E_bounds
real(dp), dimension(:), intent(in) :: dE
!real(dp), intent(in) :: E_max
!real(dp), intent(in) :: E_min
real(dp), dimension(:), intent(in) :: old_phi
real(dp), dimension(:), intent(out) :: E_bounded
real(dp), dimension(:), intent(out) :: phi_bounded
integer, intent(in) :: int_size
integer :: n_max
intrinsic :: findloc
integer :: x(1)
real(dp) :: trial
!print *,E_min, E_max
E_bounded = E_bounds(n_max:n_max+int_size-1)
phi_bounded = old_phi(2*n_max:(2*n_max)+(2*int_size)-1)
!print *, 2*n_max, (2*n_max)+79
!print *, ((2*n_max)+80)-(2*n_max)
!I = pack([(j, j=1, size(E_bounds))],E_bounds==E_max)
!trial = E_bounds(1)-(dE(12)*12.0_dp)
!print *,trial
!x = findloc(E_bounds, value=trial)

!print *, E_bounded,phi_bounded
!print *, E_bounds_new
!print *, E_bounds([5,15])



end subroutine det_bounds

subroutine update_bounds(E_bounds, phi_old, E_bounds_old, dE, phi_bounded_old, E_bounded, phi_bounded, step, n_max, int_size,no_steps, updated)
use quadrature
use functions
use f90_kind
implicit none

real(dp), dimension(:), intent(in) :: E_bounds
real(dp), dimension(:), intent(in) :: phi_old
real(dp), dimension(:), intent(in) :: E_bounds_old
real(dp), dimension(:), intent(in) :: dE
real(dp), dimension(:), intent(in) :: phi_bounded_old
real(dp), dimension(:), intent(out) :: E_bounded
real(dp), dimension(:), intent(out) :: phi_bounded
integer, intent(in) :: step
integer, intent(in) :: n_max
integer, intent(in) :: int_size
integer, intent(in) :: no_steps
integer, intent(in) :: updated
integer :: interval
integer :: divider
integer :: rest_length
integer :: slowing_fact
integer :: fill
real(dp), dimension(size(phi_old)) :: phi_proj
integer :: bounded_f
integer :: new_interval

slowing_fact = 1.0_dp
rest_length = n_max - int_size
divider = no_steps / 10
!print *, (step/200)*(rest_length/divider)/slowing_fact
!interval = n_max + (step/100)*(rest_length / divider)*slowing_fact

new_interval = n_max + updated*5
interval = n_max + (updated-1)*5
print *, interval, size(E_bounds)-int_size
if (interval>(size(E_bounds)-int_size-40)) then
    interval = size(E_bounds)-1-int_size
    new_interval = size(E_bounds)-1-int_size
    print *, "larger"
endif
!print *, interval, interval+int_size-1
!print *, "E bounds", interval, interval+int_size-1
print *, "phi bounds", 2*interval, (2*interval)+(2*int_size)-1
E_bounded = E_bounds(interval:interval+int_size-1)
!print *, E_bounded
do fill=1,size(phi_proj)
    bounded_f = fill-2*interval
    if (fill<2*interval) then
        phi_proj(fill) = 0.0_dp
    else if (fill>2*interval .AND. fill<((2*interval)+(2*int_size))) then
        phi_proj(fill) = phi_bounded_old(bounded_f)
    else
        phi_proj(fill) = 0.0_dp
    end if

enddo
!print *, phi_proj
!print *, "int: ", interval
!print *, "new int: ", new_interval
!print *, phi_proj
!phi_proj = phi_bounded_old
phi_bounded = phi_proj((2*new_interval):(2*new_interval)+((2*int_size)-1))
!phi_bounded = phi_bounded_old
!E_bounded = E_bounds(n_max
end subroutine    


subroutine build_G_band(n,E_bounds,dE,G_band,kl_G,ku_G, step)
use f90_kind
use heterogeneous
implicit none

integer :: n
real(dp), dimension(:), intent(in) :: E_bounds
real(dp), dimension(:), intent(in) :: dE
real(dp), dimension(:,:), allocatable, intent(inout) :: G_band
integer :: kl_G
integer :: ku_G
integer, intent(in) :: step

real(dp), parameter :: penalty = 2.0_dp
character(len=1) :: mat
integer :: gr,i,j,row,col,row_band,col_band,no_grps,n_min,limit,n_max,E_step
real(dp) :: E_low,E_high,E_g,S_A,S_E,S_low,S_high,min_dE
real(dp) :: T_low,T_high,T_avg
real(dp) :: A_group(2,2)
real(dp) :: A_low_E(2,2)
real(dp) :: A_high_E(2,2)

if (allocated(G_band)) deallocate(G_band)
allocate(G_band(2*kl_G+ku_G+1,n))
G_band = 0.0_dp

mat="W"
limit=10
no_grps = size(dE)
E_step = NINT(real(step/4))

if (E_step<limit+1) then
    n_min = 1
else 
    n_min = E_step-limit
end if

if (E_step>(size(dE)-limit)) then
    n_max = size(dE)
else 
    n_max = E_step+limit
end if

!print *,"nmin = ", n_min
!print *, "nmax = ", n_max


do gr=1,no_grps 
!do gr=n_min,n_max 
  E_low  = E_bounds(gr+1)
  E_high = E_bounds(gr)
  E_g = (E_low + E_high)/2.0_dp

  T_low  = T(E_low,mat) 
  T_high = T(E_high,mat)
  T_avg = (T_low + T_high)/2.0_dp
  ! Here the factor 1/2 is incorporated

  T_low  = T_low  / 2.0_dp
  T_high = T_high / 2.0_dp
  T_avg  = T_avg  / 2.0_dp

  S_low  = S(E_low,mat)
  S_high = S(E_high,mat)

  S_A    = (S_high + S_low) / 2.0_dp
  S_E    = (S_high - S_low) / 2.0_dp

  ! Init

  A_group  = 0.0_dp
  A_high_E = 0.0_dp
  A_low_E  = 0.0_dp

  ! Average equation

  A_group(1,1)  =  A_group(1,1)  + S_low
  A_group(1,2)  =  A_group(1,2)  - S_low
  A_high_E(1,1) =  A_high_E(1,1) - S_high
  A_high_E(1,2) =  A_high_E(1,2) + S_high

  ! Slope equation

  A_group(2,1)  = A_group(2,1)  - S_low + 2.0_dp * S_A
  A_group(2,2)  = A_group(2,2)  + S_low + 2.0_dp / 3.0_dp * S_E
  A_high_E(2,1) = A_high_E(2,1) - S_high
  A_high_E(2,2) = A_high_E(2,2) + S_high
            
  ! Volume term
      
  A_group(2,2) = A_group(2,2) + T_avg * 4.0_dp / dE(gr)

  ! Penalty term on high-E side (g-1/2)

  if (gr /= 1) then
  !if (gr /= n_min) then
    min_dE = min(dE(gr),dE(gr-1))
    A_group(1,1)  = A_group(1,1)  + penalty * T_high / min_dE
    A_group(1,2)  = A_group(1,2)  + penalty * T_high / min_dE
    A_group(2,1)  = A_group(2,1)  + penalty * T_high / min_dE
    A_group(2,2)  = A_group(2,2)  + penalty * T_high / min_dE

    A_high_E(1,1) = A_high_E(1,1) - penalty * T_high / min_dE
    A_high_E(1,2) = A_high_E(1,2) + penalty * T_high / min_dE
    A_high_E(2,1) = A_high_E(2,1) - penalty * T_high / min_dE
    A_high_E(2,2) = A_high_E(2,2) + penalty * T_high / min_dE
  endif

  ! Penalty term on low-E side (g+1/2)

  if (gr /= no_grps) then
  !if (gr /= n_max) then  
    min_dE = min(dE(gr),dE(gr+1))
    A_group(1,1) = A_group(1,1) + penalty * T_low / min_dE
    A_group(1,2) = A_group(1,2) - penalty * T_low / min_dE
    A_group(2,1) = A_group(2,1) - penalty * T_low / min_dE
    A_group(2,2) = A_group(2,2) + penalty * T_low / min_dE

    A_low_E(1,1) = A_low_E(1,1) - penalty * T_low / min_dE
    A_low_E(1,2) = A_low_E(1,2) - penalty * T_low / min_dE
    A_low_E(2,1) = A_low_E(2,1) + penalty * T_low / min_dE
    A_low_E(2,2) = A_low_E(2,2) + penalty * T_low / min_dE
  endif

  ! consistency/symmetry terms on high-E side (g-1/2)

  if (gr /= 1) then
  !if (gr /= n_min) then
    A_group(1,2)  = A_group(1,2)  - T_high/dE(gr)
    A_group(2,1)  = A_group(2,1)                    - T_high/dE(gr)
    A_group(2,2)  = A_group(2,2)  - T_high/dE(gr)   - T_high/dE(gr)

    A_high_E(1,2) = A_high_E(1,2) - T_high/dE(gr-1)
    A_high_E(2,1) = A_high_E(2,1)                   + T_high/dE(gr)
    A_high_E(2,2) = A_high_E(2,2) - T_high/dE(gr-1) - T_high/dE(gr)
  endif

  ! consistency/symmetry terms on low-E side  (g+1/2)

  if (gr /= no_grps) then
   !if (gr /= n_max) then
    A_group(1,2) = A_group(1,2) + T_low/dE(gr)
    A_group(2,1) = A_group(2,1)                  + T_low/dE(gr)
    A_group(2,2) = A_group(2,2) - T_low/dE(gr)   - T_low/dE(gr)

    A_low_E(1,2) = A_low_E(1,2) + T_low/dE(gr+1)
    A_low_E(2,1) = A_low_E(2,1)                  - T_low/dE(gr)
    A_low_E(2,2) = A_low_E(2,2) - T_low/dE(gr+1) - T_low/dE(gr)
  endif

  do i=1,2
    row = (gr-1)*2 + i
    do j=1,2
      ! Fill diag block into G

      col = (gr-1)*2 + j
      row_band = kl_G + ku_G + 1 + row - col
      col_band = col

      G_band(row_band,col_band) = A_group(i,j)

      ! Fill left block into G

      if (gr > 1) then
        col = (gr-2)*2 + j
        row_band = kl_G + ku_G + 1 + row - col
        col_band = col

        G_band(row_band,col_band) = A_high_E(i,j)
      endif

      ! Fill right block into G

      if (gr < no_grps) then
        col = (gr)*2 + j
        row_band = kl_G + ku_G + 1 + row - col
        col_band = col

        G_band(row_band,col_band) = A_low_E(i,j)
      endif
    enddo
  enddo
enddo

end subroutine build_G_band



end module fem
