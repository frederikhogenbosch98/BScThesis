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
real(dp) :: M11,M22,M33

kl_M = 0
ku_M = 0

no_grps = size(dE)
N = 3 * no_grps

if (allocated(M_band)) deallocate(M_band)
allocate(M_band(2*kl_M+ku_M+1,N))

M_band = 0.0_dp

do gr=1,no_grps
  M11 = dE(gr)
  M22 = dE(gr) / 3.0_dp
  M33 = dE(gr) / 5.0_dp
!print *, M11
  row = (gr-1)*3 + 1
  col = (gr-1)*3 + 1
!  print *, 'col M11: ', col
  row_band = kl_M + ku_M + 1 + row - col
  col_band = col

  M_band(row_band,col_band) = M11

  row = (gr-1)*3 + 2
  col = (gr-1)*3 + 2
!    print *, 'col M22: ', col
  row_band = kl_M + ku_M + 1 + row - col
  col_band = col

  M_band(row_band,col_band) = M22

  row = (gr-1)*3 + 3
  col = (gr-1)*3 + 3
!  print *, 'col M33: ',col
  row_band = kl_M + ku_M + 1 + row - col
  col_band = col
  M_band(row_band,col_band) = M33

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
real(dp) :: rhs(3),E_low,E_high,Jac,M11,M22,M33,E_qp
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
  M33 = dE(gr) / 5.0_dp 
  rhs = 0.0_dp
  Jac = dE(gr) / 2.0_dp

  do qp=1,nqp
    E_qp = E_low + ((points(qp) + 1.0_dp) / 2.0_dp) * dE(gr)

    call calc_shape_fun_E(E_qp,E_low,E_high,3,.false.,fun_E)
    rhs = rhs + Jac * weights(qp) * gaussian(E_qp,A,mu,sigma) * fun_E
     
    
  enddo

  phi(3*(gr-1)+1) = rhs(1) / M11
  phi(3*(gr-1)+2) = rhs(2) / M22
  phi(3*(gr-1)+3) = rhs(3) / M33
enddo
!print *, phi
end subroutine project_gaussian_on_dg


subroutine det_bounds(E_bounds, dE, phi, E_bounded, phi_bounded, n_max, int_size)

use quadrature
use functions
use f90_kind
implicit none

real(dp), dimension(:), intent(in) :: E_bounds
real(dp), dimension(:), intent(in) :: dE
!real(dp), intent(in) :: E_max
!real(dp), intent(in) :: E_min
real(dp), dimension(:), intent(in) :: phi
real(dp), dimension(:), intent(out) :: E_bounded
real(dp), dimension(:), intent(out) :: phi_bounded
integer, intent(in) :: int_size
integer :: n_max
intrinsic :: findloc
integer :: x(1)
real(dp) :: trial
!print *,E_min, E_max
E_bounded = E_bounds(n_max+1:n_max+1+int_size)
phi_bounded = phi(3*n_max+1:(3*n_max+1)+(3*int_size)-1)
!print *, phi(480:520)
!phi_bounded = phi(480:529)
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

subroutine project_phi(interval, int_size, phi_bounded_old, phi_proj)

use functions
use f90_kind
implicit none

integer, intent(in) :: interval
integer, intent(in) :: int_size
real(dp), dimension(:), intent(in) :: phi_bounded_old
real(dp), dimension(:), intent(out) :: phi_proj
integer :: fill, bounded_f


do fill=1,size(phi_proj)
    bounded_f = fill-3*interval
    if (fill<3*interval) then
        phi_proj(fill) = 0.0_dp
    else if (fill>3*interval .AND. fill<((3*interval)+(3*int_size))) then
        phi_proj(fill) = phi_bounded_old(bounded_f)
    else
        phi_proj(fill) = 0.0_dp
    end if
enddo


end subroutine


subroutine update_bounds(E_bounds, phi_old, E_bounds_old, dE, phi_bounded_old, E_bounded, phi_bounded, step, n_max,int_size,no_steps, updated,phi_proj)
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
real(dp), dimension(:), intent(out) :: phi_proj
integer :: fill
real(dp), dimension(size(E_bounds)) :: phi_high
integer :: interval
!real(dp), dimension(size(phi_old)) :: phi_proj
integer :: bounded_f
integer :: new_interval
integer :: lp



new_interval = n_max + updated*(int_size/3)
interval = n_max + (updated-1)*(int_size/3)
!print *, interval, size(E_bounds)-int_size
if (interval>(size(E_bounds)-int_size-40)) then
    interval = size(E_bounds)-int_size
    print *, 'at the end'
    new_interval = size(E_bounds)-int_size
endif


E_bounded = E_bounds(interval+1:interval+int_size+1)

call project_phi(interval, int_size, phi_bounded_old, phi_proj)

do lp=size(E_bounds)-1,1,-1
    phi_high(lp) = phi_proj(3*(lp-1)+1)-phi_proj(3*(lp-1)+2)
enddo

write (14,*) phi_high


phi_bounded = phi_proj((3*new_interval+1):(3*new_interval+1)+((3*int_size)-1))


end subroutine    


subroutine build_G_band(n,E_bounds,dE,G_band,kl_G,ku_G)
use f90_kind
use heterogeneous
implicit none

integer :: n
real(dp), dimension(:), intent(in) :: E_bounds
real(dp), dimension(:), intent(in) :: dE
real(dp), dimension(:,:), allocatable, intent(inout) :: G_band
integer :: kl_G
integer :: ku_G

real(dp), parameter :: penalty = 2.0_dp
character(len=1) :: mat
integer :: gr,i,j,row,col,row_band,col_band,no_grps
real(dp) :: E_low,E_high,E_g,S_A,S_E,S_low,S_high,min_dE
real(dp) :: T_low,T_high,T_avg
real(dp) :: A_group(2,2)
real(dp) :: A_low_E(2,2)
real(dp) :: A_high_E(2,2)

if (allocated(G_band)) deallocate(G_band)
allocate(G_band(2*kl_G+ku_G+1,n))
G_band = 0.0_dp

mat="W"
no_grps = size(dE)

do gr=1,no_grps 
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


subroutine plot(E, phi, E_low, E_high, phi_low, phi_high, E_avg, phi_avg)

use f90_kind
implicit none
real(dp), dimension(:), intent(in) :: E
real(dp), dimension(:), intent(in) :: phi
real(dp), dimension(:), intent(out) :: E_low
real(dp), dimension(:), intent(out) :: E_high
real(dp), dimension(:), intent(out) :: phi_low
real(dp), dimension(:), intent(out) :: phi_high
real(dp), dimension(:), intent(out) :: E_avg
real(dp), dimension(:), intent(out) :: phi_avg
integer :: gr
!print *, phi

do gr=size(E)-1,1,-1
    !E_low(gr) = E(gr+1)
    !E_high(gr) = E(gr)
    phi_low(gr) = phi(2*(gr-1)+1) - phi(2*(gr-1)+2)
    phi_high(gr) = phi(2*(gr-1)+1) + phi(2*(gr-1)+2)
enddo

end subroutine

!subroutine return_p(E, E_avg, E_mid, p_func)
!use f90_kind
!implicit none

!real(dp), intent(in) :: E
!real(dp), intent(in) :: E_avg
!real(dp), intent(in) :: E_mid
!real(dp), allocatable, dimension(:), intent(out) :: p_func

!p_func(1) = 1.0_dp
!p_func(2) = 2.d0/E_avg*(E-E_mid)
!p_func(3) = ((3.0_dp/2.0_dp)*((2.0_dp*(E-E_mid))/E_avg)**2)-(1.0_dp/2.0_dp)

!print *, 'enter'
!end subroutine


subroutine phi_at_E(iter_coef, E, E_max, E_min, phi_coeff, phi_xE)

use f90_kind
use shape_fun_E
implicit none

integer, intent(in) :: iter_coef
real(dp), intent(in) :: E
real(dp), intent(in) :: E_max
real(dp), intent(in) :: E_min
real(dp), dimension(:), intent(in) :: phi_coeff
real(dp), intent(out) :: phi_xE
real(dp) :: E_avg
real(dp) :: E_mid
real(dp), dimension(3) :: p_func
integer :: phi_loc
E_avg = E_max - E_min
E_mid = E_min + E_avg/2

call return_p(E, E_avg, E_mid, p_func)

print *, E, E_avg, E_mid
print *, p_func(3)
!phi_loc = 3*NINT(E_min/dE)
phi_loc = iter_coef
phi_xE = phi_coeff(phi_loc-2)*p_func(1)+phi_coeff(phi_loc-1)*p_func(2)+phi_coeff(phi_loc)*p_func(3)
!print *,phi_coeff(phi_loc), phi_coeff(phi_loc+1)
print *, phi_loc, phi_loc-1, phi_loc-2
end subroutine


end module fem
