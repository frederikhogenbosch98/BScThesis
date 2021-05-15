module CN

contains

subroutine CN_1step(n,M_band,kl_M,ku_M,G_band,kl_G,ku_G,sol_old,sol,dx)
use f90_kind
use lapack_operations
implicit none

! Perform 1 step of Crank-Nicholson on the system
! M * dy/dx + G * y = 0
! Note: G is destroyed in the process, M is kept

real(dp), parameter :: theta = 0.5_dp

integer, intent(in) :: n
real(dp), dimension(:,:), intent(in) :: M_band
integer, intent(in) :: kl_M
integer, intent(in) :: ku_M
real(dp), dimension(:,:), intent(inout) :: G_band
integer, intent(in) :: kl_G
integer, intent(in) :: ku_G
real(dp), dimension(:), intent(in) :: sol_old
real(dp), dimension(:), intent(out) :: sol
real(dp) :: dx

integer :: info
integer, dimension(n) :: ipiv
real(dp), dimension(n) :: rhs
!print *, 'CNstep'
!print *,sol_old

rhs = 0.0_dp
call band_matvec(n,M_band,kl_M,ku_M,1.0_dp/dx,        sol_old,rhs,.false.)
call band_matvec(n,G_band,kl_G,ku_G,-(1.0_dp - theta),sol_old,rhs,.true.)

G_band(kl_G+1:,:) = theta * G_band(kl_G+1:,:)
call band_mat_add(n,G_band,kl_G,ku_G,M_band,kl_M,ku_M,1.0_dp/dx,.true.)

! Solve the system

call dgbsv(n,kl_G,ku_G,1,G_band,size(G_band,1),ipiv,rhs,n,info)
if (info/=0) then
  print *,'WARNING in solver, info=',info
endif

sol = rhs

end subroutine CN_1step

end module CN
