module lapack_operations

contains

subroutine band_matvec(n,ab,kl,ku,alpha,x,y,add_to)
use f90_kind
implicit none
!
! Perform y = (y +) alpha * ab*x (dependig on the 'add_to_y' input
! where ab is a banded matrix stored in the Lapack format
!
integer, intent(in) :: n
real(dp), dimension(:,:), intent(in) :: ab
integer, intent(in) :: kl
integer, intent(in) :: ku
real(dp), intent(in) :: alpha
real(dp), dimension(:), intent(in) :: x
real(dp), dimension(:), intent(inout) :: y
logical, intent(in) :: add_to

integer :: row_diag,diag,row,x_size,y_size,ab_size_1,ab_size_2

! Perform checks on sizes of ab and x and y

x_size = size(x)
y_size = size(y)
ab_size_1 = size(ab,1)
ab_size_2 = size(ab,2)

!print *, x_size, y_size, ab_size_1, ab_size_2

if (x_size < n) STOP 'band_matvec: x too short'
if (y_size < n) STOP 'band_matvec: y too short'
if (ab_size_1 < 2*kl + ku + 1) STOP 'band_matvec: ab too small'
if (ab_size_2 < n) STOP 'band_matvec: ab too small'

! Init

if (.not. add_to) y = 0.0_dp

row_diag = kl+ku+1

! diagonal

row = row_diag
y(1:n) = y(1:n) + alpha * ab(row,1:n) * x(1:n)

! super-diagonals

do diag=1,ku
  row = row_diag - diag
  y(1:n-diag) = y(1:n-diag) + alpha * ab(row,1+diag:n) * x(1+diag:n)
enddo

! sub-diagonals

do diag=1,kl
  row = row_diag + diag
  y(1+diag:n) = y(1+diag:n) + alpha * ab(row,1:n-diag) * x(1:n-diag)
enddo

end subroutine band_matvec


subroutine band_mat_add(n,ab1,kl1,ku1,ab2,kl2,ku2,alpha,add_to)
use f90_kind
implicit none
!
! Perform ab1 = (ab1 +) alpha * ab2 (dependig on the 'add_to_y' input
! where abx are banded matrices stored in the Lapack format
! Note that the locations of ab1 that do not overlap with ab2 are not touched
!
integer, intent(in) :: n
real(dp), dimension(:,:), intent(inout) :: ab1
integer, intent(in) :: kl1
integer, intent(in) :: ku1
real(dp), dimension(:,:), intent(in) :: ab2
integer, intent(in) :: kl2
integer, intent(in) :: ku2
real(dp), intent(in) :: alpha
logical, intent(in) :: add_to

integer :: row_diag,row_start,row_end

! Perform checks on sizes of ab and x and y

if (kl1 < kl2) STOP 'band_mat_add: kl2 > kl1'
if (ku1 < ku2) STOP 'band_mat_add: ku2 > ku1'

row_diag = kl1+ku1+1

row_start = row_diag - ku2
row_end   = row_diag + kl2

! Init

if (.not. add_to) then
  ab1(row_start:row_end,1:n) = 0.0_dp
endif

! ab1 = ab1 + alpha * ab2

ab1(row_start:row_end,1:n) = ab1(row_start:row_end,1:n) + alpha * ab2

end subroutine band_mat_add

end module lapack_operations
