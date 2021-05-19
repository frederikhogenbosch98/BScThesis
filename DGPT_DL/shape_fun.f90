module shape_fun_E

contains

subroutine calc_shape_fun_E(E,E_min,E_max,no_nod_E,scaling,fun_E)
! Note: to get the normalization correct we have divided by dE and have a factor
! 3 for the slope. This routine only to be used as multiplier for the source.
! It can not be used for integrals of the form fi_x_fj ...
! scaling=F: no scaling
! scaling=T:    scaling
use f90_kind
implicit none

real(dp), intent(in) :: E
real(dp), intent(in) :: E_min
real(dp), intent(in) :: E_max
integer, intent(in) :: no_nod_E
real(dp), allocatable, dimension(:), intent(out) :: fun_E
logical :: scaling

real(dp) :: dE,E_mid

allocate(fun_E(no_nod_E))
if (no_nod_E == 1) then
  fun_E = 1.0_dp
else if (no_nod_E == 2) then
    dE = E_max - E_min
    E_mid = E_min + dE/2.0_dp
    fun_E(1) = 1.0_dp              ! average
    fun_E(2) = (2.d0/dE)*(E-E_mid) ! slope
    if (scaling) then
        fun_E(1) = fun_E(1) * (1.0_dp/dE)
        fun_E(2) = fun_E(2) * (3.0_dp/dE)
    endif
else
    dE = E_max - E_min
    E_mid = E_min + dE/2.0_dp
    fun_E(1) = 1.0_dp
    fun_E(2) = 2.d0/dE*(E-E_mid)
    fun_E(3) = ((3.0_dp/2.0_dp)*((2.0_dp*(E-E_mid))/dE)**2)-(1.0_dp/2.0_dp)
endif

end subroutine calc_shape_fun_E

end module shape_fun_E
