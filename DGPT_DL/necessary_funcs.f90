module functions
  ! This module contains all the functions that are needed throughout
  ! the code

  use f90_kind
  use physconst
  implicit none 

contains

  real(dp) function gaussian(E,A,mu,sigma)
    ! Calculates the flux based on the amplitude, position and spread
    ! Parameters:
    ! _________________________________________________________
    ! Input: A_x - amplitude of gaussian
    !        E_x - position of gaussian
    !        Sigma_x - spread of gaussian
    !        E       - energy domain array
    ! Output: phi - flux
    ! _________________________________________________________

    real(dp), intent(in) :: E,A,mu,sigma

    gaussian = A * exp(-(E-mu)**2 / sigma**2)

  end function gaussian
  
end module functions 
