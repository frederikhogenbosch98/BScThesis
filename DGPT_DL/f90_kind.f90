!  Module defining the KIND numbers for the NAGWare f90 Compiler.
!  Copyright 1991-1999 The Numerical Algorithms Group Ltd., Oxford, U.K.
!  Malcolm Cohen, Robert Iles, June 1991
!
module f90_kind

!---------------------------------------------------------------------
! This module defines the double precision kind parameter
!---------------------------------------------------------------------

    intrinsic kind,selected_int_kind,selected_real_kind  ! we use these,
    private kind,selected_int_kind,selected_real_kind    ! but do not force
                                                         ! them on the user.
    integer, parameter :: dp      = kind(0.0d0)
end module f90_kind
