!==================================================================================
!----------------------------------------------------------------------------------
module Mod_UAD1d
    implicit none
    integer :: Nx, is, ie
    real(8) :: x1, x2, hx, Length
    real(8) :: U, D

    real(8) :: EndTime, time, dt, CFL
    integer :: MaxStep, iTimeStep

    real(8) :: tprint, dtprint
    integer :: nprint

    integer :: TimeScheme ! 0 for FTCS, conditional stable
                          ! 1 for BTCS, unconditional stable
                          ! 2 for CN.   unconditional stable
    character(len=40) :: out_path

    real(8), dimension(:), allocatable :: f, fo

end module Mod_UAD1d
