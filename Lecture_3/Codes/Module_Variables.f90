!==================================================================================
!----------------------------------------------------------------------------------
module Mod_UAD2d
    implicit none
    real(8), dimension(:, :), allocatable :: U, V, phi, phi_exat
    real(8) :: D, pi, Error

    real(8) :: EndTime, time, dt, CFL
    integer :: MaxStep, iTimeStep

    real(8) :: tprint, dtprint
    integer :: nprint

    integer :: TimeScheme ! 0 for RK-2, 
                          ! 1 for RK-3, 
                          ! 2 for CN,  unconditional stable
                          ! 3 for Lax-Wendroff, unconditional stable

    integer :: TestCase   ! 0 for Vortex, 
                          ! 1 for Gauss pluse 
    
    character(len=40) :: out_path
end module Mod_UAD2d
!==================================================================================
!----------------------------------------------------------------------------------
module Mod_Grid
    implicit none
    integer :: Nx, Ny, is, ie, js, je
    real(8) :: x1, x2, hx, xLength, y1, y2, hy, yLength
    real(8), dimension(:), allocatable :: x, y 
end module Mod_Grid
