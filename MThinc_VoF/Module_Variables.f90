!--------------------------------------------------------------------------------------
!======================================================================================
module Mod_grid
    implicit none 
    integer :: Nx, Ny, is, ie, js, je 
    real(8) :: xLength, yLength, hxi, hyi
    
    real(8) :: Start_time, End_Time, CFL, time, delta_t 
    integer :: MaxStep, iStep 

    real(8) :: tprint, dtprint 
    integer :: nprint, iproblem

    integer :: MaxIter 
    real(8) :: beta, Error, pi 

    logical :: restart 
    character(len=40) :: output_path
    integer :: outputFormat

end module Mod_grid 
!--------------------------------------------------------------------------------------
!======================================================================================
module Mod_flow
    implicit none
    real(8), parameter :: eplison = 1.0E-8
    real(8), dimension(:,:), allocatable :: u, v, phi, phio 
    real(8), dimension(:,:), allocatable :: flux_x, flux_y
    real(8), dimension(:,:), allocatable :: NormVectorX, NormVectorY, Lxx, Lyy, Lxy 

end module Mod_flow
!--------------------------------------------------------------------------------------
!======================================================================================