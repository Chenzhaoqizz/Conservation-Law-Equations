!==================================================================================
!----------------------------------------------------------------------------------
module Mod_Poisson2d
    implicit none

    integer :: MaxIter, nIter
    real(8) :: Error, MaxError

    integer :: IterScheme 

    character(len=40) :: out_path
    
    real(8), dimension(:,:), allocatable :: phi, phio
    real(8), dimension(:), allocatable :: phi1d, phio1d
    real(8), dimension(:), allocatable :: AW, AE, AS, AN, P

end module Mod_Poisson2d

!==================================================================================
!----------------------------------------------------------------------------------
module Mod_Grid
    implicit none
    
    integer :: Nx, Ny, is, ie, js, je, NIJ
    real(8) :: xLeft, xRight, yBottom, yTop, xLength, yLength
    real(8), dimension(:), allocatable :: x, y, deltaX, deltaY

    logical(4) :: Is_Uniform_Grid 

end module Mod_Grid