!--------------------------------------------------------------------------------------
!                   define variables here
!======================================================================================
module flow_var
    implicit none
    integer :: Nx, Ny, is, ie, js, je
    real(8) :: xLeft, xRight, yBottom, yTop
    real(8), allocatable :: delta_x(:), x(:), delta_y(:), y(:)
    real(8), allocatable :: u(:,:), v(:,:), vort(:,:), sf(:,:), p(:,:)

    ! variables to solve possion equation of stream function
    real(8), allocatable :: AP(:,:), AE(:,:), AW(:,:), AN(:,:), AS(:,:)

    real(8) :: CFL, Epsilon, MaxIter
    real(8) :: Start_time, End_Time, delta_t, time
    real(8) :: uTop, Re
    integer :: iflux_splitting, ischeme_inv, nprint
    integer :: isave, istep, imesh, iPossion_Equation

    real(8), parameter :: pi = 4.d0 * atan(1.d0)
    
contains
    subroutine alloc_data(Nx, Ny)
        implicit none
        integer, intent(in) :: Nx, Ny

        print*, 'Allocating data... '
    
        allocate(delta_x(Nx), x(Nx), delta_y(Ny), y(Ny))
        allocate(u(Nx,Ny), v(Nx,Ny), vort(Nx,Ny), sf(Nx,Ny), p(Nx,Ny))
        allocate(AP(Nx,Ny), AE(Nx,Ny), AW(Nx,Ny), AN(Nx,Ny), AS(Nx,Ny))

    end subroutine
end module flow_var
!--------------------------------------------------------------------------------------
