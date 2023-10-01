!==================================================================================
!        FTCS scheme to solve the unsteady advection diffusion equation
!----------------------------------------------------------------------------------
subroutine FTCS_Scheme 
    use Mod_UAD2d
    use Mod_Grid
    implicit none

    integer :: i, j
    real(8) :: Advx, Advy, Fluxx, Fluxy, AdvectionFlux
    real(8) :: Visx, Visy, Diffx, Diffy, DiffusionFlux
    real(8), dimension(Nx, Ny) :: phio

    phio = phi

    Advx = 0.5d0*dt/hx; 
    Advy = 0.5d0*dt/hy;
    Visx = D*dt/(hx)**2; 
    Visy = D*dt/(hy)**2;

    do i = is+1, ie-1
        do j = js+1, je-1
            Fluxx = AdvectionFlux(u(i+1,j), phio(i+1,j), u(i-1,j), phio(i-1,j)) 
            Fluxy = AdvectionFlux(v(i,j+1), phio(i,j+1), v(i,j-1), phio(i,j-1)) 
            Diffx = DiffusionFlux(phio(i+1,j), phio(i,j), phio(i-1,j))
            Diffy = DiffusionFlux(phio(i,j+1), phio(i,j), phio(i,j-1))

            phi(i,j) = phio(i,j) - Advx * Fluxx - Advy * Fluxy &
                + Visx * Diffx + Visy * Diffy
        enddo
    enddo


end subroutine FTCS_Scheme
