!==================================================================================
!        RK3 scheme to solve the unsteady advection diffusion equation
!----------------------------------------------------------------------------------
subroutine RK3_Scheme
    use Mod_UAD2d
    use Mod_Grid
    implicit none

    integer :: i, j
    real(8) :: Advx, Advy, Fluxx, Fluxy, AdvectionFlux 
    real(8) :: Visx, Visy, Diffx, Diffy, DiffusionFlux 
    real(8), dimension(Nx, Ny) :: phio, phi1

    phio = phi
    Advx = 0.5d0*dt/hx; Advy = 0.5d0*dt/hy;
    Visx = D*dt/(hx)**2; Visy = D*dt/(hy)**2;

    ! RK3 step 1
    do i = is+1, ie-1 
        do j = js+1, je-1 
            Fluxx = AdvectionFlux( u(i+1,j), phio(i+1,j), u(i-1,j), phio(i-1,j) )
            Fluxy = AdvectionFlux( v(i,j+1), phio(i,j+1), v(i,j-1), phio(i,j-1) )
            Diffx = DiffusionFlux( phio(i+1,j), phio(i,j), phio(i-1,j) )
            Diffy = DiffusionFlux( phio(i,j+1), phio(i,j), phio(i,j-1) )

            phi(i,j) = phio(i,j) & 
                - Advx * Fluxx - Advy * Fluxy & 
                + Visx * Diffx + Visy * Diffy 
        enddo
    enddo

    phi1 = phi

    ! RK3 step 2
    do i = is+1, ie-1 
        do j = js+1, je-1 
            Fluxx = AdvectionFlux( u(i+1,j), phi1(i+1,j), u(i-1,j), phi1(i-1,j) )
            Fluxy = AdvectionFlux( v(i,j+1), phi1(i,j+1), v(i,j-1), phi1(i,j-1) )
            Diffx = DiffusionFlux( phi1(i+1,j), phi1(i,j), phi1(i-1,j) )
            Diffy = DiffusionFlux( phi1(i,j+1), phi1(i,j), phi1(i,j-1) )
            phi(i,j) = 0.75d0 * phio(i,j) + 0.25d0 * phi1(i,j) + 0.25d0 * ( & 
                - Advx * Fluxx - Advy * Fluxy & 
                + Visx * Diffx + Visy * Diffy )
        enddo
    enddo

    phi1 = phi

    ! RK3 step 3
    do i = is+1, ie-1 
        do j = js+1, je-1 
            Fluxx = AdvectionFlux( u(i+1,j), phi1(i+1,j), u(i-1,j), phi1(i-1,j) )
            Fluxy = AdvectionFlux( v(i,j+1), phi1(i,j+1), v(i,j-1), phi1(i,j-1) )
            Diffx = DiffusionFlux( phi1(i+1,j), phi1(i,j), phi1(i-1,j) )
            Diffy = DiffusionFlux( phi1(i,j+1), phi1(i,j), phi1(i,j-1) )
            phi(i,j) = (1.d0/3.d0) * phio(i,j) + (2.d0/3.d0) * phi1(i,j) + (2.d0/3.d0) * ( &
                - Advx * Fluxx - Advy * Fluxy & 
                + Visx * Diffx + Visy * Diffy ) 
        enddo
    enddo

    return
end subroutine RK3_Scheme

!==================================================================================
!        RK2 scheme to solve the unsteady advection diffusion equation
!----------------------------------------------------------------------------------
subroutine RK2_Scheme
    use Mod_UAD2d
    use Mod_Grid
    implicit none

    integer :: i, j
    real(8) :: Advx, Advy, Fluxx, Fluxy, AdvectionFlux 
    real(8) :: Visx, Visy, Diffx, Diffy, DiffusionFlux 
    real(8), dimension(Nx, Ny) :: phio, phi1

    phio = phi
    phi1 = phi
    Advx = 0.5d0*dt/hx; Advy = 0.5d0*dt/hy;
    Visx = D*dt/(hx)**2; Visy = D*dt/(hy)**2;

    ! RK2 step 1
    do i = is+1, ie-1 
        do j = js+1, je-1 
            Fluxx = AdvectionFlux( u(i+1,j), phio(i+1,j), u(i-1,j), phio(i-1,j) )
            Fluxy = AdvectionFlux( v(i,j+1), phio(i,j+1), v(i,j-1), phio(i,j-1) )
            Diffx = DiffusionFlux( phio(i+1,j), phio(i,j), phio(i-1,j) )
            Diffy = DiffusionFlux( phio(i,j+1), phio(i,j), phio(i,j-1) )

            phi(i,j) = phio(i,j) & 
                - Advx * Fluxx - Advy * Fluxy & 
                + Visx * Diffx + Visy * Diffy 
        enddo
    enddo

    ! RK2 step 2
    do i = is+1, ie-1 
        do j = js+1, je-1 
            Fluxx = AdvectionFlux( u(i+1,j), phi(i+1,j), u(i-1,j), phi(i-1,j) )
            Fluxy = AdvectionFlux( v(i,j+1), phi(i,j+1), v(i,j-1), phi(i,j-1) )
            Diffx = DiffusionFlux( phi(i+1,j), phi(i,j), phi(i-1,j) )
            Diffy = DiffusionFlux( phi(i,j+1), phi(i,j), phi(i,j-1) )
            phi1(i,j) = phi(i,j) & 
                - Advx * Fluxx - Advy * Fluxy & 
                + Visx * Diffx + Visy * Diffy 
        enddo
    enddo
    phi = 0.5d0 *( phio + phi1 )

    return 
end subroutine RK2_Scheme

!==================================================================================
! Advection terms
!----------------------------------------------------------------------------------
real(8) function AdvectionFlux(UR, phiR, UL, phiL)
    implicit none
    real(8), intent(in) :: UR, phiR, UL, phiL

    AdvectionFlux = UR*phiR - UL*phiL
end function AdvectionFlux

!==================================================================================
! Diffusion terms
!----------------------------------------------------------------------------------
real(8) function DiffusionFlux(phiR, phiC, phiL)
    implicit none
    real(8), intent(in) :: phiR, phiC, phiL

    DiffusionFlux = phiR - 2.d0*phiC + phiL
end function DiffusionFlux
!==================================================================================