!==================================================================================
! subroutine that uses BTCS scheme to solve 
! the unsteady advection diffusion equation
!----------------------------------------------------------------------------------
Subroutine BTCS_Scheme
    use Mod_UAD1d
    implicit none
    integer :: M 
    real(8) :: beta, lamb
    real(8), dimension(Nx-2) :: Dia, Lower, Up, RM, x

    fo = f; 
    M = Nx - 2;

    lamb = 0.5d0*U*dt/hx;
    beta = D*dt/(hx*hx);

    ! Assemble the coefficient matrix
    Dia(1:M) = 1.d0 + 2.d0*beta;

    Lower(1) = 0.d0;
    Lower(2:M) = -(beta+lamb);

    Up(1:M-1) = lamb - beta;
    Up(M) = 0.d0;

    ! Assemble the Right hand side(RHS.)
    RM(1) = f(2) + f(1)*(lamb+beta);
    RM(2:M-1) = f(3:Nx-2);
    RM(M) = f(Nx-1) + f(Nx)*(beta-lamb);

    call TDMA(M,Dia,Lower,Up,RM,x)
    f(2:Nx-1) = x(1:M);
    
    call SetupBC

    return
End Subroutine BTCS_Scheme
!----------------------------------------------------------------------------------
