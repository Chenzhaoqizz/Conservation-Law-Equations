!==================================================================================
! save the computational result to files 
!----------------------------------------------------------------------------------
subroutine Error_Func
    use Mod_UAD2d
    use Mod_Grid
    implicit none
    integer :: i, j
    real(8) :: t

    Error = 0.d0
    if ( TestCase == 0 ) then 
        t = -2.d0*D*pi*pi*time
        do j = js, je 
            do i = is, ie 
                phi_exat(i,j) = cos(pi*x(i)) * cos(pi*y(j)) * exp(t)
                Error = Error + ( phi(i,j) - phi_exat(i,j) )**2
            enddo 
        enddo
    elseif ( TestCase == 1 ) then 
        t = 4.d0*time + 1.d0
        do j = js, je 
            do i = is, ie 
                phi_exat(i,j) = exp( ( (x(i)-u(i,j)*time-1.d0)**2 & 
                    + (y(j)-v(i,j)*time-1.d0)**2 ) / ( -D*t ) ) / t
                Error = Error + ( phi(i,j) - phi_exat(i,j) )**2
            enddo 
        enddo
    endif

    Error = sqrt(Error) / ( (Nx-1)*(Ny-1) )
    ! write(*,20) Error
    ! 20 Format('Error = ', 3X, 1e15.8) 

    return 
end subroutine Error_Func