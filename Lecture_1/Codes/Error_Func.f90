!==================================================================================
!     subroutine that calculates the error between 
!     analytical solution and numerical solution 
!----------------------------------------------------------------------------------
Subroutine Error_Func
    use Mod_UAD1d
    implicit none
    integer :: j
    real(8) :: x, c1, c2, error
    real(8), dimension(Nx) :: f_ana
    
    c1 = sqrt(4.d0*time + 1.d0)
    c2 = D*c1**2
    error = 0.d0
    
    do j = is+1, ie-1
        x = hx*(j-1)
        f_ana(j) = exp( -(x-1.d0-U*time)**2/c2 ) / c1 
        error = error + ( f(j) - f_ana(j) )**2
    enddo

    error = hx*sqrt(error)
	
	do j = is+1, ie-1
		write(*,66) j, f_ana(j)
        66 	format("Analytic solution (",1x,1I6,") = ",2x,1e13.5)
	enddo

    write(*,*) "Error = ", error
	
    return
End Subroutine Error_Func
	