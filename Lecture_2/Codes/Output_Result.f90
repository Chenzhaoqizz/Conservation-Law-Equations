!==================================================================================
! save the computational result to files 
!----------------------------------------------------------------------------------
Subroutine Output_Result
    use Mod_Grid
    use Mod_Poisson2d
    implicit none
    integer :: i, j, k, nfile
    real(8) :: phi_exat, err
    character(len=30) :: outfile, strs

    nfile = 001
    err = 0.d0
    strs = ".dat"
    write(outfile,'("field", i3.3, a4)') nfile, strs
    outfile = trim(out_path)//"/"//trim(outfile)

    open(8,file=outfile,status='unknown')

    write(8,'(a,f12.5,a)') 'TITLE=" Time= ',0.d0, ' "'
    write(8,*) ' VARIABLES = "X" "Y" "phi_exat" "phi_num" '
    write(8,10) Nx, Ny
    10 format('ZONE T="EF" I=',i5, ',J=',i5, ',F=POINT')

    ! transformation between 2d array and 1d array
    ! do i = is, ie
    !     do j = js, je
    !         k = i + Nx * (j - 1)
    !         phi(i,j) = phi1d(k)
    !     enddo
    ! enddo
    
    do j = js, je 
        do i = is, ie 
            k = j + Ny * (i - 1)
            phi_exat = cos( x(i) - y(j) ) * exp( x(i) - y(j) )
            err = err + abs( phi1d(k) - phi_exat ) / abs( phi_exat + 1e-32 )
            write(8,'(2f10.5,2e15.7)') x(i), y(j), phi_exat, phi1d(k)
        enddo 
    enddo

    close(Unit = 8) ! close file after writen

    Err = Err / dble( (Nx-0)*(Ny-0) ) 
    write(*,20) Err
    20 Format('Final error = ', 3X, 1e15.8) 

End Subroutine Output_Result