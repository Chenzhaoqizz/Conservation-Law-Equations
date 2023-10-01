!----------------------------------------------------------------------------------
! save the computational result to files 
!==================================================================================
Subroutine output
    use flow_var
    implicit none
    integer :: i, j
    character(len=30) :: outfile, strs

    strs = ".dat"
    write(outfile,'("field", I4.4, a4)') nprint, strs
    ! outfile = trim(out_path)//"/"//trim(outfile)
    outfile = trim("out")//"/"//trim(outfile)
    open(Unit = 8, file=outfile, status='unknown')

    write(8,*) 'Variables = "X" "Y" "Ux" "Uy" "Stream function" "Vorticity" '
    write(8,10) time, Nx, Ny 
    10 format('ZONE solutiontime=', 1PG15.7e2, ', I=', I5, 2X, ', J=', I5, 2X)

    do i = is+1, ie-1
        do j = js+1, je-1
            u(i,j) = ( sf(i,j+1) - sf(i,j-1) ) / ( y(j+1) - y(j-1) )
            v(i,j) = - ( sf(i+1,j) - sf(i-1,j) ) / ( x(i+1) - x(i-1) )
        enddo
    enddo
    
    do j = js, je 
        do i = is, ie
            write(8,'(2f10.6,4e15.7)') x(i), y(j), U(i,j), V(i,j), sf(i,j), vort(i,j)
        enddo 
    enddo

    close(Unit = 8) ! close file after writen

End Subroutine output