!======================================================================================
! sunroutine to output result
!--------------------------------------------------------------------------------------
Subroutine output(prefix)
    use Mod_grid
    use Mod_flow
    implicit none
    integer :: i, j
    character(len=3) :: prefix
    character(len=30) :: outfile, strs

    strs = ".dat"
    write(outfile,'(a3, I4.4, a4)') prefix, nprint, strs
    outfile = trim("Result")//"/"//trim(outfile)
    ! outfile = trim(outfile)
    open(Unit = 8, file=outfile, status='unknown')

    write(8,*) 'Variables = "X" "Y" "VoF" "Ux" "Uy" '
    write(8,10) time, Nx, Ny 
    10 format('ZONE solutiontime=', 1PG15.7e2, ', I=', I5, 2X, ', J=', I5, 2X)
    
    do j = js, je 
        do i = is, ie
            write(8,'(2f10.6,4e15.7)') (i-1)*hxi, (j-1)*hyi, phi(i,j), U(i,j), V(i,j)
        enddo 
    enddo

    close(Unit = 8) ! close file after writen

End Subroutine output
