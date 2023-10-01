!==================================================================================
! save the computational result to files 
!----------------------------------------------------------------------------------
Subroutine Output_Result
    use Mod_UAD2d
    use Mod_Grid
    implicit none
    integer :: i, j
    character(len=30) :: outfile, strs

    strs = ".dat"
    write(outfile,'("field", i4.4, a4)') nprint, strs
    outfile = trim(out_path)//"/"//trim(outfile)
    open(8,file=outfile,status='unknown')

    write(8,*) 'VARIABLES = "X" "Y" "Ux" "Uy" "phi_exat" "phi_num" '
    write(8,10) time, Nx, Ny 
    10 format('ZONE solutiontime=', 1PG15.7e2, ', I=', I5, 2X, ', J=', I5, 2X)
    
    do j = js, je 
        do i = is, ie 
            write(8,'(2f10.6,4e15.7)') x(i), y(j), U(i,j), V(i,j), phi_exat(i,j), phi(i,j)
        enddo 
    enddo

    close(Unit = 8) ! close file after writen

End Subroutine Output_Result