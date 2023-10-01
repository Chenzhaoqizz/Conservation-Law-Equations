!==================================================================================
! save the computational result to files 
!----------------------------------------------------------------------------------
Subroutine Output_Result
    use Mod_UAD1d
    implicit none
    integer :: j, nfile
    real(8) :: x
    character(len=30) :: outfile, strs

    nfile = nprint
    strs = ".dat"
    write(outfile,'("field", i3.3, a4)') nfile, strs
    outfile = trim(out_path)//"/"//trim(outfile)

    open(8,file=outfile,status='unknown')
    do j = is, ie
        x = j*hx
        write(8,10) x, f(j)
    enddo
    10 format(2x,2e13.5)

End Subroutine Output_Result