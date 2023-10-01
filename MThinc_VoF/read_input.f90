!======================================================================================
! sunroutine to read input parameters
!--------------------------------------------------------------------------------------
Subroutine read_input
    use Mod_grid
    implicit none 

    open(10,file='input.in',form='formatted')
    write(*,*) 'Reading paraeters...'
    read(10,*)
    read(10,*) Nx, Ny
    write(*,*) 'Grid number is (', Nx, Ny, '). ' 
    read(10,*)
    read(10,*) CFL 
    read(10,*)
    read(10,*) Start_time, End_Time, delta_t
    read(10,*)
    read(10,*) nprint
    read(10,*)
    read(10,*) iproblem
    close(10)

End Subroutine read_input


