!==============================================================================
!        define variables here 
!------------------------------------------------------------------------------
subroutine read_input
    use flow_var
    implicit none

    open(10,file='input.in',form='formatted')
    write(*,*) 'Reading paraeters...'
    read(10,*)
    read(10,*) Nx, Ny
    write(*,100) Nx, Ny
    100 format('Grid number in x direction is ', 1I5, ', grid number in y direction is ', 1I5)
    read(10,*)
    read(10,*) CFL, Epsilon, MaxIter
    read(10,*)
    read(10,*) xLeft, xRight, yBottom, yTop
    read(10,*)
    read(10,*) Start_time, End_Time, delta_t
    read(10,*)
    read(10,*) uTop, Re
    write(*,110) uTop, Re 
    110 format('Velocity at top boundary is ', 1e13.6, ', Reynold number is ', 1e13.6)
    read(10,*)
    read(10,*) isave
    write(*,*) 'iSave = ', iSave
    read(10,*)
    read(10,*) imesh
    write(*,*) 'imesh = ', imesh
    read(10,*)
    read(10,*) iPossion_Equation
    write(*,*) 'iPossion_Equation = ', iPossion_Equation
    read(10,*)
    read(10,*) iflux_splitting
    read(10,*)
    read(10,*) ischeme_inv
    close(unit=10)

end subroutine read_input

