!--------------------------------------------------------------------------------------
!           create mesh and set initial conditions here
!======================================================================================
subroutine Initialize
    use flow_var
    implicit none

    integer :: i, j

    is = 1; ie = Nx; 
    js = 1; je = Ny;

    time = Start_time
    istep = 0

    ! generate mesh 
    select case(imesh)
    case(1)
        call uniform_mesh
    case(2)
        ! call tanh_mesh
    case(3) 
        ! call sin_mesh
    end select

    AE = 0.d0; AW = 0.d0; 
    AN = 0.d0; AS = 0.d0;
    AP = 0.d0;

    do j = js + 1, je - 1 
        do i = is + 1, ie - 1
            AW(i,j) = 2.d0 / ( (x(i+1)-x(i-1)) * (x(i)-x(i-1)) )
            AE(i,j) = 2.d0 / ( (x(i+1)-x(i-1)) * (x(i+1)-x(i)) )
            AS(i,j) = 2.d0 / ( (y(j+1)-y(j-1)) * (y(j)-y(j-1)) )
            AN(i,j) = 2.d0 / ( (y(j+1)-y(j-1)) * (y(j+1)-y(j)) )
            AP(i,j)  = - ( AW(i,j) + AE(i,j) + AS(i,j) + AN(i,j) )
        enddo
    enddo

    ! print out mesh
    open(Unit=8, File="Grid.plt", Status='unknown')
        write(8,'(a,f12.5,a)') 'TITLE=" Time= ',0.d0, ' "'
        write(8,*) ' VARIABLES = "X" "Y" '
        write(8,10) Nx, Ny
        10 format('ZONE T="EF" I=',1I5, ',J=',1I5, ',F=POINT')
        do j = js, je 
            do i = is, ie 
                write(8,'(2f10.5,2e15.7)') x(i), y(j)
            enddo 
        enddo
    close(Unit=8)

    ! Set initial condition
    u=0.d0; v=0.d0; vort=0.d0; sf=0.d0; p=0.d0;

    ! Set boundary condition of velocity at top boundary, u = const
    u(:, je) = uTop
    vort(:, je) = 2.d0*( sf(:, je) - sf(:,je-1) ) / ( y(je) - y(je-1) )**2 & 
        - 2.d0 * uTop / ( y(je) - y(je-1) ) 

end subroutine Initialize
!======================================================================================
subroutine uniform_mesh
    use flow_var
    implicit none
    integer :: i, j

    do i = is, ie 
        delta_x(i) = (xRight-xLeft)/real(ie-is)
        x(i) = delta_x(i)*real(i-is) + xLeft
    enddo

    do j = js, je
        delta_y(j) = (yTop-yBottom)/real(je-js)
        y(j) = delta_y(j)*real(j-js) + yBottom
    enddo

end subroutine uniform_mesh
!======================================================================================
