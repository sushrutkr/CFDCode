subroutine writepostproc()
    
    use global_variables
    use boundary_conditions
    use immersed_boundary
    character(len=50) :: fname, no, ext

    !modifying domain data
    x(1) = 0 
    x(nx) = lx
    y(1) = 0
    y(ny) = ly

    !modifying the velocity data
    !Bottom Wall
    j = 1
    DO i = 1,nx
        ! u(i,j) = u_bc_s
        v(i,j) = v_bc_s
    END DO
    u(:,1) = (u(:,1) + u(:,2))/2

    !top wall data
    j = ny
    DO i=1,nx
        ! u(i,j) = u_bc_n
        v(i,j) = v_bc_n
    END DO
    u(:,ny) = (u(:,ny) + u(:,ny-1))/2
    
    !left wall data
    i=1
    DO j=2,ny-1
        u(i,j) = u_bc_w
        v(i,j) = v_bc_w
    END DO

    !right wall data
    ! i=nx
    ! DO j=2,ny-1
    !     u(i,j) = u_bc_e
    !     v(i,j) = v_bc_e
    ! END DO
    u(nx,:) = (u(nx,:) + u(nx-1,:))/2
    v(nx,:) = (v(nx,:) + v(nx-1,:))/2


    vor = 0
    ! DO j=2,ny-1
    !     DO i = 2,nx-1
    !         vor(i,j) = ((v(i+1,j) - v(i-1,j))/(2*dx)) - ((u(i,j+1) - u(i,j-1))/(2*dy)) 
    !     END DO
    ! END DO

    DO j=1,ny
        DO i = 1,nx
            velmag(i,j) = (u(i,j)**2 + v(i,j)**2)**0.5 
        END DO
    END DO
        
    !Writing Files For Post Processing
    ext = '.dat'
    fname = 'Data/data.'
    write(no, "(I7.7)") iter/1000
    fname = TRIM(ADJUSTL(fname))//no
    fname = TRIM(ADJUSTL(fname))//TRIM(ADJUSTL(ext))
    open(12, file='data.dat', status='unknown')
    WRITE(12,*) 'TITLE = "Post Processing Tecplot"'
    WRITE(12,*) 'VARIABLES = "X", "Y", "U", "V","P","Velocity Magnitude",  "Vorticity", "iblank", "Ghost"'
    WRITE(12,*) 'ZONE T="BIG ZONE", I=',nx,', J=',ny,', DATAPACKING=POINT'
    
    DO j=1,ny
        DO i = 1,nx
            WRITE(12,*) x(i), y(j), U(i,j), V(i,j),P(i,j), SQRT(u(i,j)**2 + v(i,j)**2) , vor(i,j), iblank_cc(i,j), ghost(i,j)
        END DO
    END DO
    close(12)
end subroutine

SUBROUTINE data_write()
    use global_variables
    use immersed_boundary
    character(len=50) :: fname, no, ext
    vor = 0
    ! DO j=2,ny-1
    !     DO i = 2,nx-1
    !         vor(i,j) = ((v(i+1,j) - v(i-1,j))/(2*dx)) - ((u(i,j+1) - u(i,j-1))/(2*dy)) 
    !     END DO
    ! END DO

    ext = '.dat'
    fname = 'Data/data.'
    write(no, "(I7.7)") write_flag/1000
    fname = TRIM(ADJUSTL(fname))//no
    fname = TRIM(ADJUSTL(fname))//TRIM(ADJUSTL(ext))

    open(13, file=fname, status='unknown')
    WRITE(13,*) 'TITLE = "Post Processing Tecplot"'
    WRITE(13,*) 'VARIABLES = "X", "Y", "U", "V", "P",' !"Vorticity", "iblank", "Ghost", "Velocity Magnitude"'
    WRITE(13,*) 'ZONE T="BIG ZONE", I=',nx,', J=',ny,', DATAPACKING=POINT'
    
    DO j=1,ny
        DO i = 1,nx
            WRITE(13,*) i, j, U(i,j), V(i,j), P(i,j)!, vor(i,j), iblank_cc(i,j), ghost(i,j), SQRT(u(i,j)**2 + v(i,j)**2)
        END DO
    END DO
    close(13)

END SUBROUTINE