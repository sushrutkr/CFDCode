SUBROUTINE readdata()
    USE global_variables
    ! READing Input Data
    OPEN(2, file='input.dat', status='old')
    READ(2,*) 
    READ(2,*)
    READ(2,*) restart, re_time
    READ(2,*)
    READ(2,*)
    READ(2,*)
    READ(2,*) nx, ny
    READ(2,*)
    READ(2,*)
    READ(2,*) lx, ly
    READ(2,*)
    READ(2,*)
    READ(2,*)
    READ(2,*)
    READ(2,*) w_AD, w_PPE, AD_itermax, PPE_itermax, solvetype_AD, solvetype_ppe
    READ(2,*)
    READ(2,*)
    READ(2,*)
    READ(2,*)
    READ(2,*) errormax, tmax, dt, Re, mu
    READ(2,*)
    READ(2,*)
    READ(2,*)
    READ(2,*) write_inter
    CLOSE(2)

    ! nx = nx + 2
    ! ny = ny + 2
    ! dx = lx/(nx) !!!!!!!!!!! TO BE COMMENTED LATER ON
    ! dy = ly/(ny)
    nxf = nx + 1
    nyf = ny + 1
    nx = nx + 2
    ny = ny + 2
END SUBROUTINE

SUBROUTINE domain_init()
    USE global_variables
    use files

    real :: dummy

    
    ! dxc = lx/(nx-2)
    ! dyc = ly/(ny-2)
    
    ! Reading Face Coordinates
    open(x_face, file='xgrid.dat', status='old')
    do i = 1, nxf
        read(x_face,*) dummy, xf(i)
    enddo
    close(x_face)

    open(y_face, file='ygrid.dat', status='old')
    do i = 1, nyf
        read(y_face,*) dummy, yf(i)
    enddo
    close(y_face)

    ! Setting up the cell centers
    DO i = 2,nx-1
        x(i) = (xf(i-1)+xf(i))/2.0 
    END DO
    
    DO j = 2,ny-1
        y(j) = (yf(j-1)+yf(j))/2.0
    END DO


    x(1) = -x(2)
    y(1) = -y(2)
    x(nx) = xf(nxf) + (xf(nxf) - x(nx-1)) 
    y(ny) = yf(nyf) + (yf(nyf) - y(ny-1))
    

    do j=2,ny-1
        do i=2,nx-1
            dx(i,j) = xf(i) - xf(i-1)
            dy(i,j) = yf(j) - yf(j-1)
        enddo
    enddo
    dx(1,:) = dx(2,:)
    dx(nx,:) = dx(nx-1,:)
    dx(:,1) = dx(:,2)
    dx(:,ny) = dx(:,ny-1)

    dy(1,:) = dy(2,:)
    dy(nx,:) = dy(nx-1,:)
    dy(:,1) = dy(:,2)
    dy(:,ny) = dy(:,ny-1)

    open(gridfile, file='grid.dat', status='unknown')
    write(*,*) 'Writing Grid Data'
    WRITE(gridfile,*) 'TITLE = "Post Processing Tecplot"'
    WRITE(gridfile,*) 'VARIABLES = "X", "Y", "dxi", "dyi"'
    WRITE(gridfile,*) 'ZONE T="BIG ZONE", I=',nx,', J=',ny,', DATAPACKING=POINT'

    do j=1,ny
        do i=1,nx
            WRITE(gridfile,*) i,j,dx(i,j),dy(i,j)
        enddo
    enddo
    close(gridfile)

END SUBROUTINE

SUBROUTINE flow_init()

    USE global_variables 
    USE immersed_boundary
    USE boundary_conditions

    character(len=50) :: fname, no, ext  
    real :: dump 
    
    IF (restart .EQ. 0) THEN
        u = 0
        v = 0
        uf = 0
        vf = 0
        p = 0
        t = 0
        CALL set_dirichlet_bc()
        CALL set_neumann_bc()

    ELSE 
        t = re_time
        write_flag = 1000*re_time
        ext = '.dat'
        fname = 'data.'
        WRITE(no, "(I7.7)") re_time
        fname = TRIM(ADJUSTL(fname))//no
        fname = TRIM(ADJUSTL(fname))//TRIM(ADJUSTL(ext))
        open(3, file=fname, status='unknown')
        READ(3,*) 
        READ(3,*) 
        READ(3,*) 
        
        DO j=1,ny
            DO i = 1,nx
                READ(3,*) x(i), y(j), u(i,j), v(i,j), p(i,j), dump, dump, dump
            END DO
        END DO
        close(3)
    END IF

    CALL set_SSM_bc()
    
    DO j=1,ny-2
        DO i=2,nx-2
            uf(i,j) = iblank_fcu(i,j)*(u(i+1,j+1) + u(i,j+1))/2
        END DO
    END DO

    DO j = 2,ny-2
        DO i = 1,nx-2
            vf(i,j) = iblank_fcv(i,j)*(v(i+1,j+1) + v(i+1,j))/2
        END DO 
    END DO

    ! Comment This out when appying top and bottom wall bc
    ! uf(nx,:) = iblank_fcu(nx,:)*(2*u(nx,2:ny-1) - uf(nx-1,:))
    ! vf(:,ny) = iblank_fcv(:,ny)*(2*v(2:nx-1,ny) - vf(:,ny-1))

END SUBROUTINE

