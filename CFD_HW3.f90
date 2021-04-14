program CFD_HW3

    implicit none

    ! Variables
    INTEGER ::  i, j, k, i2d, iter,nx, ny, itermax
    !INTEGER, PARAMETER :: nx=100, ny=50 
    REAL :: errormax, errorx, errory, error, tmax, dt, Re
    REAL :: r, rx, ry, cflx, cfly, A, B, C,  t, w, dx, dy
    REAL, PARAMETER :: mu = 0.01, Lx = 2, Ly = 1, Vt = 0.25, x0 =0.5, y0 = 0.5, r0 = 0.5
    REAL, PARAMETER :: n_steps = 2
    REAL, ALLOCATABLE, DIMENSION(:) :: x, bx, un, ukx, bcx, by, vn, vky, bcy
    REAL, ALLOCATABLE, DIMENSION(:) :: y
    REAL, ALLOCATABLE, DIMENSION(:,:) :: U, V, Uk, Vk, Ukp1, Vkp1, vor
    ! REAL, ALLOCATABLE, DIMENSION(:) ::  !nx-2
    REAL, DIMENSION(2) :: err
    
    CALL readdata(nx, ny, dx, dy, w, errormax,tmax, dt, Re,itermax)
    open(11, file='log.dat', status='old')
    
    ! Starting Simulation    
    write(11,*) " CODE BEGUN.........................."
    write(11,*)
    write(11,*) " Simulation Parameters are nx, ny, dx, dy, w, errormax,tmax, dt, Re "
    write(11,*) nx, ny, dx, dy, w, errormax,tmax, dt, Re
    write(11,*)
    write(11,*) "Initialising the domain .........................."
    write(11,*)
    
    allocate(x(nx), y(ny))
    allocate(U(nx,ny)) 
    allocate(Uk(nx,ny)) 
    allocate(Ukp1(nx,ny)) 
    allocate(V(nx,ny)) 
    allocate(Vk(nx,ny)) 
    allocate(Vkp1(nx,ny))
    allocate(vor(nx,ny))
    allocate(bx(nx-2), by(nx-2), un(nx-2), vn(nx-2), ukx(nx-2), vky(nx-2), bcx(nx-2), bcy(nx-2))
    
    ! Initialising the variables and arrays
    rx = (mu*dt)/(dx**2)
    ry = (mu*dt)/(dy**2)
    
    A = 1 + 2*rx
    B = -rx
    C = -ry
    
    errormax = 1E-10 
    
    ! Initialising the Domain with a vortex at (0.5,0.5) and Applying Boundary Conditions
    x(1) = 0
    y(1) = 0
    
    DO i = 2,nx
        x(i) = x(i-1) + dx
    END DO
    
    DO j = 2,ny
        y(j) = y(j-1) + dy
    END DO  
    
    U(1,:) = 1
    U(nx,:) = 1
    U(:,1)  = 1
    U(:,ny) = 1
    
    V(1,:) = 0
    V(nx,:) = 0
    V(:,1)  = 0
    V(:,ny) = 0
    
    DO j=2,ny-1
        DO i = 2,nx-1
            r = SQRT((x(i) - x0)**2 + (y(j) - y0)**2)
            U(i,j) = 1 - Vt*(y(j) - y0)*EXP((1-(r/r0)**2)/2)
            V(i,j) = Vt*(x(i) - x0)*EXP((1-(r/r0)**2)/2)
        END DO
    END DO
    
    Uk = 1
    Vk = 1
    Ukp1 = U
    Vkp1 = V
    
    w = 1
    
    t = 0
    
    DO WHILE (t<tmax)
        errorx = 1
        errory = 1
        error = 1
        Uk = U
        Vk = V
        
        iter = 0
        DO k=1,itermax !WHILE (error>errormax)
            DO j=2,ny-1
                IF (errorx>errormax) THEN
                    ! x-Advection Difussion Iterations 
                    un = 0
                    ukx = 0
                    bx = 0 
                    bcx = 0
    
                    DO i=1,nx-2
                        i2d = i+1
                        cflx = (U(i2d,j) * dt)/dx
                        cfly = (V(i2d,j) * dt)/dx
                        un(i) = U(i2d,j) + (-cflx/2)*(U(i2d+1,j) - U(i2d-1,j)) + (-cfly/2)*(U(i2d,j+1) - U(i2d,j-1))
                        ukx(i) = -C*(Ukp1(i2d,j+1) - 2*Ukp1(i2d,j) + Ukp1(i2d,j-1))
                    END DO
    
                    bcx(1) = -B*U(1,j) 
                    bcx(nx-2) = -B*U(nx, j)
                    bx = un + ukx + bcx
                    
                    CALL TDMA(nx, bx, A, B, B)
                    
                    Ukp1(2:nx-1,j) = bx(:)
                    Ukp1(:,j) = Uk(:,j)*(1-w) + w*Ukp1(:,j)
                END IF
                
                IF (errory>errormax) THEN
                    ! y-Advection Difussion Iterations
                    vn = 0
                    vky = 0
                    by = 0 
                    bcy = 0
    
                    DO i=1,nx-2
                        i2d = i+1
                        cflx = (U(i2d,j) * dt)/dx
                        cfly = (V(i2d,j) * dt)/dx
                        vn(i) = V(i2d,j) + (-cflx/2)*(V(i2d+1,j) - V(i2d-1,j)) + (-cfly/2)*(V(i2d,j+1) - V(i2d,j-1))
                        vky(i) = -C*(Vk(i2d,j+1) - 2*Vk(i2d,j) + Vkp1(i2d,j-1))
                    END DO
    
                    bcy(1) = -B*V(1,j) 
                    bcy(nx-2) = -B*V(nx, j)
                    by = vn + vky + bcy
                    
                    CALL TDMA(nx, by, A, B, B)
                    
                    Vkp1(2:nx-1,j) = by(:)
                    Vkp1(:,j) = Vk(:,j)*(1-w) + w*Vkp1(:,j)                    
                END IF    
            END DO
            errorx = MAXVAL(ABS(Ukp1) - ABS(Uk))
            Uk(:,:) = Ukp1(:,:)
            errory = MAXVAL(ABS(Vkp1) - ABS(Vk))
            Vk(:,:) = Vkp1(:,:)
            err(1) = errorx
            err(2) = errory
            error = MAXVAL(err) 
            iter = iter + 1
        END DO
        U(:,:) = Ukp1(:,:)
        V(:,:) = Vkp1(:,:)
        write(11,*) 'At t = ', t,' Iteration = ', iter, ' Error = ',error
        t = t+dt
    END DO
    close(11)
    
    vor = 0
    DO j=2,ny-1
        DO i = 2,nx-1
            vor(i,j) = ((V(i+1,j) - V(i-1,j))/(2*dx)) - ((U(i,j+1) - U(i,j-1))/(2*dy)) 
        END DO
    END DO
    
    !Writing Files For Post Processing
    open(12, file='Vortex Convect.dat', status='old')
    WRITE(12,*) 'TITLE = "Vortex Convect"'
    WRITE(12,*) 'VARIABLES = "X", "Y", "u", "v", "Vorticity"'
    WRITE(12,*) 'ZONE T="BIG ZONE", I=',nx,', J=',ny,', DATAPACKING=POINT'
    DO j=1,ny
        DO i = 1,nx
            write(12,*) x(i), y(j), U(i,j), V(i,j), vor(i,j)
        END DO
    END DO
    
    close(12)


    deallocate(x)
    deallocate(y)
    deallocate(U)
    deallocate(Ukp1)
    deallocate(Uk) 
    deallocate(Vkp1)
    deallocate(Vk)
    deallocate(V)
    deallocate(vor)
    deallocate(bx, by, un, vn, ukx, vky, bcx, bcy)

end program CFD_HW3




subroutine readdata(nx, ny, dx, dy, w, errormax,tmax, dt, Re, itermax)
    integer:: nx,ny, itermax
    real :: dx, dy, w, errormax,tmax, dt, Re
    ! Reading Input Data
    open(2, file='input.dat', status='old')
    read(2,*) 
    read(2,*)
    read(2,*)
    read(2,*)
    read(2,*) nx, ny
    read(2,*)
    read(2,*)
    read(2,*) dx, dy
    read(2,*)
    read(2,*)
    read(2,*)
    read(2,*)
    read(2,*) w, itermax
    read(2,*)
    read(2,*)
    read(2,*)
    read(2,*)
    read(2,*) errormax, tmax, dt, Re
    close(2)
end subroutine