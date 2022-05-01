! My CFD code 
! Compile - ifort modules.f90 UTIL_PRE_SIM.f90 UTIL_POST_PROC.f90 TDMA.f90 UTIL_SOLVE.f90 UTIL_IMMERSED_BNDRY.f90 UTIL_BC.f90 UTIL_STATS.f90 MAIN.f90
! Read File Flag
!       2 - read input, 3 - readrestart, 11 - log, 12 - data.dat, 13 -midwrite, 14 - drag vs time
program CFDCode
    USE global_variables
    USE immersed_boundary
    USE stats

    CALL cpu_time(start)
    CALL readdata()
    CALL read_body()

    open(11, file='log.dat', status='unknown')
    write(11,*) " CODE BEGUN.........................."
    write(11,*)
    write(11,*) " Simulation Parameters are nx, ny, dx, dy, w, errormax,tmax, dt, Re "
    write(11,*) nx, ny, dx, dy, w_PPE, errormax,tmax, dt, Re
    write(11,*)
    write(11,*) "Initialising the domain .........................."
    write(11,*)  
    allocate(iblank_cc(nx,ny), ghost(nx,ny), iblank_fcu(nx,ny-2), iblank_fcv(nx-2,ny))
    allocate(x(nx), y(ny))
    allocate(u(nx,ny), uk(nx,ny), ukp1(nx,ny)) 
    allocate(v(nx,ny), vk(nx,ny), vkp1(nx,ny)) 
    allocate(p(nx,ny), pk(nx,ny), uf(nx,ny-2), vf(nx-2,ny), bp(nx-2))
    allocate(vor(nx,ny), velmag(nx,ny))
    allocate(bx(nx-2), by(nx-2), un(nx-2), vn(nx-2), ukx(nx-2), vky(nx-2), bcx(nx-2), bcy(nx-2))
    allocate(p_coeff_drag(nx,ny))

    open(14, file="drag_history.dat", status='unknown')
    write(14,*) " Drag around the body"
    write(14,*)
    write(14,*) "Time, Pressure Drag, Form Drag"
    write(14,*) 
    
    ! Defining Coefficient for [A]x=b   
    AAD = 1 + ((2*dt)/(Re*(dx**2)))
    BAD = -dt/(Re*(dx**2))

    APPE = -2/(dx**2)
    BPPE = 1/(dx**2)

    CALL domain_init()   
    CALL calc_in_cells()
    CALL flow_init()
    write(11,*) "Initiating Calculations .........................."
    write(11,*)
    write(11,*) "time,  Error_x-AD,  Error_y-AD,  Error_PPE"
        
    IF (restart .EQ. 0) THEN
        t = 0
    ELSE
        t = re_time
    END IF
    DO WHILE (t<tmax)
        CALL set_dirichlet_bc()
        CALL set_neumann_bc()
        CALL ADSolver()
        CALL PPESolver()
        CALL vel_correct() 

        write(11,*) t, errorx, errory, errorppe

        write_flag = write_flag + 1
        IF (write_flag .EQ. 275) THEN
            CALL data_write()
            
        END IF

        IF (mod(write_flag,write_inter) .EQ. 0) THEN
            CALL data_write()
        END IF

        CALL calc_drag

        write(14,*) t, force_drag, coeff_drag
        t = t+dt
    END DO

    write(11,*)
    write(11,*) "Calculations Done  .............................."
    write(11,*)
    write(11,*) "Writing file for Post-Processing  ................"
    write(11,*)
    
    CALL writepostproc()
    
    deallocate(iblank_cc, ghost, iblank_fcu, iblank_fcv)
    deallocate(x,y)
    deallocate(u,uk,ukp1)
    deallocate(v,vk,vkp1)
    deallocate(vor,velmag)
    deallocate(p, pk, bp, vf, uf)
    deallocate(bx, by, un, vn, ukx, vky, bcx, bcy)
    deallocate(p_coeff_drag)
    call cpu_time(finish)

    write(11,*) "Total Simulation Time : ", (finish-start)/60 , "mins"
    close(11)

end program CFDCode
