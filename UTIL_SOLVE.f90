subroutine ADSolver()
    use global_variables
    use immersed_boundary
    use ieee_arithmetic
    use files

    REAL, DIMENSION(2) :: err
    REAL, DIMENSION(nx,ny) :: sx, sy 
    real, dimension(nx-2) :: d
    real, dimension(nx-3) :: l,upp
    REAL, DIMENSION(nx,ny) :: coeff_ad, coeff_dy2_p1, coeff_dx2_p1, coeff_dy2_m1, coeff_dx2_m1

    !CALL set_SSM_bc()
    error = 1
    errorx = 1
    errory = 1
    uk(:,:) = u(:,:)
    vk(:,:) = v(:,:)
    ukp1(:,:) = u(:,:)
    vkp1(:,:) = v(:,:) 
    sx = 0
    sy = 0
    DO j = 2,ny-1
        DO i = 2,nx-1
            ! Method A
            ! sx(i,j) = (u(i,j) - (dt/(2*dx))*(u(i+1,j)*uf(i+1,j-1) - u(i-1,j)*uf(i-1,j-1)) &
            !                   - (dt/(2*dy))*(u(i,j+1)*vf(i-1,j+1) - u(i,j-1)*vf(i-1,j-1)))

            ! sy(i,j) = (v(i,j) - (dt/(2*dx))*(v(i+1,j)*uf(i+1,j-1) - v(i-1,j)*uf(i-1,j-1)) &
            !                   - (dt/(2*dy))*(v(i,j+1)*vf(i-1,j+1) - v(i,j-1)*vf(i-1,j-1)))

            ! Method B
            sx(i,j) = (u(i,j) - (dt/dx(i,j))*(((1/(dx(i,j) + dx(i+1,j)))*(u(i+1,j)*dx(i,j) + u(i,j)*dx(i+1,j))*uf(i,j-1)) &
                                             -((1/(dx(i,j) + dx(i-1,j)))*(u(i,j)*dx(i-1,j) + u(i-1,j)*dx(i,j))*uf(i-1,j-1))) &
                              - (dt/dy(i,j))*(((1/(dy(i,j) + dy(i,j+1)))*(u(i,j+1)*dy(i,j) + u(i,j)*dy(i,j+1))*vf(i-1,j)) &
                                             -((1/(dy(i,j) + dy(i,j-1)))*(u(i,j)*dy(i,j-1) + u(i,j-1)*dy(i,j))*vf(i-1,j-1))))

            sy(i,j) = (v(i,j) - (dt/dx(i,j))*(((1/(dx(i,j) + dx(i+1,j)))*(v(i+1,j)*dx(i,j) + v(i,j)*dx(i+1,j))*uf(i,j-1)) &
                                             -((1/(dx(i,j) + dx(i-1,j)))*(v(i,j)*dx(i-1,j) + v(i-1,j)*dx(i,j))*uf(i-1,j-1))) &
                              - (dt/dy(i,j))*(((1/(dy(i,j) + dy(i,j+1)))*(v(i,j+1)*dy(i,j) + v(i,j)*dy(i,j+1))*vf(i-1,j)) &
                                             -((1/(dy(i,j) + dy(i,j-1)))*(v(i,j)*dy(i,j-1) + v(i,j-1)*dy(i,j))*vf(i-1,j-1))))
            ! sy(i,j) = (v(i,j) - (dt/(dx))*(0.5*(v(i+1,j) + v(i,j))*uf(i,j-1) - 0.5*(v(i,j) + v(i-1,j))*uf(i-1,j-1)) &
            !                   - (dt/(dy))*(0.5*(v(i,j+1) + v(i,j))*vf(i-1,j) - 0.5*(v(i,j) + v(i,j-1))*vf(i-1,j-1))) 

            !Method C
            ! sx(i,j) = u(i,j) - (dt/(2*dx))*(u(i+1,j)*(uf(i+1,j-1) + uf(i,j-1))*0.5 - u(i-1,j)*(uf(i-1,j-1) + uf(i,j-1))*0.5) &
            !                  - (dt/(2*dy))*(u(i,j+1)*(vf(i-1,j+1) + vf(i-1,j))*0.5 - u(i,j-1)*(vf(i-1,j-1) + vf(i,j-1))*0.5)

            ! sy(i,j) = v(i,j) - (dt/(2*dx))*(v(i+1,j)*(uf(i+1,j-1) + uf(i,j-1))*0.5 - v(i-1,j)*(uf(i-1,j-1) + uf(i,j-1))*0.5) &
            !                  - (dt/(2*dy))*(v(i,j+1)*(vf(i-1,j+1) + vf(i-1,j))*0.5 - v(i,j-1)*(vf(i-1,j-1) + vf(i,j-1))*0.5)    
                             
            !Ji's Methods
            ! sx(i,j) = u(i,j) - (dt/(2*dx))*(u(i+1,j) - u(i-1,j))*(0.5*(uf(i+1,j-1) + uf(i-1,j-1))) &
            !                  - (dt/(2*dy))*(u(i,j+1) - u(i,j-1))*(0.5*(vf(i-1,j+1) + vf(i-1,j-1))) 

            ! sy(i,j) = v(i,j) - (dt/(2*dx))*(v(i+1,j) - v(i-1,j))*(0.5*(uf(i+1,j-1) + uf(i-1,j-1))) &
            !                  - (dt/(2*dy))*(v(i,j+1) - v(i,j-1))*(0.5*(vf(i-1,j+1) + vf(i-1,j-1)))       
        END DO
    END DO
    ! open(13, file='source_term.dat', status='unknown')
    ! WRITE(13,*) 'TITLE = "Post Processing Tecplot"'
    ! WRITE(13,*) 'VARIABLES = "X", "Y", "U", "V", "Sx","Sy"' !"Vorticity", "iblank", "Ghost", "Velocity Magnitude"'
    ! WRITE(13,*) 'ZONE T="BIG ZONE", I=',nx,', J=',ny,', DATAPACKING=POINT'
    
    ! DO j=1,ny
    !     DO i = 1,nx
    !         WRITE(13,*) i, j, U(i,j), V(i,j), sx(i,j), sx(i,j)!, vor(i,j), iblank_cc(i,j), ghost(i,j), SQRT(u(i,j)**2 + v(i,j)**2)
    !     END DO
    ! END DO
    ! close(13)

    ! coeff = 1 + ((2*dt)/(Re*(dx**2))) + ((2*dt)/(Re*(dy**2)))
    coeff_ad = 1
    coeff_dx2_m1 = 1
    coeff_dx2_p1 = 1 
    coeff_dy2_m1 = 1 
    coeff_dy2_p1 = 1
    
    do j = 2,ny-1
        do i=2,nx-1
            coeff_ad(i,j) = 1 + (dt/Re)*((2/(dx(i,j)*(dx(i,j) + dx(i+1,j)))) + (2/(dx(i,j)*(dx(i,j) + dx(i-1,j))))) &
                              + (dt/Re)*((2/(dy(i,j)*(dy(i,j) + dy(i,j+1)))) + (2/(dy(i,j)*(dy(i,j) + dy(i,j-1)))))

            coeff_dx2_m1(i,j) = (dt/Re)*(2/(dx(i,j)*(dx(i,j) + dx(i-1,j))))
            coeff_dx2_p1(i,j) = (dt/Re)*(2/(dx(i,j)*(dx(i,j) + dx(i+1,j))))
            coeff_dy2_m1(i,j) = (dt/Re)*(2/(dy(i,j)*(dy(i,j) + dy(i,j-1))))
            coeff_dy2_p1(i,j) = (dt/Re)*(2/(dy(i,j)*(dy(i,j) + dy(i,j+1))))
        enddo
    enddo

    ! open(gridfile, file='coeff_ad.dat', status='unknown')
    ! WRITE(gridfile,*) 'TITLE = "Post Processing Tecplot"'
    ! WRITE(gridfile,*) 'VARIABLES = "X", "Y", "coeff", "coeff_dx2_m1", "coeff_dx2_p1", "coeff_dy2_m1", "coeff_dy2_p1"'
    ! WRITE(gridfile,*) 'ZONE T="BIG ZONE", I=',nx,', J=',ny,', DATAPACKING=POINT'

    ! do j=1,ny
    !     do i=1,nx
    !         WRITE(gridfile,*) i,j,coeff_ad(i,j), coeff_dx2_m1(i,j), coeff_dx2_p1(i,j), coeff_dy2_m1(i,j), coeff_dy2_p1(i,j)
    !     enddo
    ! enddo
    ! close(gridfile)
    
    
    iter = 0
    ! CALL set_velocity_BC()
    SELECT CASE(solvetype_AD)
        CASE(1)
            DO WHILE(error>errormax)
                !CALL set_velocity_BC()
                DO j = 2,ny-1
                    DO i =2,nx-1
                        ukp1(i,j) = sx(i,j) + (coeff_dx2_p1(i,j)*ukp1(i+1,j) + coeff_dx2_m1(i,j)*ukp1(i-1,j)) &
                                            + (coeff_dy2_p1(i,j)*ukp1(i,j+1) + coeff_dy2_m1(i,j)*ukp1(i,j-1))
                        ukp1(i,j) = (iblank_cc(i,j)*ukp1(i,j))/coeff_ad(i,j) !
                        ukp1(i,j) = (1-w_AD)*uk(i,j) + w_AD*ukp1(i,j)

                        vkp1(i,j) = sy(i,j) + (coeff_dx2_p1(i,j)*vkp1(i+1,j) + coeff_dx2_m1(i,j)*vkp1(i-1,j)) &
                                            + (coeff_dy2_p1(i,j)*vkp1(i,j+1) + coeff_dy2_m1(i,j)*vkp1(i,j-1))
                        vkp1(i,j) = (iblank_cc(i,j)*vkp1(i,j))/coeff_ad(i,j) !
                        vkp1(i,j) = (1-w_AD)*vk(i,j) + w_AD*vkp1(i,j)            
                    END DO
                END DO
                iter = iter + 1
                ! errorx = 0
                ! errory = 0
                !CALL calc_residual(ukp1,uk,errorx,nx,ny) 
                errorx = MAXVAL(ABS(ukp1 - uk))
                uk(:,:) = ukp1(:,:)
                errory = MAXVAL(ABS(vkp1- vk))
                !CALL calc_residual(vkp1,vk,errory,nx,ny)
                vk(:,:) = vkp1(:,:)
                err(1) = errorx
                err(2) = errory
                error = MAXVAL(err) 
                ! if(ieee_is_nan(maxval(ukp1))) then
                !     print *, 'AD Nan', t, iter
                ! end if
                ! if(ieee_is_nan(maxval(vkp1))) then
                !     print *, 'AD Nan', t, iter
                ! end if
            END DO
            !CALL set_velocity_BC()
            !print *, solvetype_AD, 'AD', error, 'error', t

        CASE(2)
            ! DO WHILE(error>errormax)
            !     DO j=2,ny-1
            !         ! x-Advection Difussion Iterations 
            !         un = 0
            !         ukx = 0
            !         bx = 0 
            !         bcx = 0

            !         DO i=1,nx-2
            !             i2d = i+1
            !             un(i) = sx(i2d,j)
            !             ukx(i) = (dt/(Re*(dy**2)))*(ukp1(i2d,j+1) - 2*ukp1(i2d,j) + ukp1(i2d,j-1))
            !         END DO

            !         bcx(1) = -BAD*u(1,j) 
            !         bcx(nx-2) = -BAD*u(nx,j)
            !         bx = un + ukx + bcx
                    
            !         d(:) = AAD 
            !         l(:) = BAD 
            !         upp(:) = BAD
            !         CALL TDMA_vary(nx, bx, d, l, upp)
                    
            !         ukp1(2:nx-1,j) = bx(:)
            !         ukp1(:,j) = uk(:,j)*(1-w_AD) + w_AD*ukp1(:,j)
                    
            !         ! y-Advection Difussion Iterations
            !         vn = 0
            !         vky = 0
            !         by = 0 
            !         bcy = 0

            !         DO i=1,nx-2
            !             i2d = i+1
            !             vn(i) = sy(i2d,j)
            !             vky(i) = (dt/(Re*(dy**2)))*(vkp1(i2d,j+1) - 2*vkp1(i2d,j) + vkp1(i2d,j-1))
            !         END DO

            !         bcy(1) = -BAD*v(1,j) 
            !         bcy(nx-2) = -BAD*v(nx, j)
            !         by = vn + vky + bcy

            !         d(:) = AAD 
            !         l(:) = BAD 
            !         upp(:) = BAD
            !         CALL TDMA_vary(nx, by, d, l, upp)
                    
            !         vkp1(2:nx-1,j) = by(:)
            !         vkp1(:,j) = vk(:,j)*(1-w_AD) + w_AD*vkp1(:,j)                     
            !     END DO
            !     errorx = MAXVAL(ABS(ukp1) - ABS(uk))
            !     uk(:,:) = ukp1(:,:)
            !     errory = MAXVAL(ABS(vkp1) - ABS(vk))
            !     vk(:,:) = vkp1(:,:)
            !     err(1) = errorx
            !     err(2) = errory
            !     error = MAXVAL(err) 
            ! END DO    
    END SELECT 
    
    u(:,:) = ukp1(:,:)
    v(:,:) = vkp1(:,:)

    ! open(13, file="AFTER_AD_SOLVE.dat", status='unknown')
    ! WRITE(13,*) 'TITLE = "Post Processing Tecplot"'
    ! WRITE(13,*) 'VARIABLES = "X", "Y", "U", "V", "P",' !"Vorticity", "iblank", "Ghost", "Velocity Magnitude"'
    ! WRITE(13,*) 'ZONE T="BIG ZONE", I=',nx,', J=',ny,', DATAPACKING=POINT'
    
    ! DO j=1,ny
    !     DO i = 1,nx
    !         WRITE(13,*) i, j, U(i,j), V(i,j), P(i,j)!, vor(i,j), iblank_cc(i,j), ghost(i,j), SQRT(u(i,j)**2 + v(i,j)**2)
    !     END DO
    ! END DO
    ! close(13)
end subroutine

subroutine PPESolver()
    use global_variables
    use immersed_boundary
    use ieee_arithmetic
    use files
    REAL, DIMENSION(nx,ny) :: ps, pkp1
    REAL :: factor, f1, f2
    REAL, DIMENSION(nx,ny) :: coeff_ppe, coeff_dy2_p1, coeff_dx2_p1, coeff_dy2_m1, coeff_dx2_m1
    
    ps = 0
    pkp1 = p
    errorppe = 1
    pk(:,:) = p(:,:)

    ! Calculating face velocities 
    DO j=1,ny-2
        DO i=2,nx-1
            uf(i,j) = iblank_fcu(i,j)*((1/(dx(i,j+1) + dx(i+1,j+1)))*(u(i+1,j+1)*dx(i,j+1) + u(i,j+1)*dx(i+1,j+1)))
        END DO
    END DO

    DO j = 2,ny-2
        DO i = 1,nx-2
            vf(i,j) = iblank_fcv(i,j)*((1/(dy(i+1,j+1) + dy(i+1,j)))*(v(i+1,j+1)*dy(i+1,j) + v(i+1,j)*dy(i+1,j+1)))
        END DO 
    END DO

    CALL mass_correct()

    DO j=2,ny-1
        DO i=2,nx-1
            ps(i,j) = ((1/dt)*(((uf(i,j-1) - uf(i-1,j-1))/dx(i,j)) + ((vf(i-1,j) - vf(i-1,j-1))/dy(i,j))))
        END DO 
    END DO

    ! uf(nx,:) = 2*u(nx,2:ny-1) - uf(nx-1,:)
    ! vf(:,ny) = 2*v(2:nx-1,ny) - vf(:,ny-1)


    iter = 0
    ! coeff_abs = -((2/(dx**2)) + (2/(dy**2)))
    ! dx2 = (1/(dx**2))
    ! dy2 = (1/(dy**2))
    coeff_ppe = 1
    coeff_dx2_m1 = 1
    coeff_dx2_p1 = 1 
    coeff_dy2_m1 = 1 
    coeff_dy2_p1 = 1

    do j = 2,ny-1
        do i=2,nx-1
            coeff_ppe(i,j) = -1*(((2/(dx(i,j)*(dx(i,j) + dx(i+1,j)))) + (2/(dx(i,j)*(dx(i,j) + dx(i-1,j)))))&
                                +((2/(dy(i,j)*(dy(i,j) + dy(i,j+1)))) + (2/(dy(i,j)*(dy(i,j) + dy(i,j-1))))))

            coeff_dx2_m1(i,j) = (2/(dx(i,j)*(dx(i,j) + dx(i-1,j))))
            coeff_dx2_p1(i,j) = (2/(dx(i,j)*(dx(i,j) + dx(i+1,j))))
            coeff_dy2_m1(i,j) = (2/(dy(i,j)*(dy(i,j) + dy(i,j-1))))
            coeff_dy2_p1(i,j) = (2/(dy(i,j)*(dy(i,j) + dy(i,j+1))))
        enddo
    enddo

    ! open(gridfile, file='coeff_ppe.dat', status='unknown')
    ! WRITE(gridfile,*) 'TITLE = "Post Processing Tecplot"'
    ! WRITE(gridfile,*) 'VARIABLES = "X", "Y", "coeff", "coeff_dx2_m1", "coeff_dx2_p1", "coeff_dy2_m1", "coeff_dy2_p1"'
    ! WRITE(gridfile,*) 'ZONE T="BIG ZONE", I=',nx,', J=',ny,', DATAPACKING=POINT'

    ! do j=1,ny
    !     do i=1,nx
    !         WRITE(gridfile,*) i, j,coeff_ppe(i,j), coeff_dx2_m1(i,j), coeff_dx2_p1(i,j), coeff_dy2_m1(i,j), coeff_dy2_p1(i,j)
    !     enddo
    ! enddo
    ! close(gridfile)

    CALL set_pressure_BC()
    DO WHILE(errorppe>errormax)
        DO j=2,ny-1
            DO i=2,nx-1
                ! coeff = coeff_abs + dx2*ghost(i+1,j) + dx2*ghost(i-1,j) &
                !                   + dy2*ghost(i,j+1) + dy2*ghost(i,j-1)

                ! p(i,j) = ps(i,j)
                ! p(i,j) = p(i,j) - ((iblank_cc(i+1,j)*p(i+1,j) + iblank_cc(i-1,j)*p(i-1,j))/(dx**2)) &
                !                 - ((iblank_cc(i,j+1)*p(i,j+1) + iblank_cc(i,j-1)*p(i,j-1))/(dy**2))
                ! p(i,j) = iblank_cc(i,j)*p(i,j)/coeff

                ! p(i,j) = ps(i,j)
                p(i,j) = ps(i,j) - (p(i+1,j)*coeff_dx2_p1(i,j) + p(i-1,j)*coeff_dx2_m1(i,j)) &
                                 - (p(i,j+1)*coeff_dy2_p1(i,j) + p(i,j-1)*coeff_dy2_m1(i,j))

                p(i,j) = iblank_cc(i,j)*p(i,j)/coeff_ppe(i,j)
                p(i,j) = (1-w_PPE)*pk(i,j) + w_PPE*p(i,j)

                pkp1(i,j) = p(i,j)
                factor = ghost(i,j)*(iblank_cc(i,j-1) + iblank_cc(i,j+1) + iblank_cc(i-1,j) + iblank_cc(i+1,j)) &
                            + (1-ghost(i,j))
                f1 = (1-ghost(i,j))*p(i,j)
                f2 = ghost(i,j)*((iblank_cc(i,j-1)*p(i,j-1) + iblank_cc(i,j+1)*p(i,j+1) &
                                                            + iblank_cc(i-1,j)*p(i-1,j) &
                                                            + iblank_cc(i+1,j)*p(i+1,j))/factor)
                p(i,j) =  f1 + f2 

                if(ieee_is_nan(p(i,j))) then
                    print *, 'PPE Iter Nan', t, i, j, ghost(i,j), p(i,j), pkp1(i,j), &
                                         f1,f2, factor, p(i,j-1), p(i,j+1), p(i+1,j), p(i-1,j) !iblank_cc(i,j-1)*p(i,j-1),&
                                         !iblank_cc(i,j+1)*p(i,j+1), iblank_cc(i-1,j)*p(i-1,j),&
                                         !iblank_cc(i+1,j)*p(i+1,j)
                end if
            END DO
        END DO
        !errorppe = 0
        errorppe = MAXVAL(ABS(p - pk))
        !CALL calc_residual(p,pk,errorppe,nx,ny)
        pk(:,:) = p(:,:)
        iter = iter + 1
        ! if(ieee_is_nan(maxval(p))) then
        !     print *, 'PPE Nan', t, iter
        ! end if
        !print *, 'PPE', errorppe, iter
        !CALL set_pressure_BC()
    END DO
    CALL set_pressure_BC()
    !print *, 'PPE', t, iter, coeff, coeff_abs

    ! open(13, file="AFTER_PPE_SOLVE.dat", status='unknown')
    ! WRITE(13,*) 'TITLE = "Post Processing Tecplot"'
    ! WRITE(13,*) 'VARIABLES = "X", "Y", "U", "V", "P",' !"Vorticity", "iblank", "Ghost", "Velocity Magnitude"'
    ! WRITE(13,*) 'ZONE T="BIG ZONE", I=',nx,', J=',ny,', DATAPACKING=POINT'
    
    ! DO j=1,ny
    !     DO i = 1,nx
    !         WRITE(13,*) i, j, U(i,j), V(i,j), P(i,j)!, vor(i,j), iblank_cc(i,j), ghost(i,j), SQRT(u(i,j)**2 + v(i,j)**2)
    !     END DO
    ! END DO
    ! close(13)
end subroutine

subroutine vel_correct()
    use global_variables
    use immersed_boundary
    !Cell Center Velocity Correction
    DO j=2,ny-1
        DO i=2,nx-1
            u(i,j) = iblank_cc(i,j)*(u(i,j) - ((2*dt)/(2*dx(i,j)+dx(i+1,j)+dx(i-1,j)))*(p(i+1,j) - p(i-1,j)))
            v(i,j) = iblank_cc(i,j)*(v(i,j) - ((2*dt)/(2*dy(i,j)+dy(i,j+1)+dy(i,j-1)))*(p(i,j+1) - p(i,j-1)))
        END DO 
    END DO

    !Face Center Velocity Correction
    DO j=1,ny-2
        DO i=2,nx-1
            uf(i,j) = iblank_fcu(i,j)*(uf(i,j) - ((2*dt)/(dx(i+1,j+1)+dx(i,j+1)))*(p(i+1,j+1) - p(i,j+1)))
        END DO
    END DO

    DO j = 2,ny-2
        DO i = 1,nx-2
            vf(i,j) = iblank_fcv(i,j)*(vf(i,j) - ((2*dt)/(dy(i+1,j+1)+dy(i+1,j)))*(p(i+1,j+1) - p(i+1,j)))
        END DO 
    END DO

    ! open(13, file="AFTER_VEL_CORRECT.dat", status='unknown')
    ! WRITE(13,*) 'TITLE = "Post Processing Tecplot"'
    ! WRITE(13,*) 'VARIABLES = "X", "Y", "U", "V", "P",' !"Vorticity", "iblank", "Ghost", "Velocity Magnitude"'
    ! WRITE(13,*) 'ZONE T="BIG ZONE", I=',nx,', J=',ny,', DATAPACKING=POINT'
    
    ! DO j=1,ny
    !     DO i = 1,nx
    !         WRITE(13,*) i, j, U(i,j), V(i,j), P(i,j)!, vor(i,j), iblank_cc(i,j), ghost(i,j), SQRT(u(i,j)**2 + v(i,j)**2)
    !     END DO
    ! END DO
    ! close(13)
end subroutine

SUBROUTINE mass_correct()
    USE global_variables
    USE immersed_boundary
    REAL :: left_sum, right_sum, mass_error

    left_sum = 0
    left_sum = 0
    DO j = 1,ny-2
        left_sum  = left_sum  + dy(2,j+1)*uf(1,j)
        right_sum = right_sum + dy(2,j+1)*uf(nx-1,j) 
    END DO

    mass_error = left_sum - right_sum
    mass_error = mass_error/ly

    DO j = 1,ny-2
        uf(nx-1,j) = uf(nx-1,j) + mass_error
    END DO
    
END SUBROUTINE

SUBROUTINE calc_residual(mat_old, mat_new, res, nx, ny)
    integer, intent(in) :: nx, ny
    real ,dimension(nx,ny), intent(in) :: mat_old, mat_new 
    real, intent(inout) :: res 
    real :: sum_new, sum_old, total

    sum_new = 0
    sum_old = 0  

    total = (nx-2)*(ny-2)
    do j = 2,ny-1
        do i = 2,nx-1
            sum_new = sum_new + mat_new(i,j)
            sum_old = sum_old + mat_old(i,j)
        end do
    end do

    sum_new = abs(sum_new/total )
    sum_old = abs(sum_old/total)

    res = abs(sum_new - sum_old)
END SUBROUTINE

SUBROUTINE calc_poisson_residual()
    use global_variables
    use immersed_boundary
    real, dimension(nx,ny) :: errormat
    real :: errorsum
    errormat(:,:) = 0

    do j = 2,ny-1
        do i = 2,nx-1
            ! errormat(i,j) = ps(i,j) - ((p(i+1,j) -2*p(i,j) + p(i-1,j))/(dx**2)) + ((p(i,j+1) -2*p(i,j) + p(i,j-1))/(dy**2))

            ! errormat(i,j) = ((p(i+1,j) -2*p(i,j) + p(i-1,j))/(dx**2)) + ((p(i,j+1) -2*p(i,j) + p(i,j-1))/(dy**2))
            ! errormat(i,j) = errormat(i,j) - ((1/dt)*(((uf(i,j-1) - uf(i-1,j-1))/dx) + ((vf(i-1,j) - vf(i-1,j-1))/dy)))
        end do 
    end do

    errorsum = 0

    do j = 2,ny-1
        do i = 2,nx-1
            errorsum = errorsum + errormat(i,j) 
        end do
    end do
    
    errorppe = (errorsum)/((nx-2)*(nx-2))
END SUBROUTINE