SUBROUTINE calc_drag()
    use global_variables
    use immersed_boundary
    use stats

    p_coeff_drag = 0
    p_coeff_lift = 0
    force_drag = 0
    force_lift = 0

    do j=2,ny-1
        do i=2,nx-1
            p_coeff_drag(i,j) = iblank_cc(i,j)*(ghost(i+1,j)*p(i,j) - ghost(i-1,j)*p(i,j))*dy(i,j)
            force_drag = force_drag + p_coeff_drag(i,j)
            p_coeff_lift(i,j) = iblank_cc(i,j)*(ghost(i,j+1)*p(i,j) - ghost(i,j-1)*p(i,j))*dx(i,j)
            force_lift = force_lift + p_coeff_lift(i,j)
        end do
    end do

    ! force_drag = 0
    ! force_lift = 0
    ! do j=1,ny
    !     do i=1,nx
    !         force_drag = force_drag + (dy(i,j)*p_coeff_drag(i,j))
    !         force_lift = force_lift + (dx(i,j)*p_coeff_lift(i,j))
    !         !print *, force_drag
    !     end do
    ! end do

    coeff_drag = 2*force_drag
    coeff_lift = 2*force_lift
END SUBROUTINE