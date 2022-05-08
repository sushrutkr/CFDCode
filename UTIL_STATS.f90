SUBROUTINE calc_drag()
    use global_variables
    use immersed_boundary
    use stats

    p_coeff_drag = 0
    do j=1,ny
        do i=1,nx
            p_coeff_drag(i,j) = iblank_cc(i,j)*(ghost(i+1,j)*p(i,j) - ghost(i-1,j)*p(i,j))
        end do
    end do

    force_drag = 0
    do j=1,ny
        do i=1,nx
            force_drag = force_drag + (dy(i,j)*p_coeff_drag(i,j))
            !print *, force_drag
        end do
    end do

    coeff_drag = 2*force_drag
END SUBROUTINE