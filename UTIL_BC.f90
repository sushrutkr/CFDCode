MODULE boundary_conditions
    REAL, PARAMETER :: u_bc_w = 1, u_bc_e = 0, u_bc_n = 0, u_bc_s = 0
    REAL, PARAMETER :: v_bc_w = 0, v_bc_e = 0, v_bc_n = 0, v_bc_s = 0 
    REAL, PARAMETER :: p_bc_w = 0, p_bc_e = 0, p_bc_n = 0, p_bc_s = 0
END MODULE boundary_conditions

SUBROUTINE set_dirichlet_bc()
    USE global_variables
    USE boundary_conditions
    USE immersed_boundary

    !Dirichlet
    !for cell center velocities
    u(1,:) = -u(2,:) + u_bc_w*2
    ! u(nx,:) = -u(nx-1,:) + u_bc_e*2
    ! u(:,1) = -u(:,2) + u_bc_s*2
    ! u(:,ny) = -u(:,ny-1) + u_bc_n*2

    v(1,:) = -v(2,:) + v_bc_w*2
    ! v(nx,:) = -v(nx-1,:) + v_bc_e*2
    v(:,1) = -v(:,2) + v_bc_s*2
    v(:,ny) = -v(:,ny-1) + v_bc_n*2

    !for face velocities
    ! uf(nx-1,:) = u_bc_e
    uf(1,:) = u_bc_w
    vf(:,1) = v_bc_s
    vf(:,ny-1) = v_bc_n

    uf(nx,:) = 2*u(nx,2:ny-1) - uf(nx-1,:)
    vf(:,ny) = 2*v(2:nx-1,ny) - vf(:,ny-1)
END SUBROUTINE

SUBROUTINE set_neumann_bc()
    USE global_variables
    USE boundary_conditions
    USE immersed_boundary

    !Neumann
    p(nx,:) = p(nx-1,:)
    p(1,:) = p(2,:)
    p(:,ny) = p(:,ny-1)
    p(:,1) = p(:,2)

    u(nx,:) = u(nx-1,:)
    u(:,1) = u(:,2)
    u(:,ny) = u(:,ny-1)

    v(nx,:) = v(nx-1,:)
    ! v(:,1) = v(:,2)
    ! v(:,ny) = v(:,ny-1)
    uf(nx,:) = 2*u(nx,2:ny-1) - uf(nx-1,:)
    vf(:,ny) = 2*v(2:nx-1,ny) - vf(:,ny-1)
END SUBROUTINE

SUBROUTINE set_pressure_BC()
    USE global_variables
    USE boundary_conditions
    USE immersed_boundary


    !To satify dp/dn = 0 at every iteration    
    p(nx,:) = p(nx-1,:)
    p(1,:) = p(2,:)
    p(:,ny) = p(:,ny-1)
    p(:,1) = p(:,2)

END SUBROUTINE

SUBROUTINE set_velocity_BC()
    USE global_variables
    USE boundary_conditions
    USE immersed_boundary

    u(nx,:) = u(nx-1,:)
    u(:,1) = u(:,2)
    u(:,ny) = u(:,ny-1)

    v(nx,:) = v(nx-1,:)
    ! v(:,1) = v(:,2)
    ! v(:,ny) = v(:,ny-1)
END SUBROUTINE

SUBROUTINE set_SSM_bc()
    use global_variables
    use immersed_boundary

    u = u*iblank_cc
    v = v*iblank_cc
    ! p = p*iblank_cc

    uf = uf*iblank_fcu
    vf = vf*iblank_fcv
END SUBROUTINE