MODULE global_variables
    INTEGER:: i,j,i2d,iter,AD_itermax,PPE_itermax,solvetype_ppe,solvetype_AD,iterx,itery
    INTEGER :: nx,ny,nxf,nyf
    INTEGER :: restart, re_time, write_inter, write_flag
    REAL :: errorx, errory, errorppe, errormax
    REAL :: lx, ly
    REAL :: w_AD, w_PPE, tmax, dt, t, Re, mu
    REAL :: cflx, cfly, rx, ry, AAD, BAD, APPE, BPPE
    REAL :: x0, y0, r0, r, vt
    REAL :: start, finish
    REAL, ALLOCATABLE, DIMENSION(:) :: x, y, xf, yf
    REAL, ALLOCATABLE, DIMENSION(:) :: bx, un, ukx, bcx, by, vn, vky, bcy, bp
    REAL, ALLOCATABLE, DIMENSION(:,:) :: u, v, uk, vk, ukp1, vkp1, vor, p, pk, uf, vf,velmag
    REAL, ALLOCATABLE, DIMENSION(:,:) :: dx, dy    
END MODULE global_variables

MODULE boundary_conditions
    REAL, PARAMETER :: u_bc_w = 1, u_bc_e = 0, u_bc_n = 0, u_bc_s = 0
    REAL, PARAMETER :: v_bc_w = 0, v_bc_e = 0, v_bc_n = 0, v_bc_s = 0 
    REAL, PARAMETER :: p_bc_w = 0, p_bc_e = 0, p_bc_n = 0, p_bc_s = 0
END MODULE boundary_conditions

MODULE immersed_boundary
    REAL :: a, b
    REAL :: xcent, ycent
    INTEGER :: ncinside, shape
    REAL, ALLOCATABLE, DIMENSION(:,:) :: iblank_cc, iblank_fcv, iblank_fcu, ghost
END MODULE immersed_boundary

MODULE STATS
    REAL :: force_drag, coeff_drag
    REAL, ALLOCATABLE, DIMENSION(:,:) :: p_coeff_drag
END MODULE

module files 
    integer, parameter :: x_face = 31, y_face = 32, gridfile = 33
end module
