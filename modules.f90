module global_variables
    integer:: i,j,i2d,iter,ad_itermax,ppe_itermax,solvetype_ppe,solvetype_ad,iterx,itery
    integer :: nx,ny,nxf,nyf
    integer :: restart, re_time, write_inter, write_flag
    real :: errorx, errory, errorppe, errormax
    real :: lx, ly
    real :: w_ad, w_ppe, tmax, dt, t, re, mu
    real :: cflx, cfly, rx, ry, aad, bad, appe, bppe
    real :: x0, y0, r0, r, vt
    real :: start, finish
    real, parameter :: pi = 3.141592653589793
    real, allocatable, dimension(:) :: x, y, xf, yf
    real, allocatable, dimension(:) :: bx, un, ukx, bcx, by, vn, vky, bcy, bp
    real, allocatable, dimension(:,:) :: u, v, uk, vk, ukp1, vkp1, vor, p, pk, uf, vf,velmag
    real, allocatable, dimension(:,:) :: dx, dy    
end module global_variables

module boundary_conditions
    real, parameter :: u_bc_w = 1, u_bc_e = 0, u_bc_n = 0, u_bc_s = 0
    real, parameter :: v_bc_w = 0, v_bc_e = 0, v_bc_n = 0, v_bc_s = 0 
    real, parameter :: p_bc_w = 0, p_bc_e = 0, p_bc_n = 0, p_bc_s = 0
end module boundary_conditions

module immersed_boundary
    real :: a, b
    real :: xcent, ycent, aoa
    integer :: ncinside, shape
    real, allocatable, dimension(:,:) :: iblank_cc, iblank_fcv, iblank_fcu, ghost
end module immersed_boundary

module stats
    real :: force_drag, coeff_drag, force_lift, coeff_lift
    integer :: nprobes
    real, allocatable, dimension(:,:) :: probe_loc
    integer, allocatable, dimension(:) :: probe_index_x, probe_index_y
    real, allocatable, dimension(:,:) :: p_coeff_drag, p_coeff_lift
end module

module files 
    integer, parameter :: x_face = 31, y_face = 32, gridfile = 33
    integer, parameter :: probe_input_file = 41, probe_output_file = 42
end module
