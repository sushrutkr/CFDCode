subroutine writepostproc()
    
    use global_variables

    vor = 0
    DO j=2,ny-1
        DO i = 2,nx-1
            vor(i,j) = ((v(i+1,j) - v(i-1,j))/(2*dx)) - ((u(i,j+1) - u(i,j-1))/(2*dy)) 
        END DO
    END DO

    DO j=1,ny
        DO i = 1,nx
            velmag(i,j) = (u(i,j)**2 + v(i,j)**2)**0.5 
        END DO
    END DO
    
    
    !Writing Files For Post Processing
    open(12, file='data.dat', status='unknown')
    WRITE(12,*) 'TITLE = "Post Processing Tecplot"'
    WRITE(12,*) 'VARIABLES = "X", "Y", "u", "v","Velocity Magnitude", "P", "Vorticity"'
    WRITE(12,*) 'ZONE T="BIG ZONE", I=',nx,', J=',ny,', DATAPACKING=POINT'
    
    DO j=1,ny
        DO i = 1,nx
            write(12,*) x(i), y(j), U(i,j), V(i,j),velmag(i,j), P(i,j), vor(i,j) 
        END DO
    END DO
    close(12)

    open(13, file="data_ufacevel.dat", status='unknown')
    WRITE(13,*) 'TITLE = "Post Processing Tecplot"'
    WRITE(13,*) 'VARIABLES = "X", "Y", "uf"'
    WRITE(13,*) 'ZONE T="BIG ZONE", I=',nx,', J=',ny-2,', DATAPACKING=POINT'
    DO j=1,ny-2
        DO i = 1,nx
            write(13,*) i, j, uf(i,j)
        END DO
    END DO
    ClOSE(13)

    open(14, file="data_vfacevel.dat", status='unknown')
    WRITE(14,*) 'TITLE = "Post Processing Tecplot"'
    WRITE(14,*) 'VARIABLES = "X", "Y", "vf"'
    WRITE(14,*) 'ZONE T="BIG ZONE", I=',nx-2,', J=',ny,', DATAPACKING=POINT'
    DO j=1,ny
        DO i = 1,nx-2
            write(14,*) i, j, vf(i,j)
        END DO
    END DO
    ClOSE(14)

end subroutine