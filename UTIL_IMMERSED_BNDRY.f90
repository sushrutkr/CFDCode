SUBROUTINE read_body()
    use immersed_boundary

    OPEN(3, file='canonical_body_in.dat',status='old')
    READ(3,*)
    READ(3,*)
    READ(3,*) shape
    READ(3,*)
    READ(3,*)
    READ(3,*) a, b 
    READ(3,*)
    READ(3,*)
    READ(3,*) xcent, ycent
    CLOSE(3)
END SUBROUTINE

SUBROUTINE calc_in_cells()
    USE global_variables
    USE immersed_boundary
    
    iblank_cc = 1
    iblank_fcu = 1
    iblank_fcv = 1

    ncinside = 0 
    aoa = (pi/180)*30
    
    SELECT CASE(shape)

    case(0)
        write(11,*) 'No Body Added'
        iblank_fcu = 1
        iblank_fcv = 1
        iblank_cc = 1

    case(1)
        write(11,*) 'Circle added with radius = ', a
        DO j = 1,ny
            DO i = 1,nx
                ! I am not dumb to not use same equation for circle or ellipse it's there to reduce chances of error
                IF (((x(i) - xcent)**2 + (y(j) - ycent)**2) .LE. a**2) THEN
                    iblank_cc(i,j) = 0
                    iblank_fcu(i,j-1) = 0
                    iblank_fcu(i-1,j-1) = 0
                    iblank_fcv(i-1,j) = 0
                    iblank_fcv(i-1,j-1) = 0
                    ncinside = ncinside + 1
                END IF
            END DO
        END DO

    case(2)
        write(11,*) 'Square added with hardcoded parameters'
        DO j = 1,ny
            DO i = 1,nx
                IF (x(i) > 2.5 .AND. x(i)<3.5) THEN
                    IF (y(j) > 2 .AND. y(j) < 3) THEN
                        iblank_cc(i,j) = 0
                        iblank_fcu(i,j-1) = 0
                        iblank_fcu(i-1,j-1) = 0
                        iblank_fcv(i-1,j) = 0
                        iblank_fcv(i-1,j-1) = 0
                        ncinside = ncinside + 1
                    END IF
                END IF
            END DO
        END DO


    case(3)
        write(11,*) 'Ellipse added with major axis length = ', a, ' & minor axis length = ', b, ' at an AoA of ', aoa
        DO j = 1,ny
            DO i = 1,nx
                IF ((((x(i)-xcent)*cos(aoa) - (y(j)-ycent)*sin(aoa))**2/a**2 &
                    +((x(i)-xcent)*sin(aoa) + (y(j)-ycent)*cos(aoa))**2/b**2) .LE. 1) THEN
                    iblank_cc(i,j) = 0
                    iblank_fcu(i,j-1) = 0
                    iblank_fcu(i-1,j-1) = 0
                    iblank_fcv(i-1,j) = 0
                    iblank_fcv(i-1,j-1) = 0
                    ncinside = ncinside + 1
                END IF
            END DO
        END DO

    end select     

    CALL calc_ghost()
END SUBROUTINE

SUBROUTINE calc_ghost()
    USE global_variables
    USE immersed_boundary

    ghost = 0 

    DO j = 1,ny 
        DO i =1,nx
            IF (iblank_cc(i,j) .EQ. 0) THEN 
                IF (iblank_cc(i,j-1) .EQ. 1) THEN
                    ghost(i,j) = 1
                ELSEIF (iblank_cc(i-1,j) .EQ. 1) THEN
                    ghost(i,j) = 1
                ELSEIF (iblank_cc(i,j+1) .EQ. 1) THEN
                    ghost(i,j) = 1
                ELSEIF (iblank_cc(i+1,j) .EQ. 1) THEN
                    ghost(i,j) = 1
                ! ELSE 
                !     ghost(i,j) = 0
                END IF            
            END IF
        END DO 
    END DO
END SUBROUTINE

    
