program test_nan
    use ieee_arithmetic
    implicit none
    real :: x
    x = -2
    x = SQRT(x)
    print *, x
    if(ieee_is_nan(x)) then
        print *, 'AD Nan'
    end if
  end program test_nan