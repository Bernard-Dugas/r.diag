program test16
! Checks that the current FORTRAN compiler
! support quad-math operation natively or not ...
real(16) :: var=123456789012345678901234567890.0_16
print *,'var =',var,' = 123456789012345678901234567890.0_16'
stop
end program test16
