subroutine messys( string )
#  if defined (__HOS_AIX__)
   USE XLFUTILITY_EXTNAME
   character * (*) string
   character * 512 ligne
   integer ier
   ier = ierrno( )
   write(ligne,'("SYSERROR no. ",I5)') ier
#  else
#  if defined (__INTEL_COMPILER_UPDATE)
   USE IFCORE, only: GERROR
#  endif
   character * (*) string
   character * 512 ligne
   call gerror( ligne )
#  endif
   write(6,'(A)') string
   write(6,'(A)') ligne

   return

end subroutine messys

function CHKENDI() result(outval)

   use ISO_C_BINDING  

   implicit none

   ! Determines the endian-ess of a computer
   ! Version re-coded on Nov 6, 2014 (after M. Valin)

   integer :: outval

   type(C_PTR) :: ptr

   integer (kind=2), dimension(:), pointer :: I2
   integer (kind=4),                target :: I4

   I4 = 1 ; ptr = c_loc(I4) ; call c_f_pointer( ptr,I2,[2] )

   outval = 1 - I2(1) ! = 0 (Little) or 1 (Big) endian
   ! if ( I2(1) == 1 ) outval = 0 ! Little endian machine
   ! if ( I2(1) == 0 ) outval = 1 ! Big endian machine

   return

end function CHKENDI

function ME32O64() result(outval)

   implicit none

   integer :: outval
 
   ! Tell calling program whether the current environment is
   ! configured with default 32-bit integers (2) or with
   ! default 64-bit integers (1)

   ! Original version by A.J. Stacey -MARCH 23,1992-
   ! Re-coded (following M. Valin) by B. Dugas, Nov 2014

   integer, dimension(2) :: INUM

   outval = -1

   if (loc( INUM(2) ) - loc( INUM(1) ) == 4) then
      outval = 2 ! 32-bits integers
   else if (loc( INUM(2) ) - loc( INUM(1) ) == 8) then
      outval = 1 ! 64-bits integers
   end if
 
   return

end function ME32O64

function INTEFLT() result(outval)

   implicit none

   integer :: outval
 
   ! This routine checks if integer word size is the same
   ! as the real word size. This is subtly different from
   ! the output of "ME32O64".

   ! Original version by A.J. Stacey - March 1992
   ! Re-coded (following M. Valin) by B. Dugas, Nov 2014

   real,    dimension(2) :: ANUM
   integer, dimension(2) :: INUM

   if ((loc( INUM(2) ) - loc( INUM(1) )) == &
       (loc( ANUM(2) ) - loc( ANUM(1) ))) then
      outval = 1 ! default integer and real sizes are identical
   else
      outval = 2 ! default integer and real sizes are different
   end if
 
   return

end function INTEFLT
