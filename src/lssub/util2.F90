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

