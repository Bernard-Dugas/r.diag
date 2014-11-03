!
!---------------------------------- LICENCE BEGIN -------------------------------
! R.DIAG - Diagnostic tool kit for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This code is free software; you can redistribute it and/or modify it 
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------
!
      subroutine handle_err2(status,name)

      implicit none

!!!   Version R.DIAG par B. Dugas

      include 'netcdf.inc'
      character*(*) name
      integer status

      if(status .ne. nf_noerr) then
        write(6,1111)nf_strerror(status),'status=', status
        print *,'detected error - stopped in ',trim(name)
        call xit('handle_err',-1)
      endif
!-----------------------------------------------------------------------
 1111 format(/,a80,a7,i4,/)

      end
