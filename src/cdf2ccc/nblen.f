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
      integer function nblen (string)
!
!     given a character string, nblen returns the length of the string
!     to the last non-blank character, presuming the string is left-
!     justified, i.e. if string = '   xs  j   ', nblen = 8.
!
!     called non-library routines: none
!     language: standard fortran 77
!
      integer ls,i
      character*(*) string, blank*1, null*1
      data blank /' '/
!
      null = char(0)
      nblen = 0
      ls = len(string)
      if (ls .eq. 0) return
      do i = ls, 1, -1
         if (string(i:i) .ne. blank .and. string(i:i) .ne. null) go to 2
      enddo
      return
    2 nblen = i
      return
!-----------------------------------------------------------------------
      end
