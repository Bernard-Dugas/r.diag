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
      subroutine def_date2( date,aaaa,mm,jj,hh,mn,ss,operation )

      implicit none

      integer date,aaaa,mm,jj,hh,mn,ss
      character*(*)   operation
      character*(128) string
      
******
*
*AUTEUR Guy Bergeron   avril 2004
*
*     Permet de definir la "date" dans le format datetimestamp a partir
*     des valeur aaaa, mm, jj, hh, mn, ss et ff si operation='encode'
*     ou l'inverse si operation='decode'.
*
*REVISIONS
*
*     B.Dugas mai 2008 :
*     - Ajouter les arguments "mn" et "ss"
*     - Changer de nom a def_date2
*     - coder/decoder avec newdate
*     B.Dugas avril 2008 : Utliser "call xit" en cas d'erreur     
*
******

      logical  ok
      integer  part1,part2

      integer  newdate,datchek
      external newdate

*-----------------------------------------------------------------------
*     notez que part1=aaaammjj, part2=hhmnsscc

      if (operation.eq.'decode')then

         datchek = newdate( date, part1,part2, -3 )

         if (datchek.ne.0) then
            datchek = newdate( date, part1,part2, -5 )
            if (datchek.ne.0) call                 xit('def_date2', -1 )
         endif

         write(string,2222) part1
         read (string,2223) aaaa,mm,jj

         write(string,2222) part2
         read (string,2224) hh,mn,ss

      else if (operation.eq.'encode')then

         write(string,2225) hh,mn,ss
         read (string,2222) part2

         write(string,2223) aaaa,mm,jj
         read (string,2222) part1

         datchek = newdate( date, part1,part2, +3 )

         if (datchek.ne.0) then
            datchek = newdate( date, part1,part2, +5 )
            if (datchek.ne.0) call                 xit('def_date2', -2 )
         endif

      else 

         write(6,6000)
         call xit('def_date2',-3)

      endif
*-----------------------------------------------------------------------
 2222 format(i8.8)
 2223 format(i4.4,2i2.2)
 2224 format(3i2.2)
 2225 format(3i2.2,'00')
 6000 format(/,"DEF_DATE2 : 'operation' mal definie ",/)
*-----------------------------------------------------------------------
      end
