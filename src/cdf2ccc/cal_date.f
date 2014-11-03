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
      subroutine cal_date(string,dtime,iyear,imonth,iday,ihour)  

      implicit none

      include 'cdf2ccc.h'

      integer iyear,imonth,iday,ihour
      real*8  dtime
      character*(*) string
 
******
*
*AUTEUR Guy Bergeron  mai  2005
*
*     Definir la date a partir d'une valeur numerique et de l'attribut
*     "units" de la variable "time".
*
*NOTA: Ce sous-programme fonctionne uniquement pour 
*
*      time:units = "months since aaa-mm-dd hh:mm:s.s"
*      time:units = "months since aaa-mm-dd"
*      time:units = "months since aaa-mm"
*     (ou la meme chose en remplacant "months" par "years") 
*
*      et
*
*      time:units = "Season DJF MAM JJA SON"
*
*REVISIONS
*
*     Bernard Dugas mai 2008 :
*     - On suppose que "string" est toujours en minuscule
*     - Ajouter du code pour traiter l'option "years since"
*
******

      character*128 ref,date(3),dummy
      character*1 delim

      integer i,n,ni,nt,nlen,mlen
      integer aaaa,mm,dd,hh
      data dd,hh /1,0/
      data date /'01','01','01'/

*-----------------------------------------------------------------------

*     Identifier l'origine et la frequence a partir de "time units":

      n=1
      nt=len(string)
      delim=' '
      call get_string(string(n:nt),nt,delim,dummy,nlen)
      if (nlen.ne.0) then

         ref=string(1:nlen)
         ni = nt-nlen-n-1
         n=nlen+n+1
      end if

      if (ref.eq.'months' .or.
     .    ref.eq.'years') then

         call get_string(string(n:nt),ni,delim,dummy,mlen)
         if (nlen.ne.0) then

            ni = nt-mlen-n-1
            n=n+mlen+1
         endif

         call get_string(string(n:nt),ni,delim,dummy,mlen)

         delim='-'
         do i=1,3
            if (i.eq.3) delim=' '
            call get_string(string(n:nt),ni,delim,dummy,mlen)

            if (mlen.gt.0) then
               nlen=len_trim(dummy)
               if (nlen.ne.mlen)then
                  call get_string(string(n:nt),ni,':',dummy,mlen)
                  nlen=len_trim(dummy)
                  if (nlen.eq.mlen)date(i)=dummy(1:nlen)
               else
                  date(i)=dummy(1:nlen)
               endif
               n=n+mlen+1            
               ni = nt-n
            else
               goto 910
            endif
 910     enddo

*     Traduire date en integer :

         read(date(1),2000)aaaa
         read(date(2),2000)mm
         read(date(3),2000)dd

*     Evaluer la date :

         if (ref.eq.'months') then
            imonth= mm+nint(dtime)
            iyear = int((imonth-1)/12)
            imonth= imonth-(iyear*12)
            iyear = aaaa+iyear
         else
            iyear = aaaa+nint(dtime)
            imonth= mm
         endif

         iday  = dd
         ihour = hh

      else if (ref.eq.'season') then

*     Donnees de season de PRUDENCE  :

         if (dtime.gt.3.0)then          ! CNRM :dtime est nieme jour de l'annee
            iyear =1
            imonth=int(dtime/30.+.5)
         else                           ! SMHI :dtime est un indice (0,1,2,3) 
            iyear =1
            imonth=1+int(dtime*3)
         endif
         iday  = 15
         ihour = hh
      endif

*-----------------------------------------------------------------------
 2000   format(i4)
*-----------------------------------------------------------------------
      end
