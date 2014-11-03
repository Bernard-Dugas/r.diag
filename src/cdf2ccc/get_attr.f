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
      subroutine get_attr(LENAME,nbr,mult,add,vname)

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'workmem.h'

      character*(*) lename,vname
      integer nbr
      real*8 mult,add

******
*
* AUTEUR Guy Bergeron     Juin 2003
*
*     Lecture des attributs associes a la variables "lename". "lename" peut 
*     etre soit le nom netCDF (chaine de caracteres) ou l'etiquette CCCma
*     (chaine de 4 caracteres en majuscule). En plus des attributs, le sous-
*     programme retourne "vname" qui est le nom netCDF. La lecture se fait 
*     dans le fichier "attributs_netcdf.dat".
*
* REVISIONS
*
* B. Dugas, jul 2012
* - Message informatif lorsque la definition d'un attribut
*   n'est pas  terminee par une caractere ';'
* - Les variables string et dummy passent a 600 caracteres
* - La variable cvalue passe a 512 caracteres
*
* B. Dugas, fev 2009
* - Ajouter le support des donnees de type nf_byte
*
* B. Dugas: aout 2007 :
* - Utliser numero I/O attunit plutot que "99"
* - Rembobinner attunit a l'entree et ne pas le fermer en quittant
*
******

      integer  i,n,nn,ncdf,nccc,ni,nt,nlen,mlen,vlen,idummy
      integer k,kk,nk(3),klen(3)

      logical       message
      character*1   delim
      character*4  cccname
      character*600 string,dummy
      character*512 cvalue(3)


      
      logical ok(2)
*-----------------------------------------------------------------------

      rewind attunit
      call get_name(attunit,lename,vname,cccname,mult,add,ok)

      if (.not.ok(1)) then
         call getnam( attunit,string )
         write(6,*) 'GET_ATTR : ',
     .          'Definir les attributs de "',trim(lename),
     .          '" dans "'//trim( string )//'"'
         call                                         xit('get_attr',-1)
      endif

      if(.not.ok(2)) return

 210  read(attunit,1002,end=901)string
      nt=len_trim(string)
      delim='='
      call def_name(string,nt,delim,dummy,nlen)
      if (nlen.ne.0) then

         message = .true.

         nbr=nbr+1
         attr(nbr)%name=dummy(1:nlen)          ! nom de l'attribut

         n=1
         kk=0
 301     k=index(string(n:nt),';')
         if(k.gt.0)then
            kk=kk+1
            n=n+k+1          
            nk(kk)=n
            goto 301
         else if (kk == 0 .and. message) then
            message = .false.
            cvalue(1) = 'Error reading attribut <'//attr(nbr)%name
            cvalue(1) = trim( cvalue(1) )//'>, no ending ;'
            klen(1) = len_trim( cvalue(1) )
            write(6,'(/" attribut ",A," non delimite par un ;"/)')
     .            trim( attr(nbr)%name )
         endif

         mlen=index(string,'=')+1
         delim=';'

         do i=1,kk
            ni=nt-mlen
            call get_string(string(mlen:nt),ni,delim,dummy,klen(i))
            cvalue(i)=dummy(1:klen(i))
            mlen=nk(i)
         enddo

         if (kk.gt.1) then

            read(cvalue(kk),1004) attr(nbr)%type
            attr(nbr)%len=kk-1

            do i=1,kk-1 

               if(attr(nbr)%type .eq. nf_byte) then
                  read(cvalue(i),1004) attr(nbr)%i1value(i)

               else if(attr(nbr)%type .eq. nf_short) then
                  read(cvalue(i),1004) attr(nbr)%i2value(i)

               else if(attr(nbr)%type .eq. nf_float) then
                  read(cvalue(i)(1:klen(i)),1005) attr(nbr)%rvalue(i)

               else if(attr(nbr)%type .eq. nf_double) then
                  read(cvalue(i)(1:klen(i)),1005) attr(nbr)%dvalue(i)

               else if(attr(nbr)%type .eq. nf_char .and. kk.lt.3) then
                  attr(nbr)%len=klen(i) 
                  attr(nbr)%cvalue=cvalue(i)(1:klen(i))

               else
                  call                                xit('def_cdf',-1)
               endif

            enddo
         else

            attr(nbr)%type=nf_char
            attr(nbr)%len=klen(1) 
            attr(nbr)%cvalue=cvalue(1)(1:klen(1))

         endif

         if (index(string,']').eq.0) goto 210

      end if

      if(.false.) then                                                   !debug
                                                                         !debug
         if (attr(nbr)%type .eq. nf_byte) write(6,9002)                  !debug
     .   attr(nbr)%name(1:nlen),(attr(nbr)%i1value(i),i=1,attr(nbr)%len) !debug
                                                                         !debug
         if (attr(nbr)%type .eq. nf_short) write(6,9002)                 !debug
     .   attr(nbr)%name(1:nlen),(attr(nbr)%i2value(i),i=1,attr(nbr)%len) !debug
                                                                         !debug
         if (attr(nbr)%type .eq. nf_float) write(6,9003)                 !debug
     .   attr(nbr)%name(1:nlen),(attr(nbr)%rvalue(i),i=1,attr(nbr)%len)  !debug
                                                                         !debug
         if (attr(nbr)%type .eq. nf_char) write(6,9004)                  !debug
     .   attr(nbr)%name(1:nlen),attr(nbr)%cvalue(1:attr(nbr)%len)        !debug
                                                                         !debug
 9002 format('get_attr :',A,' = ',2i8)                                   !debug
 9003 format('get_attr :',A,' = ',2f12.5)                                !debug
 9004 format('get_attr :',A,' = ',A)                                     !debug
                                                                         !debug
      endif                                                              !debug


 901  return

*-----------------------------------------------------------------------
 1002 format(a)
 1004 format(i10)
 1005 format(e12.5)
*-----------------------------------------------------------------------
      end
      
