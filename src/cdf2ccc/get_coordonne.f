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
      subroutine get_coordonne (NCID)

      implicit none
      
      include 'cdf2ccc.h'      
      include 'dimmem.h'
      include 'infomem.h'
      include 'varmem.h'
      include 'workmem.h'

      integer ncid

******
*
*AUTEUR Guy Bergeron         juillet  2003
*
*     Extraire les variables coordonnees du fichier netCDF
*
*REVISIONS
*
*     Bernard Dugas, mai 2013 :
*     - Verifier la monotonicite de la coordonnee 'J' si elle est definie
*       ici et aussi, verifier que la variable INVJ est bien definie
*     - Aussi verifier les noms de variables 'longitude' et 'latitude'
*     Bernard Dugas, janv 2012 :
*     - Definir ydid, xdid lorsqu'on trouve les coordonnees correspondantes
*     Bernard Dugas, sept 2008 :
*     - definir les valeurs par defaut des coordoonnees
*       temporelle et verticale
*     - meme si coord(xid) et coord(yid) ne sont pas dans
*       la liste de variables, on tente de les lire
*       indirectement
*
******

      integer err
      integer ii,i,id,status,ngatts
      integer ndim,nvar,ngatt,unlimd
      logical linvj,ok

*-----------------------------------------------------------------------

!     valeur par defaut pour coordoonnee temporelles

      do i=1,dim(timedid)%len 
c         dcoordonne(i,tid) = i
      enddo

!     valeur par defaut pour coordoonnee verticales

      do i=1,dim(zdid)%len
         dcoordonne(i,zid) = dim(zdid)%len-i+1
      enddo

!     lire ce qui est present

      do i=1,ncoord
      do id=1,nlist
         if (coord(i)%name.eq.list(id)%name) then
          call get_vard(ncid,id,list(id)%type,dim(list(id)%dimid(1))%len
     .                                                 ,dcoordonne(1,i))
        endif
      enddo
      enddo

!     meme si coord(xid) et coord(yid) ne sont pas dans la liste
!     de variables, on peut tenter de les lire indirectement

      if (coord(xid)%nattr == -1) then

         do id=1,nlist
            if (list(id)%ndim ==  2          .and.
     .         (list(id)%name == 'longitude' .or.
     .          list(id)%name == 'lon'     ) ) then
               if (dim(list(id)%dimid(1))%name == coord(xid)%name  .and.
     .             dim(list(id)%dimid(2))%name == coord(yid)%name) then
                  call get_vard( ncid,id,list(id)%type,
     .                           dim(list(id)%dimid(1))%len,
     .                           dcoordonne(1,xid) )
                  coord(xid)%nattr = -2 ! signal que dcoordonne(:,xid) est defini
                  xdid = list(id)%dimid(1)
                  exit
               endif
            endif
         enddo

      endif

      if (coord(yid)%nattr == -1) then

         do id=1,nlist
            if (list(id)%ndim ==  2         .and.
     .         (list(id)%name == 'latitude' .or.
     .          list(id)%name == 'lat'    ) ) then
               if (dim(list(id)%dimid(1))%name == coord(xid)%name  .and.
     .             dim(list(id)%dimid(2))%name == coord(yid)%name) then
                  call get_vard( ncid,id,list(id)%type,
     .                           dim(list(id)%dimid(1))%len*
     .                           dim(list(id)%dimid(2))%len,
     .                           dval )
                  do i=1,dim(list(id)%dimid(2))%len
                     ii=(i-1)*dim(list(id)%dimid(1))%len+1
                     dcoordonne(i,yid) = dval(ii)
                  enddo
                  linvj = (dcoordonne(1,yid) > dcoordonne(2,yid))
                  do i=2,dim(list(id)%dimid(2))%len-1
                     ! Verifier la monotonicite de ces valeurs
                     ok = (dcoordonne(i,yid) > dcoordonne(i+1,yid))
                     if (linvj .neqv. ok) exit
                  enddo
                  if (invj .neqv. linvj) then
                     write(6,6001) char(7),linvj,trim( list(id)%name )
                     invj = linvj
                  endif
                  coord(yid)%nattr = -2 ! signal que dcoordonne(:,yid) est defini
                  ydid = list(id)%dimid(2)
                  exit
               endif
            endif
         enddo

      endif

*-----------------------------------------------------------------------
 6001 format(/A,'*** La cle "-invj" est re-definie a ',L1,
     .      ' apres lecture de la variable ',A/)

      end
