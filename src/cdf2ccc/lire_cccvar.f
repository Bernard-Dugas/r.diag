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
      subroutine lire_cccvar3(ncid,funit,id,cccname,ilev,itime,ibuf,mss)

      implicit none

      include 'cdf2ccc.h'
      include 'dimmem.h'
      include 'infomem.h'
      include 'varmem.h'   
      include 'workmem.h'  

      logical  mss
      integer  ncid,funit,id,itime,ibuf(8),ilev
      character*4 cccname

******
*
*AUTEUR Guy Bergeron         juillet 2003
*
*     Lire dans le fichier CCCma la variable ID associee au temps ITIME.
*     Realise la conversion des unites CCCma->netCDF 
*     
*        value=(dval - add)/mult
*
*REVISIONS
*
*     Bernard Dugas juin 2013 :
*     - Tenir compte de fill_ccc_def dans le calul de variable
*     - Ajouter la variable de sortie mss qui indique que des
*       valeurs manquantes sont presentes dans les donnees
*     - Renommer a lire_cccvar3
*     Bernard Dugas juin 2012 :
*     - Renommer a lire_cccvar2 (ajout de l'argument ilev)
*     - Seulement faire les appels a def_level
*       et get_ip1_string lorsque ilev > 0
*     Bernard Dugas novembre 2009 :
*     - Ne plus faire appel a decodate2
*     - Si tlbl est vrai ibuf(2) contient maintenant un DateTimeStamp
*     Bernard Dugas fevrier 2009 :
*     - Faire appel a decodate2
*     - Tenir compte des minutes et secondes dans la
*       recherche des donnees avec les fichiers CMC/RPN
*     Bernard Dugas juillet 2008 :
*     - Tenir compte de delta_type dans le calcul de nstep
*     Bernard Dugas janvier 2008 :
*     - La variable lablvl est allouee comme champs automatique
*     - Tenir compte de la nouvelle definition des IP1 (CMC/RPN)
*     Bernard Dugas 3 aout 2007 : Tenir compte de la cle invj
*     Guy Bergeron avril 2004 : Declaration de coordonne en REAL*8
*
******

      character(len=12) ip1out
      character(len=18) dateout
      integer i,j,k,ij,ijk,ni,nj,nstep,nlev 
      real    rcste,grav,rgas,value,xtime
      integer stamp,part1,part2, err
      integer dateo,npas,deet
      integer*8, parameter :: million = 1 00 00 00
      integer*8 nstep8
      real*8  rip2,toler
      logical ok

      LOGICAL           rpn_info,rpn_debug
      COMMON  /ZZVERBO/ rpn_info
      COMMON  /ZZDEBUG/          rpn_debug

      integer lablvl(maxlev)

      integer  newdate,gethigh
      external newdate,gethigh
 
*-----------------------------------------------------------------------
      if(tlbl)then         

         call  decodate(ncid,dcoordonne(itime,tid),nstep)

      else                 

         ! Si ibuf(2)= nbre de dt

         if      (delta%type == 's') then
            nstep=nint(dcoordonne(itime,tid)        /dt)
         else if (delta%type == 'm') then
            nstep=nint(dcoordonne(itime,tid)*   60.0/dt)
         else if (delta%type == 'h') then
            nstep=nint(dcoordonne(itime,tid)* 3600.0/dt)
         else if (delta%type == 'd') then
            nstep=nint(dcoordonne(itime,tid)*86400.0/dt)
         endif

      endif

      if (level_desc.eq."Sea Level" .or. level_desc.eq."Surface") then 
         nlev=1
         lablvl(1)=1
      else
         nlev=ilev
         if (ilev > 0) then
            call def_level(lablvl,'encode')
         else
            nlev = 1
            lablvl(1)=-1
         endif
      endif

      ip1out = ' ' ; mss = .false.

      do k=1,nlev

 100     call recget(funit,' ',-1,cccname,lablvl(k),ibuf,maxpk,ok)

         if (.not.ok) then

            if (tlbl .and. rpn_info) then
               call PDATE( dateout,nstep )
               if (ilev > 0) call get_ip1_string( lablvl(k),ip1out )
               write(6,6000) trim(dateout),trim(cccname//' '//ip1out)
            endif
            call                         xit('lire_cccvar',-2)

         else if (ibuf(2) /= nstep) then

            if (ccc_pktyp(1:2) == 'PK' .or.
     .         (ccc_pktyp(1:2) == 'SQ' .and. delta%type /= 's'
     .                                 .and. delta%type /= 'm')) then

               goto 100

            else if (ccc_pktyp(1:2) == 'SQ') then

               dateo  = gethigh('DATEO',ibuf )
               npas   = gethigh('NPAS', ibuf )
               deet   = gethigh('DEET', ibuf )
               rip2   = ( npas*( DBLE( deet )/ 60) )/60

               call incdatr( stamp,dateo,rip2 )

               if (stamp /= ibuf(2)) goto 100

            endif

         endif
         
         call recup2(dval,ibuf) 

         ni=dim(xdid)%len
         nj=dim(ydid)%len

         do j=1,nj
            if (invj) then
               do i=1,ni
                  ij = i + (nj-j)*(ni+dim(xdid)%duplic) !inversion des latitudes
                  ijk= i + (j-1)*ni + (nlev-k)*ni*nj 
                  if (abs(dval(ij)-fill_ccc) >= fill_toler) then
                     variable(ijk,id)=(dval(ij)- var(id)%add)
     .                               /    var(id)%mult
                  else
                     variable(ijk,id)=fill_ccc ; mss = .true.
                  endif
               enddo
            else
               do i=1,ni
                  ij = i + (j-1)*(ni+dim(xdid)%duplic) !pas d'inversion des latitudes
                  ijk= i + (j-1)*ni + (nlev-k)*ni*nj      
                  if (abs(dval(ij)-fill_ccc) >= fill_toler) then
                     variable(ijk,id)=(dval(ij)- var(id)%add)
     .                               /    var(id)%mult
                  else
                     variable(ijk,id)=fill_ccc ; mss = .true.
                  endif
               enddo
            endif
         enddo
      enddo
*-----------------------------------------------------------------------
 6000 FORMAT(" Lire_cccvar n'a pas trouve... ",A,1x,A)

      end
