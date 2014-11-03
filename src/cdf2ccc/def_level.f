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
      subroutine def_level(lablvl,operation)
      
      implicit none

      include 'cdf2ccc.h'
      include 'dimmem.h'
      include 'infomem.h'
      include 'varmem.h'
      include 'workmem.h'     

      integer lablvl(maxlev)
      character*(*) operation

******
*
* AUTEUR    Guy Bergeron   aout 2003
*
*     Definir les valeurs des niveaux "coord(zid)%vlaue" en fonction des 
*     etiquettes de niveaux "lablvl" (operation='decode') ou l'inverse 
*     (operation='encode'). L'encodage et le decodage sont fait de facon 
*     consistante avec lvcode et lvdcode (ou avec convpr pour ce qui est
*     des fichiers CMC/RPN)
*
*     En mode 'encode' la variable lablvl est une sortie du sous-programme
*     et est exactement la valeur de sortie de lvcode. En mode 'decode' la 
*     variable lablvl n'est ni une sortie ni un entre au sous-programme et
*     correspond a la sortie de lvdcode.
*
* REVISIONS
*
* B.Dugas aout '12 :
* - Ajouter le code de support du fichier valeurs_remplacement
* - Le cas level_desc='model_level_number' est traite comme
*   une instance de niveaux arbitraires.
* B.Dugas juillet '12 :
* - Ajouter le support de 'Log Pressure Hybrid Levels' (VKIND=5002)
* B.Dugas decembre '11 :
* - Ajouter 'Hybrid Height' (pour HadGEM2-ES)
* B.Dugas oct '08 :
* - Traitement des types de niveaux 'Arbitrary Levels', 
*   'Top of Atmosphere' et 'Soil Layers' (pour tous, kind=3)
* B. Dugas, ete/automne/hiver 2007 :
* - Ajouter le support des coordonnees de pression hybride et de hauteur
* - Ajouter le support des fichiers standards CMC/RPN.
*
******

      integer k,nlev

      logical, save :: decode_called = .false.
      integer, dimension(:), allocatable, save :: level_store

      integer  iun,ier
      integer, save :: nbrlev
      integer, external :: fnom,fclos
      character(128) valeurs_remplacement
      integer, allocatable,save :: replace(:)
      logical, save :: done=.false.,ex=.false.

      real    rcste,grav,rgas
      integer ip1,kind,mode 
      real*4  level(maxlev),rlev
*-----------------------------------------------------------------------
      call defcphy(grav,rgas)
      rcste=-1.*rgas*tmoyen/grav


      if (operation.eq.'encode') then           ! lablvl = fctn(coord(zid)%value)

         nlev=dim(zdid)%len

         if (decode_called .and. ccc_pktyp(1:2).eq.'SQ') then
             lablvl(1:nlev) = level_store(1:nlev)
             return
         endif

         mode = +2 ; kind = -1

         if (level_desc.eq."Surface"   .or. level_desc.eq."10 m"
     .  .or. level_desc.eq."Sea Level" .or. level_desc.eq."2 m" ) then
            
            if (ccc_pktyp(1:2).eq.'SQ') then
               if (level_desc.eq."Surface"
     .        .or. level_desc.eq."Sea Level") then
                  rval(1)=1.0 ; kind=1
               else
                  kind=0
                  if (level_desc.eq."10 m") rval(1)=10.0
                  if (level_desc.eq. "2 m") rval(1)= 2.0
               endif
               call convpr( lablvl(1),rval(1),kind,mode )
            else
               lablvl(1)=1
            endif

            return

         else if (level_desc.eq.'Gal-Chen Levels') then

            if(tmoyen.le.0.0) then
               write(6,*) " TMOYEN n'est pas definie "
               call                                  xit('def_level',-1)
            endif

            kind=21

            do k=1,nlev
               rval(k) = exp(dcoordonne(nlev-k+1,zid)/rcste)
            enddo

         else

            if (ccc_pktyp(1:2).eq.'SQ') then
               if      (coord(zid)%nattr == -1           .or.
     .                 level_desc == 'Arbitrary Levels'  .or.
     .                 level_desc == 'model_level_number'.or.
     .                 level_desc == 'Top of Atmosphere' .or.
     .                 level_desc == 'Soil Layers'     ) then
                  kind=3
               elseif (level_desc.eq.'Height') then
                  kind=0
               elseif (level_desc == 'Sigma Levels') then
                  kind=1
               elseif (level_desc == 'Pressure Levels') then
                  kind=2
               elseif (level_desc == 'Hybrid Height') then
                  kind=4
               elseif (level_desc == 'Hybrid Levels'
     .            .or. level_desc == 'Log Pressure Hybrid Levels') then
                  kind=5
               elseif (level_desc == ' ') then
                  kind=1
               endif
            endif

            do k=1,nlev
               rval(k) = dcoordonne(nlev-k+1,zid)
            enddo

         endif
         
         if (ccc_pktyp(1:2) == 'SQ' .and. kind >= 0) then

            do k=1,nlev
               rlev = rval(k)*coord(zid)%mult+coord(zid)%add
               call convpr( lablvl(k),rlev,kind,mode )
            enddo

         else

            do k=1,dim(zdid)%len               
               level(k) = (rval(k)*coord(zid)%mult+coord(zid)%add)/1000.0
            enddo
         
            call lvcode(lablvl,level,nlev)    

            if (chklvl .or.     ! assurer la consistance
     .          dim(coord(zid)%dimid(1))%len.eq.1) then
               do k=1,dim(zdid)%len 
                  if(lablvl(k).lt.0) lablvl(k)=int(1000.*level(k)+0.5)
               enddo
            endif

         endif

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Le code qui suit permets de courcircuiter tous les calculs
         ! precedents de lablvl. Ce vecteur contient les valeurs codees
         ! de la coordonnee verticale qui seront ecrites dans le fichier
         ! de sortie.
         !
         ! Des nouvelles valeurs sont lues dans un fichier texte, une
         ! valeur par ligne, les dix premiers caracteres seulement etant
         ! consideres, les espaces etant ignores ('null', suivant le
         ! format BN). Le fichier doit se trouver dans le repertoire
         ! de travail et porter un nom construit avec nom de la
         ! coordonnee verticale + '_remplacement.txt'. Il doit aussi
         ! y avoir au moins autant de lignes que de niveaux verticaux.
         !
         ! Si une de ces conditions n'est pas remplies, la substitution
         ! ne s'effectura pas. Aucune verification n'est effectuee 
         ! quand a la pertinence des valeurs lues. (BD, aout 2012)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (.not.done) then
            ! Une seule tentative de lecture est efectuee
            done=.true. ; nbrlev = dim(zdid)%len ; iun = 0
            valeurs_remplacement=trim(dim(zdid)%name)//
     .                           '_remplacement.txt'
            inquire( file=valeurs_remplacement,exist=ex ) 
            if (ex) then
               if (
     .           fnom( iun,valeurs_remplacement,'FTN+FMT',0 ) == 0) then
                 allocate( replace(nbrlev) )
                 read(iun,'(BN,I10)',end=999,err=999) replace(1:nbrlev)
                 write(6,'(/A/)')
     .                  ' !!! Lu les niveaux de remplacement sur '//
     .                           trim( valeurs_remplacement )//' !!!'
                 lablvl(1:nbrlev) = replace(1:nbrlev)
 999             ier = fclos( iun )
               endif
            endif
         else if (ex) then
            lablvl(1:nbrlev) = replace(1:nbrlev)
         endif

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      else if (operation.eq.'decode') then     ! coord(zid)%value = fctn(lablvl)

         nlev=dim(zdid)%len

         if (ccc_pktyp(1:2).eq.'SQ') then

            if (decode_called) deallocate( level_store )
                                 allocate( level_store(nlev) )

            level_store(1:nlev) = lablvl(1:nlev)

            mode=-2

            do k=1,nlev
               call convpr( lablvl(k),rlev,kind,mode )
               rval(nlev-k+1) = (rlev-coord(zid)%add)/coord(zid)%mult
            enddo
            
         else

            if (lablvl(1).gt.0) chklvl=.true.

            call lvdcode(level,lablvl,nlev)      

            do k=1,nlev
               rval(nlev-k+1)=(level(k)-coord(zid)%add)/coord(zid)%mult 
            enddo

         endif

         if (level_desc.eq.'Gal-Chen Levels') then 

            do k=1,dim(zdid)%len
               dcoordonne(k,zid)=rcste*alog(rval(k))
            enddo

         else if (level_desc.eq.'2 m') then 

            dcoordonne(1,zid)=2.0
         else if (level_desc.eq.'10 m') then 

            dcoordonne(1,zid)=10.0
         else if (level_desc.ne."Surface") then 

            do k=1,dim(zdid)%len
               dcoordonne(k,zid)=rval(k)
            enddo
         endif

         decode_called = .true.

      endif
*-----------------------------------------------------------------------
      end
