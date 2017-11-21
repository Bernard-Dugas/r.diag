#     if !defined (taille_entete)
#         define   taille_entete 32
#     endif
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
      subroutine inq_file2( FUNIT,PHIS_UNIT, nvar,ndim,maxlens, &
                                             max1ds,nlev,ntime )
      
      use diag_toc, only: GetTOC, LirTOC
      implicit none
      
      include 'cdf2ccc.h'

      integer funit,phis_unit,ndim,nvar,maxlens,max1ds,nlev,ntime

!*****
!
!AUTEUR Guy Bergeron    aout 2003
!      
!     Defini le nombre de variable(navr),le nombre de dimension (ndim)
!     maxlen,max1d et nlev
!
! NOTA :
!     Le type derive "infvar" est utilise pour faire la liste des variables du
!     fichier d'entre (CCCma). Il n'y a que name et var_ok qui sont affectees
!     pour le moment. Eventuellement, il pourrait y avoir aussi ndim, dimid
!     
!         
!REVISIONS
!
!  B. Dugas mai 2017 : 
!  - Correction de l'initialisation
!    de la variable locale unique_L
!  B. Dugas mai 2017 : 
!  - Convertir en fichier .F90 pour traiter
!    le macro taille_entete avec s.f90
!  B. Dugas avril 2017 :
!  - Ajouter support explicite de VKIND=5003,5004,5005
!  - Ajouter une sortie d'erreur -13 si une grille GEM4
!    de type 'U' (Ying/Yang) est detectee
!  B. Dugas janvier 2017 :
!  - Utiliser GETSAMPLZ pour definir IP3, i.e. le nombre 
!    d'echantillons correspondant a une moyenne temporelles
!    a partir des informations dans la section haute de IBUF
!  B. Dugas mars 2015 :
!  - Correctif lie ax appels GET_TOC/VGRID.
!  B. Dugas juillet 2013 :
!  - Verifier la presence de valeurs manquantes dans
!    les fichiers CMC/RPN via un appel a MISPAR
!  B. Dugas mai 2013 :
!  - Correction du declencheur de l'erreur -8 (inq_file2)
!  B. Dugas octobre 2012 :
!  - Ajouter la possibilite d'utiliser les parametres de la
!    coordonnee hybride lus en arguments lorsque HY est manquant
!  B. Dugas aout 2012 :
!  - BugFix associe au traitement de level_desc='Hybrid Levels'
!  - Definir infvar(:)%unique_L=T lorsque nlev=1
!    et que les label(1,:) sont identiques
!  - Pour VKIND=5002, verifier le deuxieme jeux de niveaux si on les
!    trouve pas tous dans le premier avant de faire un appel a XIT(-9)
!  B. Dugas juillet 2012 :
!  - Effectuer plus de tests de coherence
!  - Ajouter le support de 'Log Pressure Hybrid Levels' (VKIND=5002)
!    de meme que tenir compte des cas 'gem2' et 'gem3' (voir r.diag)
!  - Sauver tous les IP3 et DATEO des variables dans infvar
!  B. Dugas juin 2012 :
!  - Plus grand usage de INFVAR (plus de contenus)
!  - Meilleure identifications des niveaux et des
!    pas-de-temps associes a chaque variable
!  - Deplacer le traitement des coordonnees verticales
!    de RPN_PARAMS a INQ_FILE pour les fichiers CMC/RPN.
!  - Ajouter le support des TocToc (!!) via des appels DIAG_TOC
!  - Ajouter les calculs des vecteurs AP et B ('Hybrid Levels')
!  - Tenir compte de variables logiques GEM2 et GEM3
!  B. Dugas mai 2012 : Tenir compte de time_bnds dans RPN_PARAMS
!    lors de la definition de la variable ladate.
!  B. Dugas decembre 2011 :  Ajouter kind=4 => 'Hybrid Height'.
!  Bernard Dugas juillet 2011 :  Definir DATE_CONVERSION_FACTOR
!    pour les fichier de type CCC, selon la valeur de la cle tlbl.
!  Bernard Dugas juin 2008 :  On construit une "date" tenant compte
!    des minutes et secondes pour distinguer les differents pas-de-temps
!    avec la fonction grosse_date. Mais IBUF(2) est utilise pour CCCma.
!  Bernard Dugas avril 2008 :
!  - Enlever l'include de 'machtyp.h'.
!  - Conserver le facteur de compaction dans infvar(:)%npack
!  - Ajouter l'argument PHIS_UNIT pour le fichier de champs de montagne
!    et changer le nom de la routine a inq_file2
!  Bernard Dugas 3 mai 2007 : Ajouter le support des fichiers RPN/CMC
!  Anne Frigon 10 Aout 2004 : Correction de SPEC au lieu de spec pour IF IBUF(1)
!  Guy Bergeron juillet 2004 : Hybrid height coordiante
!  Guy Bergeron avril   2004 : introduction de machtyp.h 
!
!*****

      integer, parameter :: head = taille_entete

      integer(8) :: nstep, DATE_CONVERSION_FACTOR
      integer(8), allocatable :: cstep(:,:)
      integer,    allocatable :: label(:,:),ip3(:,:),dateo(:,:)

      integer ibuf(head)
      integer ni,nj,nk,nt,dummy,jm,jt,kind,ig1,itime,FirstVkind
      integer i,j,count,cvar, nlevs,ntim,ig2
      integer :: firstV=-1,firstT=-1
      real(8) :: RVALUE,REPSIL

      integer, parameter :: max_steps = 10000 ! nombre maximal de pas-de-temps

      ! Parametres definissant les differentes saveurs de la
      ! coordonnee hybride utilisee dans les fichiers CMC/RPN 
      real     :: hyb_pt0,hyb_pref0,hyb_r0
      real(8)  :: ptop8,pref8,r1,r2 
      real     :: cofa(max_levs)
      real     :: cofb(max_levs)
      real(8), pointer :: ca_m(:)
      real(8), pointer :: cb_m(:)
      real(8), pointer :: ca_t(:)
      real(8), pointer :: cb_t(:)
      integer, pointer :: vipm(:)
      integer, pointer :: vipt(:)
      real     :: eta (max_vars),valeur
      integer  :: TOC_NBR
      integer  :: name(max_vars),len(max_dims)

      character(len=4) stat_phis_unit
      character(len=4) type,grtyp
      logical, save :: DoIt=.true.
      logical ok,ok2,error,unique_L

      integer maxwrd
      common /maxwrd/ maxwrd

      real    hi,lo
      integer length,nwds,nblen,opack, rk,ipm

      character(len=4), external :: gethic
      integer,          external :: gethigh,getsamplz
      integer(8),       external :: grosse_date

!-----------------------------------------------------------------------
      nullify( ca_m,cb_m,ca_t,cb_t,vipm,vipt )

      allocate( label(max_levs,max_vars),cstep(0:max_steps,max_vars) )
      if (ccc_pktyp(1:2) == 'SQ') &
        allocate( ip3(max_levs,max_vars),dateo(max_steps,max_vars) )

      cstep = -1
      opack = -64

!NOTA: Les variables ni,nj,nk sont ici des indices et non des dimensions

      do i=1,max_dims
         len(i)=1
      enddo

      ndim=0
      count=0

      if (ccc_pktyp(1:2).ne.'SQ') then
         call GET_DCF( DATE_CONVERSION_FACTOR )
         if (tlbl) then
            if (DATE_CONVERSION_FACTOR == 0_8 ) &
                DATE_CONVERSION_FACTOR  = 1000000
         else
            if (DATE_CONVERSION_FACTOR == 0_8 ) &
                DATE_CONVERSION_FACTOR  = -1
         endif
         call SET_DCF( DATE_CONVERSION_FACTOR )
      endif

      call suivant( funit,ibuf,error,ok )

      if(cdf2_mode.eq.'cdf2rpn') &
         call rpn_params( funit,phis_unit,ibuf )

      grtyp = gethic('GRTYP',ibuf )

      cvar=1
      nvar=cvar
      nlev=1
      maxpk=0

      write(type,9901) ibuf(1)

      if (type /= 'SPEC') then
         if (count == 0) then
            ni=ndim+1
            nj=ndim+2
            ndim=ndim+2
         endif
      else
         call                                      xit('inq_file2', -1 )
      endif

      ndim=ndim+1 ; nk=ndim
      ndim=ndim+1 ; nt=ndim

      label(nlev,cvar)=ibuf(4) ; name(cvar)=ibuf(3)
      cstep(1,cvar)=grosse_date( ibuf )
      if (ccc_pktyp(1:2) == 'SQ') then
         ip3  (1,cvar) = getsamplz( rk,hi,lo,ipm, ibuf )
         dateo(1,cvar) = gethigh( 'DATEO',ibuf )
         if (.not.(miss_ccc_def .or. fill_ccc_def)) then
            call MISPAR( ok2,RVALUE,REPSIL )
            if (ok2) then
               fill_ccc_def = .true.
               fill_ccc     = RVALUE
               fill_toler   = ABS( RVALUE*REPSIL )
            endif
         endif
      endif

      write(infvar(cvar)%name,9901)ibuf(3)  
      if (npack == 999) then
         infvar(cvar)%npack = ibuf(8)
      else
         infvar(cvar)%npack = npack
      endif
      infvar(cvar)%len(ni) = ibuf(5)
      infvar(cvar)%len(nj) = ibuf(6)
      infvar(cvar)%len(nk) = 1
      infvar(cvar)%len(nt) = 1
      infvar(cvar)%ndim    = ndim
      infvar(cvar)%var_ok = .true.

 100  continue

!     definissons la longueur (maxpk) du buffer d'entre (ibuf)

      call lblchk(length,nwds,opack,ibuf)    
      if(maxpk.lt.length-head) maxpk=length-head

      if(error.or..not.ok)then
         if(count.eq.0)then
            call                                   xit('inq_file2', -2 )
         else
            call precede(funit,-1)
            goto 500
         endif
      endif
!
      len(ni)=max( ibuf(5),len(ni) )
      len(nj)=max( ibuf(6),len(nj) )
!
      do cvar=1,nvar
         if (name(cvar) == ibuf(3)) exit
      enddo
      ! cvar contient maintenant l'indice d'une variable
      ! connue ou bien nvar+1 si celle-ci est nouvelle
      if(cvar > nvar) then
         nvar=cvar
         name(nvar)=ibuf(3)
         write(infvar(nvar)%name,9901)ibuf(3)      
         infvar(nvar)%len(ni) = ibuf(5)
         infvar(nvar)%len(nj) = ibuf(6)
         infvar(nvar)%npack   = ibuf(8)
         infvar(nvar)%var_ok  = .true.
      endif

      if(cdf2_mode                        .eq.'cdf2cc' .and. &
         label(infvar(cvar)%len(nk),cvar) .eq. 1       .and. &
         ibuf(4)                          .ne. 1     )       &
         label(infvar(cvar)%len(nk),cvar)  =   ibuf(4)         

      ! Nouveau niveau ?
      ok=.false.
      do i=1,infvar(cvar)%len(nk)
         if (label(i,cvar) == ibuf(4)) ok=.true.
      enddo
      if(.not.ok) then
         infvar(cvar)%len(nk)=infvar(cvar)%len(nk)+1
         label(infvar(cvar)%len(nk),cvar)=ibuf(4)
      endif
      len(nk)=max( infvar(cvar)%len(nk),len(nk) )

      ! Nouveau pas-de-temps ?
      ok=.false.
      nstep = grosse_date( ibuf )
      do i=1,infvar(cvar)%len(nt)
         if (cstep(i,cvar) == nstep) ok=.true.
      enddo
      if(.not.ok) then
         itime=infvar(cvar)%len(nt)+1
         infvar(cvar)%len(nt)=itime
         cstep(itime,cvar) = nstep
         if (ccc_pktyp(1:2) == 'SQ') then
            ip3  (itime,cvar) = getsamplz( rk,hi,lo,ipm, ibuf )
            dateo(itime,cvar) = gethigh( 'DATEO',ibuf )
            if (.not.(miss_ccc_def .or. fill_ccc_def)) then
               call MISPAR( ok2,RVALUE,REPSIL )
               if (ok2) then
                  fill_ccc_def = .true.
                  fill_ccc     = RVALUE
                  fill_toler   = ABS( RVALUE*REPSIL )
               endif
            endif
         endif
      endif
      len(nt)=max( infvar(cvar)%len(nt),len(nt) )

      count=count+1
      call suivant( funit,ibuf,error,ok )

      goto 100

 500  continue

      max1d=1
      maxlens=len(ni)*len(nj)*len(nk)
      maxwrd=maxlens

      do i=1,ndim
         if(max1d.lt.len(i))max1d=len(i)
      enddo

      nlev=len(nk)
      ntime=len(nt)

      if(nlev == 1) then
         ! Un seul niveau defini pour toutes les variables ?
         unique_L=.true.
         do cvar=2,nvar
            if(label(1,cvar) /= label(1,1)) then
               unique_L=.false.
               exit
            endif
         enddo
      else
         unique_L=.false.
      endif

      if(.false.) write(6,6600)nvar,ndim,nlev,ntime ! Debug

      do cvar=1,nvar ! Sauver les niveaux et les pas-de-temps

         ntim = infvar(cvar)%len(nt)
         if(ntim /= ntime .and. ntim /= 1) then
            write(6,6003) ntime,ntim,infvar(cvar)%name
            call                                   xit('inq_file2', -3 )
         endif

         nlevs = infvar(cvar)%len(nk)
         if(nlevs /= nlev .and. nlevs /= 1) then
            write(6,6004) nlev,nlevs,infvar(cvar)%name
            call                                   xit('inq_file2', -4 )
         endif

         allocate( infvar(cvar)%pas(ntim), infvar(cvar)%niv(nlevs) )
         infvar(cvar)%pas  (1:ntim )=cstep(1:ntim ,cvar)
         infvar(cvar)%niv  (1:nlevs)=label(1:nlevs,cvar)
         infvar(cvar)%tid=nt ; infvar(cvar)%zid=nk
         infvar(cvar)%unique_L=unique_L

         if (ccc_pktyp(1:2) == 'SQ') then
            allocate( infvar(cvar)%ip3(ntim),infvar(cvar)%dateo(ntim) )
            infvar(cvar)%dateo(1:ntim )=dateo(1:ntim ,cvar)
            infvar(cvar)%ip3  (1:ntim )=  ip3(1:ntim ,cvar)
         endif

         if(cvar > 1) then
            if(infvar(1)%len(ni) /= infvar(cvar)%len(ni)  .or. &
               infvar(1)%len(nj) /= infvar(cvar)%len(nj)) then
               write(6,'(/A/)') 'Tailles horizontales multiples'
               call                                xit('inq_file2', -5 )
            endif
         endif

         if (ntim > 1) then ! Verifier l'unicite de la coordonnee temporelle
            if(firstT == -1) then
               firstT = cvar
            else
               do i=1,ntim
                  if(cstep(i,cvar) /= cstep(i,firstT)) then
                     write(6,'(/A/)') 'coordonnee temporelle '// &
                                      'non unique...'
                     call                          xit('inq_file2', -6 )
                  endif
               enddo
               if (ccc_pktyp(1:2) == 'SQ') then
                  do i=1,ntim
                     if(dateo(i,cvar) /= dateo(i,firstT) &
                    .or.  ip3(i,cvar) /=   ip3(i,firstT)) then
                        ! Desallouer les time_bnds automatiques
                        write(6,'(/A/)') 'dateo,ip3 non uniques...'
                        dateo = -1 ; ip3 = -1
                        infvar(cvar)%dateo = -1
                        infvar(cvar)%ip3   = -1
                        exit
                     endif
                  enddo
               endif
            endif
         endif

         if(nlevs > 1 .or. infvar(cvar)%unique_L) then
            ! Verifier l'unicite de la coordonnee verticale
            if(firstV == -1) then
               firstV = cvar
               call convpr(label(1,cvar), valeur, FirstVkind, -1 )
               do i=2,nlev
                  call convpr( label(i,cvar), valeur, vkind, -1 )
                  if (vkind /= FirstVkind) then
                     write(6,'(/A/)') 'Type de coordonnee verticale '// &
                                      'non unique...'
                     call                          xit('inq_file2', -7 )
                  endif
               enddo
            else
               do i=1,nlev
                  if(label(i,cvar) /= label(i,firstV)) then
                     write(6,'(/A/)') 'Coordonnee verticale '// &
                                      'non unique...'
                     call                          xit('inq_file2', -8 )
                  endif
               enddo
            endif
         endif

         ! 1) Definir level_desc si ce n'est deja fait et
         ! 2) Possiblement identifier les A et B de la coordonnee hybride

         if(DoIt .and. ccc_pktyp(1:2) == 'SQ' .and. &
           (nlevs > 1 .or. infvar(cvar)%unique_L)) then

            DoIt = .false. ! Passer ici une seule fois

  600       if (level_desc.eq.' ') then
               if (FirstVkind ==  0) level_desc = 'Height'
               if (FirstVkind ==  1) level_desc = 'Sigma Levels'
               if (FirstVkind ==  2) level_desc = 'Pressure Levels'
               if (FirstVkind ==  3) level_desc = 'Arbitrary Levels'
               if (FirstVkind ==  4) level_desc = 'Hybrid Height'
               if (FirstVkind ==  5) level_desc = 'Hybrid Levels'
               if (FirstVkind ==  6) level_desc = 'Theta Levels'
               if (FirstVkind == 21) level_desc = 'Gal-Chen Levels'
            endif

            if (level_desc == 'Height'          ) vkind =  0
            if (level_desc == 'Sigma Levels'    ) vkind =  1
            if (level_desc == 'Pressure Levels' ) vkind =  2
            if (level_desc == 'Arbitrary Levels') vkind =  3
            if (level_desc == 'Hybrid Height'   ) vkind =  4
            if (level_desc == 'Theta Levels'    ) vkind =  6
            if (level_desc == 'Gal-Chen Levels' ) vkind = 21

            if (level_desc == 'Log Pressure Hybrid Levels' &
           .or. level_desc == 'Hybrid Levels'   ) vkind =  5

            if (vkind /= FirstVkind) then
               level_desc = ' '
               goto 600
            endif

            if (gem2 .or. gem3) vkind = 1 ! Forcer la recherche de 'HY'

            if (vkind == 1) then ! Verifier si 'HY' existe ?
               call lirpt( funit )
               call getpt( funit , hyb_pt0,hyb_pref0,hyb_r0 )
               if (hyb_pt0  /= -1.0) then
                  level_desc = 'Hybrid Levels'
                  HYB_PT     = hyb_pt0
                  HYB_PREF   = hyb_pref0
                  HYB_R      = hyb_r0
               else if ((gem2 .or.                      &
                         gem3 .or.                      &
                         level_desc == 'Hybrid Levels') &
                         .and.                          &
                        (HYB_PT     == -1.0 .or.        &
                         HYB_PREF   == -1.0 .or.        &
                         HYB_R      == -1.0 )) then
                  write(6,'(/A/)') 'coordonnee hybride GEM, aucun HY'
                  call                             xit('inq_file2', -8 )
               endif
            endif
            if (vkind == 5) then
               level_desc = 'Hybrid Levels'
               call lirpt( funit )
               call getpt( funit , hyb_pt0,hyb_pref0,hyb_r0 )
               if (hyb_pt0 == -1.0) then ! LITPT/GETPT did not WORK.
                  call LirToc( funit, TOC_NBR )
                  if (toc_nbr > 0) then
                     if (.NOT.(grtyp == 'Z' .OR. grtyp == 'U')) then
                        ig1 = -1 ; ig2 = -1
                     else if (grtyp == 'U') then
                      write(6,'(/A/)') 'grille Ying/Yang non supportee.'
                      call                         xit('inq_file2', -13)
                     else
                        ig1 = gethigh('IG1',ibuf )
                        ig2 = gethigh('IG2',ibuf )
                     endif
                     call gettoc( funit,'VER',vkind  , ig1,ig2 )
                  else if                   &
                     (HYB_PT   == -1.0 .or. &
                      HYB_PREF == -1.0 .or. &
                      HYB_R    == -1.0 ) then
                     write(6,'(/A/)') 'coordonnee hybride, '// &
                          'aucun HY ou !!'
                     call                          xit('inq_file2', -8 )
                  endif
               else
                  HYB_PT   = hyb_pt0
                  HYB_PREF = hyb_pref0
                  HYB_R    = hyb_r0
               endif
               if (vkind == 5001) then

                  call gettoc( funit, 'NK'   ,nkm   , ig1,ig2 )
                  call gettoc( funit, 'PTOP' ,ptop8 , ig1,ig2 )
                  call gettoc( funit, 'PREF' ,pref8 , ig1,ig2 )

                  call gettoc( funit, 'RC_1' ,r1    , ig1,ig2 )

                  call gettoc( funit, 'VIPM' ,vipm  , ig1,ig2 )

                  call gettoc( funit, 'COFA' ,ca_m  , ig1,ig2 )
                  call gettoc( funit, 'COFB' ,cb_m  , ig1,ig2 )

                  HYB_PT = ptop8 ; HYB_PREF = pref8 ; HYB_R = r1

               endif
               if(vkind == 5002 .or. vkind == 5003 &
              .or.vkind == 5004 .or. vkind == 5005) then

                  level_desc = 'Log Pressure Hybrid Levels'

                  call gettoc( funit, 'NKM'  ,nkm   , ig1,ig2 )
                  call gettoc( funit, 'NKT'  ,nkt   , ig1,ig2 )
                  call gettoc( funit, 'PTOP' ,ptop8 , ig1,ig2 )
                  call gettoc( funit, 'PREF' ,pref8 , ig1,ig2 )

                  call gettoc( funit, 'RC_1' ,r1    , ig1,ig2 )
                  call gettoc( funit, 'RC_2' ,r2    , ig1,ig2 )

                  call gettoc( funit, 'VIPM' ,vipm  , ig1,ig2 )
                  call gettoc( funit, 'VIPT' ,vipt  , ig1,ig2 )

                  call gettoc( funit, 'CA_M' ,ca_m  , ig1,ig2 )
                  call gettoc( funit, 'CB_M' ,cb_m  , ig1,ig2 )
                  call gettoc( funit, 'CA_T' ,ca_t  , ig1,ig2 )
                  call gettoc( funit, 'CB_T' ,cb_t  , ig1,ig2 )

                  HYB_PT = ptop8 ; HYB_PREF = pref8
                  HYB_R  = r1    ; HYB_R2   = r2

               endif
            endif

            if(vkind == 5 .or. &
              (vkind == 1 .and. level_desc == 'Hybrid Levels')) then

               do i=1,nlev
                  call convpr( label(i,cvar), eta(i), kind, -1 )
               enddo

               ! La corrdonnee GEM3 doit etre de-normalisee
               if (gem3) eta(1:nlev) = eta(1:nlev) &
                 + ( 1.0 - eta(1:nlev) ) * hyb_pt / hyb_pref

               call genab( cofa,cofb,eta,hyb_pt,hyb_pref,hyb_r,nlev )

               cap(1:nlev) = cofa(1:nlev)
               cb (1:nlev) = cofb(1:nlev)

            endif
                  

            if (level_desc.eq.'Gal-Chen Levels') then
               call getstat( phis_unit,stat_phis_unit )
               if (stat_phis_unit.ne.'OLD') then
                  write(6,6001) '-phis "PHIS FILE NAME" ?'
                  call                             xit('lire_arg', -35)
               endif
            endif

            if (level_desc.eq.' ') then
               write(6,6002) vkind
               call                                xit('lire_arg', -8 )
            endif

            if(vkind == 5002 .or. vkind == 5003 &
           .or.vkind == 5004 .or. vkind == 5005) then

               jm = 0
               jt = 0

               do j=1,nlev ! On cherche A et B dans les niveaux dynamiques ?
                  do i=1,nkm
                     if (vipm(i) == label(j,cvar)) then
                        jm = jm+1
                        cap(jm) = exp( ca_m(i) )
                        cb (jm) =      cb_m(i)
                        cycle
                     endif
                  enddo
               enddo
               if (jm == nlev) cycle ! On a trouve tous les A et B momentums !

               do j=1,nlev ! On cherche A et B dans les niveaux thermodynamiques ?
                  do i=1,nkt
                     if (vipt(i) == label(j,cvar)) then
                        jt = jt+1
                        cap(jt) = exp( ca_t(i) )
                        cb (jt) =      cb_t(i)
                        cycle
                     endif
                  enddo
               enddo
               if (jt == nlev) cycle ! On a trouve tous les A et B thermodynamiques !

               if((jm > 0 .and. jm /= nlev)  .or. &
                  (jt > 0 .and. jt /= nlev)) then ! Il en manque !
                  write(6,'(/A/)') 'Il manque des A et/ou B hybrides'
                  call                             xit('inq_file2', -9 )
               endif

            else if(vkind == 5001) then

               jt = 0
               do j=1,nlev      ! On cherche A et B pour nos niveaux
                  do jm=1,nkm
                     if (vipm(jm) == label(j,cvar)) then
                        jt = jt+1
                        cap(jt) = ca_m(jm)
                        cb (jt) = cb_m(jm)
                        cycle
                     endif
                  enddo
               enddo

               if(jt > 0 .and. jt /= nlev) then ! Il en manque toujours ?
                  write(6,'(/A/)') 'Il manque des A et/ou B hybrides'
                  call                             xit('inq_file2', -9 )
               endif

            endif

         endif

      enddo

 1000 deallocate( label,cstep )
      if (ccc_pktyp(1:2) == 'SQ') deallocate( ip3,dateo )

      return

!-----------------------------------------------------------------------
 6001 format(/" A l'appel, definir : ",a/)
 6002 format(/' Coordonne verticale: VKIND=',I2,' non supporte.'/)
 6003 format(' INQ_file : expected ',I3,' timesteps. Only found ', &
             I3,' for ',A/)
 6004 format(' INQ_file : expected ',I3,' levels. Only found ',I3, &
             ' for ',A/)
 6600 format(' nvar= ',i5,' ndim= ',i5,' nlev= ',i5,' ntime= ',i5) !debug
 9901 format(a4)
!-----------------------------------------------------------------------

      end subroutine inq_file2
      subroutine rpn_params( funit,phis_unit,ibuf )

      use      stats_signatures

      implicit none
      
      include 'cdf2ccc.h'
      include 'ztypmem.h'

      integer, parameter :: head = taille_entete

      integer funit,phis_unit,ibuf(head)
!
!     Bernard Dugas mai 2007
!      
!     Definit certains des parametres d'entree a l'aide des
!     descripteurs de grilles contenus dans le premier
!     enregistrement d'un fichier RPN/CMC.
!
!REVISIONS
!
!  Bernard Dugas nov 2017 : Comparer TYPVAR(2:2) avec stats_signatures
!  Bernard Dugas jan 2017 : utiliser getmsamplz pour definir IP3
!  Bernard Dugas mai 2008 : ladate est maintenant en format date-time-stamp
!  Bernard Dugas nov 2008 : invj = .not.invj pour les grilles de type Z
!  Bernard Dugas mai 2009 : tenir compte de GRTYP=Z et ZTYP=N ou S (en plus de E)
!  Bernard Dugas nov 2009 : appel a SET_DCF(-1) pour cccma si tlbl est faux
!
!*****

      CHARACTER(len=4), external :: GETHIC
      INTEGER,external :: GETHIGH,NEWDATE,GETSAMPLZ

      
      real(8)      DELTAT

      character    ZTYP,GRTYP,TYPVR*2,CELLM*25

      integer      DATEO, IP1, IP3, dtpr, tmpr, &
                   datchek, i, ni, nj, rk, ipm, &
                   ZIG1, ZIG2, ZIG3, ZIG4,      &
                   IG1, IG2, IG3, IG4,          &
                   DEET, NPAS, NHOUR

      real         olat    , olon  , dlat  , dlon  , &
                   pi      , pj    , d60   , dgrw  , &
                   dlat1   , dlon1 , dlat2 , dlon2 , &
                   latproj , nhem  , dx,dy , lo,hi

!-----------------------------------------------------------------------

      if (ccc_pktyp(1:2).ne.'SQ') then

!        Fichier CCCma...
!        Tenir compte des codes d'erreurs qui ont ete
!        non verifies dans lire_arg en mode cdf2rpn

         if (ladate.lt.0) then
            write(6,6001)' -dateo "yyyymmddhh" '
            call                                   xit('lire_arg', -5 )
         endif

         if(.not.tlbl .and. dt.le.0.0) then
            write(6,6001)' -dt "pas de temps" ? '
            call                                   xit('lire_arg', -6 )
         endif
         
         if (level_desc.eq.' ') then
            write(6,6001) ' -lev "level_desc" ?'
            call                                   xit('lire_arg', -8 )
         endif

         if (project%name.eq.' '            .or. &
             project%name.eq.'rotated_pole' .or. &
             project%name.eq.'rotated_latitude_longitude') then
            write(6,6001) ' -grid "grid_desc" ?'
            call                                   xit('lire_arg',-11 )
         endif

         return

      endif

!     Fichier RPN/CMC...
!     Lire certains parametres dans la section haute de ibuf(head)

      TYPVR = GETHIC('TYPVAR', IBUF )

      DATEO = GETHIGH('DATEO', IBUF )
!CCC  IP1   = GETHIGH( 'IP1' , IBUF )
      IP1   = IBUF(4)
!CCC  IP3   = GETHIGH( 'IP2' , IBUF )
      IP3   = GETSAMPLZ( RK,HI,LO,IPM, IBUF )
      DEET  = GETHIGH('DEET' , IBUF )
      NPAS  = GETHIGH('NPAS' , IBUF )

      GRTYP = GETHIC('GRTYP' , IBUF )
      IG1   = GETHIGH( 'IG1' , IBUF )
      IG2   = GETHIGH( 'IG2' , IBUF )
      IG3   = GETHIGH( 'IG3' , IBUF )
      IG4   = GETHIGH( 'IG4' , IBUF )

!CCC  invj = .true.
      tlbl = .true.

      ni   = ibuf(5)
      nj   = ibuf(6)

!     Verifier que ladate contient un date-time-stamp valide.

      if (ladate.eq.-1) ladate = DATEO

      datchek = newdate( ladate, dtpr,tmpr, -3 )
      if (datchek.ne.0) then
         datchek = newdate( ladate, dtpr,tmpr, -5 )
         if (datchek.ne.0) call                    xit('lire_arg', -5 )
      endif

      datchek = newdate( DATEO, dtpr,tmpr, -3 )

      nhour  =  tmpr /  1000000
      DELTAT = (DBLE( DEET )/3600.)*NPAS

      ! Verifier la presence d'une signature
      ! d'operations temporelles dans TYPVAR ?

      CELLM = ' '

      if (TYPVR(2:2) == time_mean_signature) CELLM = 'time: mean'
      if (TYPVR(2:2) == variance_signature ) CELLM = 'time: variance'
      if (TYPVR(2:2) == median_signature   ) CELLM = 'time: median'
      if (TYPVR(2:2) == stdev_signature    ) CELLM = 'time: standard_deviation'
      if (TYPVR(2:2) == timmax_signature   ) CELLM = 'time: maximum'
      if (TYPVR(2:2) == timmin_signature   ) CELLM = 'time: minimum'

      if (dtsize > 0.00001 .or. &
         (IP3 > 1 .and. DELTAT > 0.01_8)) then

!        Calculer les limites temporelles associees
!        a des moyennes de plusieurs echantillons et
!        modifier ladate en consequence

         if (ladate == DATEO) then

            if (dtsize > 0.00001) then
               call incdatr( ladate, DATEO, DELTAT-dtsize )
            else if (nhour > 0) then
               if (nint( deltat/(ip3-1) ) == nhour) then
                  tmpr    = tmpr - nhour * 1000000
                  datchek = newdate( ladate, dtpr,tmpr, +3 )
               endif
            endif

            datchek = newdate( ladate, dtpr,tmpr, -3 )
            if (datchek /= 0) then
               datchek = newdate( ladate, dtpr,tmpr, -5 )
               if (datchek /= 0) call              xit('lire_arg', -6 )
            endif

         endif

         time_bnds_L = .true.

         ! Possiblement sauver ce qui se trouvait dans TYPVAR(2:2)
         if (cell_method == '?' .and. CELLM /= ' ') cell_method = CELLM

      endif

!     Toujour re-definir project%name et al.

      if (GRTYP.eq.'G') then

         project%len  = 1
         project%name = 'gaussian'

         if (IG2.eq.1) invj = .not.invj
         if (IG1.ne.0) then
            write(6,'(/A/)') ' Grilles hemispheriques G non supportees '
            call                                   xit('inq_file2',-10 )
         endif

      else &
      if (GRTYP.eq.'A' .or. GRTYP.eq.'B') then

         project%len  = 1
         project%name = 'lon/lat global'//' '//GRTYP

         if (IG2.eq.1) invj = .not.invj
         if (IG1.ne.0) then
            write(6,'(/A/)') &
                   ' Grilles hemispheriques '//GRTYP//' non supportees '
            call                                   xit('inq_file2',-11 )
         endif

      else &
      if (GRTYP.eq.'L') then

         project%len  = 5
         project%name = 'lon/lat regional'

         CALL CIGAXG( GRTYP, olat,olon,dlat,dlon, &
                             IG1, IG2, IG3, IG4 )

         project%nampar(1) = '0lon' ; project%value(1) = olon
         project%nampar(2) = '0lat' ; project%value(2) = olat
         project%nampar(3) = 'dlon' ; project%value(3) = dlon
         project%nampar(4) = 'dlat' ; project%value(4) = dlat

      else &
      if (GRTYP.EQ.'N' .or. GRTYP.eq.'S') then

         project%len  = 8
         project%name = 'polar_stereographic'

         CALL CIGAXG( GRTYP, pi , pj , d60, dgrw , &
                             IG1, IG2, IG3, IG4  )

         project%nampar(1) = 'pi'   ; project%value(1) = pi
         project%nampar(2) = 'pj'   ; project%value(2) = pj
         project%nampar(3) = 'dgrw' ; project%value(3) = dgrw
         project%nampar(4) = 'd60'  ; project%value(4) = d60
         project%nampar(5) = 'nis'  ; project%value(5) = ni
         project%nampar(6) = 'njs'  ; project%value(6) = nj

         if (GRTYP.EQ.'N') then
            latproj = +90. ; nhem = 1.
         else if (GRTYP.EQ.'S') then
            latproj = -90. ; nhem = 2.
         endif

         project%nampar(7) = 'latproj' ; project%value(7) = latproj
         project%nampar(8) = 'nhem'    ; project%value(8) = nhem

      else &
      if (GRTYP.eq.'Z') then

         allocate(  alon( ni    ) , &
                    alat(    nj ) , &
                    lonr( ni*nj ) , &
                    latr( ni*nj ) )

         invj = .not.invj ! pas d'inversion par defaut
!
!        retrieve the x- and y-directional info

         call getzref( funit, '>>',alon )
         CALL GETZREF( funit, '^^',alat )

!        Retrieve associated rotation and pole info

         CALL GETZDES( ZTYP, ZIG1,ZIG2,ZIG3,ZIG4, ni,nj )

         if (ZTYP == 'E') then

!           Rotated pole grid.

            project%len  = 5
            project%name = 'rotated_pole'

            if (ni.ne.ibuf(5)  .or. &
                nj.ne.ibuf(6)) call                XIT('inq_file2',-12 )

            call CIGAXG( ZTYP, dlat1,dlon1,dlat2,dlon2, &
                               ZIG1, ZIG2, ZIG3, ZIG4 )

            project%nampar(1) = 'dlon1' ; project%value(1) = dlon1
            project%nampar(2) = 'dlat1' ; project%value(2) = dlat1
            project%nampar(3) = 'dlon2' ; project%value(3) = dlon2
            project%nampar(4) = 'dlat2' ; project%value(4) = dlat2

         else if (ZTYP == 'N' .or. ZTYP == 'S') then

!           Polaire-stereographique codee dans une grille Z

            project%len  = 8
            project%name = 'polar_stereographic'

            CALL CIGAXG( ZTYP, pi , pj , d60, dgrw , &
                               ZIG1,ZIG2,ZIG3,ZIG4  )

            dx =  alon(2)-alon(1) ; dy =  alat(2)-alat(1)
            pi = -(alon(1)-dx)/dx ; pj = -(alat(1)-dy)/dy

            if (dx /= dy) write(6,6004) dx*d60,dy*d60
            d60 = dx*d60

!           code/decode pour confirmer les valeurs des parametres

            CALL CXGAIG( ZTYP, ZIG1,ZIG2,ZIG3,ZIG4, &
                               pi,  pj,  d60, dgrw )
            CALL CIGAXG( ZTYP, pi , pj , d60, dgrw, &
                               ZIG1,ZIG2,ZIG3,ZIG4  )

            project%nampar(1) = 'pi'   ; project%value(1) = pi
            project%nampar(2) = 'pj'   ; project%value(2) = pj
            project%nampar(3) = 'dgrw' ; project%value(3) = dgrw
            project%nampar(4) = 'd60'  ; project%value(4) = d60
            project%nampar(5) = 'nis'  ; project%value(5) = ni
            project%nampar(6) = 'njs'  ; project%value(6) = nj

            if (ZTYP.EQ.'N') then
               latproj = +90. ; nhem = 1.
            else if (ZTYP.EQ.'S') then
               latproj = -90. ; nhem = 2.
            endif

            project%nampar(7) = 'latproj' ; project%value(7) = latproj
            project%nampar(8) = 'nhem'    ; project%value(8) = nhem

         else


!           Type de vecteurs de references non-supporte

            write(6,6033) ZTYP
            call                                   xit('lire_arg',-11 )

         endif

      else

!        Projection horizontale non-supportee

         write(6,6003) GRTYP
         call                                      xit('lire_arg',-11 )

      endif

      return
!-----------------------------------------------------------------------

 6001 format(/" A l'appel, definir : ",a/)
 6003 format(/' Projection horizontale: GRTYP=',A,' non supporte.'/)
 6033 format(/' Grille de type Z avec vecteurs de references de type ', &
            A,' non supporte.'/)
 6004 format(/' Resolutions polaire-stereographiques en x et y ='/ &
                1x,f15.5,1x,f15.5,' sont differentes.'/)
 6007 format(33x,a)
 6009 format(34x,a)
 6010 format(/'Mauvaise valeur pour -',a/a)

      end subroutine rpn_params
      Integer*8 function grosse_date ( ibuf )

      implicit none

      include 'cdf2ccc.h'

      integer, parameter :: head = taille_entete

      integer  ibuf(head)

      integer  dateo,datev, npas,deet, dtpr,tmpr, datchek

      integer, external :: gethigh,newdate

      if (tlbl .or. ccc_pktyp(1:2) == 'SQ') then

         DATEO = GETHIGH('DATEO', IBUF )
         NPAS  = GETHIGH('NPAS' , IBUF )
         DEET  = GETHIGH('DEET' , IBUF )

         call incdatr( DATEV,DATEO, NPAS*(DBLE( DEET )/3600.) )
         datchek = newdate( DATEV, dtpr,tmpr, -3 )
         if (datchek /= 0) call                 xit('grosse_date', -1 )

         grosse_date = int( dtpr,8 )* 100000000 + tmpr

      else

         grosse_date = ibuf(2)

      endif

      return
      end function grosse_date
