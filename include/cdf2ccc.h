!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!AUTEUR	  Guy Bergeron  mai  2003
!
! REVISIONS
!
! B. Dugas octobre 2016 :
! - Le parametre HEAD passe de 28 a 32.
! B. Dugas juin 2013 :
! - Ajouter la variable fill_toler (tolerance d''erreur) 
! - Ajouter la variable meta-title (titre inserable dans 
!   la liste des attributs globaux du fichier NetCDF)
! B. Dugas aout 2012 :
! - Ajouter infvar(:)%unique_L=T pour tenir compte du cas ou
!   nlev=1 et que les toutes les varibles partagent la meme valeur
! - Ajouter l''argument "-cell_method" definissant le type
!   de traitement temporel qu''on subit les donnees
! B. Dugas juillet 2012 :
! - attr%cvalue passe a 512 caracteres
! - Ajouter infvar(:)%time_bnds_L lorsque les time_bnds s''appliquent
! - Ajouter le support de 'Log Pressure Hybrid Levels' (VKIND=5002)
! - Ajouter la variable gribcode (no. de la table GRIB correspondant
!   a une grille Lambert conforme conique NCEP)
! B. Dugas juin 2012 :
! - Allouer 512 caracteres pour les noms de fichiers
! - max_vars (=300 au lieu de 1000), ajouter max_levs (=200)
! - Regrouper les indicateurs de coordonnees verticales dans /VCARD/
! - Ajouter les vecteurs pas et niv au type attribut (infvar)
! - Ajouter ce qu''il faut pour la prise en charge des TocToc (!!)
! - Ajouter les variables logiques GEM2 et GEM3 (hybride normalisee)
! - Ajouter noUD, controlant l''usage de UDUNITS dans les
!   calculs des deltas temporels
! B. Dugas mai 2012 :
! - Ajouter l''option de conversion des variables asujeties
!   a des Delta T accumules ==> time_bnds (via dtsize)
! B. Dugas decembre 2011 :
! - Ajouter une autre possibilite de niveaux verticaux.
! - Ajouter CCCVX pour identifier le calendrier 360 jours.
! - Ajouter xcoord,ycoord et zcoord, les noms des
!   coordonnees en x,y,z dans le fichier netcdf.
! B. Dugas mai 2009 :
! - On peut specifier les noms des coordonnees spaciales
!   et temporelles en arguments. Cette option est requise
!   lorsqu''un fichier NetCDF utilise un nom de coordonnee
!   non reconnu par le convertisseur
! B. Dugas oct 2008 :
! - nlvl == 12 plutot que 10
! - Ajouter l''option logique non_geographique
! - Re-ordonner les declarations des blocs commons
!   /attr_com/,/card/,/scan/ et du type attribut
! B. Dugas mai 2008 :
! - Ajouter possibilite%time(ntun)
! B. Dugas automne/hiver 2007 :
! - Re-formattage (enlever les TABs)
! - Ajouter l''element npack au type derive staggered
! - Les variables miss_ccc,fill_ccc,range passent a REAL*8
! - Nouvelles declarations pour attunit, direction et ccc_pktyp
! - La variable udunit_dat est remplace par udunits.def    
! - Les variables hyb_pt, hyb_r et hyb_pref sont ajoutees dans /card/
! - Ajouter la definition du parametre head de taille d''entete (=28)
! - Les parametres (nlvl=7,ngrd=3) deviennent (nlvl=10,ngrd=5)
! - De-commenter la declaration de "range" dans le type staggered
! - Enlever la declaration de NH4TO8 (variable non utilisee)
! B. Dugas    Aout 2007 : Modifications pour support des fichiers CMC/RPN
! Anne Frigon Aout 2006 : Ajoute integer cle_nhem dans common card
!                         associe a la nouvelle cle d appel de cdf2ccc.f
!                         introduite dans lire_arg.f
! Anne Frigon Aout 2006 : Modifie 
!                         parameter (udunits_dat="/LOGICIELS/freeware/etc/udunits.dat")
!                         pour adapter a la nouvelle structure commune a rhea et a cronos
!                         parameter (udunits_dat="/usr/local/etc/udunits.dat")
! C Desrochers Mar 2005: Ajoute real fill_ccc dans common card
!                        associe a la cle d appel de cdf2ccc.f
!                        et aussi le defaut fill_ccc_def
! Anne Frigon Oct 2003 : Ajoute real miss_ccc dans common card
!                        associe a la cle d appel de cdf2ccc.f
!                        et aussi le defaut miss_ccc_def
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

!Identification de la version

      character(len=80) :: version
      character(len=20) :: vdate
      common /version_info/ version,vdate

!Mode de fonctionnement ('cdf2cc' ou bien 'cdf2rpn')

      character(8) cdf2_mode
      common /operation_mode/ cdf2_mode

!Declaration des dimensions (netCDF<-->CCCma):

      integer, parameter :: max_vars  = 300  !nombre maximal de variables < 2000
      integer, parameter :: max_dims  = 20   !nombre maximal de dimension < 100  (ass. a une variable)
      integer, parameter :: max_attrs = 50   !nombre maximum d''attributs  < 2000 (ass. a une variable)
      integer, parameter :: max_levs  = 200  ! nombre maximal de niveaux verticaux

      real,    parameter :: vlarge    = 1.d+38 ! Utilise pour definir scale et offset

!
! NOTA: 
!
!  max_dims : uniquement utilisee pour dimensionner "dimid" dans les 
!             types derives "information" et "liste".
!
!  max_attrs: uniquement utilisee pour dimensionner le "type derive" attr.
!
!  max_vars : uniquement utilisee pour dimensionner le vecteur name dans inq_file.


!
!FICHIERS :

      integer(4)     attunit       ! 
      character(512) attr_file     ! fichier de definitions (attribut_netcdf.dat)
      character(512) cccma_file    ! fichier CCCma
      character(512) netcdf_file   ! fichier netCDF
      character(512) meta_title    ! titre inserable dans les attributs globaux
      character(8)   direction     ! vers 'netcdf', 'cccma' ou bien 'rpncmc'
      character(4)   ccc_pktyp     ! type du fichier CCCma ou RPN/CMC

      common /file_comm/ attr_file,cccma_file,netcdf_file,
     .                   meta_title,direction,attunit,ccc_pktyp

!
!UDUNITS2 :

      common /booting/ udunits_dat
      character(512) udunit_def,udunits_dat  ! fichier d''initialisation UDUNITS2

      parameter (udunit_def="/opt/udunits2/share/udunits/udunits2.xml")


!
!netCDF Nombres :

      integer(4) ndims       ! nombre de dimension dans le fichier netCDF
      integer(4) nvars       ! nbre de variables dans le fichier netCDF
      integer(4) ncoord      ! nbre de coordonnees dans le fichier netCDF
      integer(4) nlist       ! nbre de variables dans la liste

      common /number/ ndims,nvars,ncoord,nlist

!
!netCDF Attributes :

      integer, parameter :: max_len = 2  ! Longueur des valeurs associees aux attributs
	
      real(8) range(max_len)     ! min et max

      common /scan/ range,infvar,delta
!

      integer(4) nattr          ! nbre d''attributs associes a une variable

      TYPE attribut
       sequence
       real(8)        dvalue(max_len)   ! valeur de type double
       real(4)        rvalue(max_len)   ! valeur de type float
       integer(4)     ivalue(max_len)   ! valeur de type integer
       integer(4)     len               ! longueur de l''attribut
       integer(4)     type              ! type de l''attribut
       integer(2)     i2value(max_len)  ! valeur de type short
       integer(1)     i1value(max_len)  ! valeur de type byte
       character(512) cvalue            ! valeur de type character
       character(80)  name              ! nom de l''attribut
       character(386) padding           ! taille totale = 1024 octets
      END TYPE attribut
      TYPE (attribut) attr(max_attrs)    

      common /attr_com/ attr,nattr    

!
!netCDF spectral :

      logical(4)     spec ! Vrais si le fichier est en coeff. sepc.

      common /inter/ spec

!
!Staggered grid :
!
! NOTA : Utilise pour les fichiers CCCma et CMC/RPN a plusieurs variables
!
      TYPE staggered
       sequence
       real(8)      range(2)          ! min et max
       integer(8),  pointer :: pas(:) ! pas de temps
       integer(4),  pointer :: niv(:) ! niveaux verticaux
       integer(4),  pointer :: ip3(:), dateo(:)
       character(4) name              ! nom de la variable
       integer(4)   npack             ! facteur de compaction
       integer(4)   ndim              ! nbre de dimensions
       integer(4)   len(max_dims)     ! longueur de la dimension
       logical(4)   var_ok            ! denote une variable a traiter
       logical(4)   unique_L          ! denote une variable unique
       logical(4)   time_bnds_L       ! variable sujete aux "time_bnds"
       integer(4)   tid,zid           ! ID des dimensions temporelle et verticale
      END TYPE staggered
      TYPE (staggered) infvar(max_vars)

!
!
      TYPE freq
       sequence
       real(8)       dval     ! valeur de la frequence d''archivage
       character     type     ! type de frequence (hour,day,month)
       character(7)  padding  ! taille totale = 16 octets
      END TYPE freq
      TYPE (freq) delta
!
!Consistance de la definition des etiquettes de niveaux CCCma:
!
      logical(4) chklvl    ! utilisee pour preserver la consitance de l''etiquette
                           ! de niveau cccma lors d''un aller retour 
                           ! (i.e. lvcode lvdcode)
      common /lvl/ chklvl
!
!CCCma Cles d''entree du programe :
!

! Divers

      logical(4)    invj          ! Inverse l'ordre de l'indice 'j' dans la sortie.
      logical(4)    cccvx         ! utiliser un calendrier a 360 jours (T/F)
      logical(4)    leap          ! tenir compte des annee bissextile (true/false)
      logical(4)    noUD          ! ne pas utiliser UDUNITS pour les calculs temporels
      logical(4)    lalo          ! sortir les longitudes et latitudes (T/F)
      logical(4)    tlbl          ! etiquette temporel cccma = aaaammddhh (T/F)
      logical(4)    miss_ccc_def  ! DEFAUT de la valeur manquante pour sortie ccc
      logical(4)    miss_ccc_oui  ! Indique qu''au moins une valeur miss_ccc a ete remplacee
      logical(4)    fill_ccc_def  ! DEFAUT de la valeur de remplissage pour sortie ccc
      logical(4)    fill_ccc_oui  ! Indique qu''au moins une valeur fill_ccc a ete remplacee
      logical(4)    non_geographique ! Indique la conversion de variables non-geographiques
      logical(4)    time_bnds_L   ! Tenir compte des limites temporelles associees aux moyennes
      integer(4)    ladate        ! date de depart de la simulation
      integer(4)    npack         ! indice de compression
      real(4)       dt            ! pas de temps en sec
      real(4)       tmoyen        ! temperature utilisee pour definir les ibuf
      real(8)       dtsize        ! Interval d''accumulation en heures
      real(8)       miss_ccc      ! valeur manquante pour sortie ccc
      real(8)       fill_ccc      ! valeur de remplissage pour sortie ccc
      real(8)       fill_toler    ! tolerance d''erreur lors des calculs de valeurs manquantes
      character(80) grid_desc     ! identification de la grille/projection
      character(80) level_desc    ! identification des niveaux verticaux
      character(80) time_desc     ! identification des unites temporelles
      integer(4)    cle_nhem      ! hemisphere pour dir=cccma et grille PS
      integer(4)    gribcode      ! table GRIB d''une Lambert conforme conique NCEP
      character(128) xcoord       ! nom de la coordonnee x dans le fichier netcdf
      character(128) ycoord       ! nom de la coordonnee y dans le fichier netcdf
      character(128) zcoord       ! nom de la coordonnee z dans le fichier netcdf
      character(128) tcoord       ! nom de la coordonnee t dans le fichier netcdf
      character(128) cell_method  ! methode utilisee avec les time_bnds

      common /card/ fill_ccc,fill_toler,miss_ccc,dtsize,fill_ccc_def,
     .              miss_ccc_def,fill_ccc_oui,miss_ccc_oui,invj,cccvx,
     .              leap,noUD,lalo,tlbl,ladate,npack,dt,tmoyen,
     .              grid_desc,level_desc,cle_nhem,gribcode,
     .              time_desc,tcoord,non_geographique,
     .              time_bnds_L,xcoord,ycoord,zcoord,
     .              cell_method

      logical(4)    gem3,gem2      ! Identifie les versions normalisees ou pas de la
      ! coordonnee GEM ETA codee en mode sigma (avec enregistrement de support 'HY')

      real(4)       htoit          ! hauteur du toit du model en metres
      integer(4)    vkind,nkm,nkt  ! \
      real(4)       hyb_pt,hyb_r2  ! parametres definissant les differentes
      real(4)       hyb_pref,hyb_r ! saveurs de la coordonnee hybride utilisee
      real(8)       cap (max_levs) ! dans les fichiers CMC/RPN
      real(8)       cb  (max_levs) ! /

      common /vcard/ cap,cb,gem2,gem3,htoit,vkind,nkm,nkt,
     .               hyb_pt,hyb_pref,hyb_r,hyb_r2

! Entete I/O :
      integer, parameter :: head = 32 ! taille de l''entete de ibuf

! Projection :

      integer, parameter :: nprj = 20

      TYPE projection
       sequence
       character(80) name         ! nom de la projection
       character(80) nampar(nprj) ! nom du parametre de projection
       integer(4)    len          ! nbre de parametres de projection
       real(4)       value(nprj)  ! valeur du parametre de projection
       character(72) padding      ! taille totale = 240 octets
      END TYPE projection
      TYPE (projection) project

      common /projects/ project

! Les possibilites

      integer, parameter :: nlvl = 14   !nbre de cas possible pour level_desc
      integer, parameter :: ngrd = 6    !nbre de cas possible pour grid_desc
      integer, parameter :: ntun = 6    !nbre de cas possible pour time_units

      TYPE possibilite
       sequence
       character(80) level(nlvl)   ! liste des level_desc possibles
       character(80) grid(ngrd)    ! liste des grid_desc possibles
       character(80) time(ntun)    ! liste des time_desc possibles
      END TYPE possibilite
      TYPE (possibilite) possible

      common /possibles/ possible


!
!UTIL Allocation de memoire

      integer(4) maxdim      ! nbr de dimensions 
      integer(4) maxvar      ! nbr de variables
      integer(4) maxlen      ! longueur maximale 3D
      integer(4) max1d       ! longueur maximale 1D
      integer(4) maxlev      ! nombre de niveaux
      integer(4) maxtime     ! nombre d''enregistrement

      integer(4) maxpk       ! longueur de ibuf

      common /mem_dim/ maxdim,maxvar,maxlen,max1d,maxlev,maxtime,maxpk
!

! 
!*
