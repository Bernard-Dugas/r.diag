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
      subroutine lire_arg ( GCLES,GDEF2,GDEF1,NBRGEN,
     +                      IONAM,IPOS,NBRUNT )

      IMPLICIT        none

***    Arguments Entree/Sortie

      INTEGER         NBRGEN,IPOS,NBRUNT
      CHARACTER(512)  GDEF1(NBRGEN),IONAM(NBRUNT+1),
     +                GDEF2(NBRGEN)
      CHARACTER(16)   GCLES(NBRGEN)

      include 'cdf2ccc.h'
      include 'ztypmem.h'

******
*
*AUTEUR Guy Bergeron     juilllet 2003
*
*      Lits les NBRGEN parametres generaux and lrd NCLE locaux.
*      Au retour, IONAM(1) contient le nom du programme et les
*      noms des fichiers qui devront etre ouverts par JCLPNT
*      sont  dans les IPOS-1 elements suivants de IONAM
*
*REVISIONS
*
*  B. Dugas avril '17 :
*   - Utiliser get_environment_variable pour verifier la valeur de
*     la variables d'environnement UDUNITS2_XML_PATH. Celle-ci a
*     pre-seance sur udunits2_def
*  B. Dugas fevrier '17 :
*   - Remplacer UDUNITS_PATH par UDUNITS2_XML_PATH
*   - Remplacer udunits.dat par udunits2.xml
*  B. Dugas juin '13 :
*   - Enlever la limitation de juste definir des valeurs manquantes
*     lorsqu'on lit des fichier NetCDF. L'operation en mode inverse
*     (ecriture) est maintenant supportee
*   - Modifier le traitement des arguments -miss_ccc et -fill_ccc
*     en ce sens que miss_ccc est declare desuet
*   - Definir la variable fill_toler qui sera utilisee lors des
*     comparaison des donnees avec fill_ccc
*   - Lire l'argument "-title" destine a la section des
*     attributs globaux du fichier NetCDF (def = ' ')
*
*  B. Dugas octobre 2012 :
*   - Ajouter la possibilite de lire les parametres de la
*     coordonnee hybride meme pour le cas cdf2_mode=cdf2rpn
*
*  B.Dugas aout '12 :
*   - Ajouter l'argument caractere -cell_method definissant
*     le type de traitement temporel qu'on subit les donnees
*   - Permettre la definition de la variable level_desc
*     (via l'argument -lev) pour un fichier source NetCDF
*
*  Bernard Dugas, juillet 2012 :
*   - La valeur par defaut de -rlonoff est de -1000. Celle-ci
*     sera traduite de facon appropriee par les differentes
*     routines selon le type de grilles rencontrees
*   - Ajouter l'argument -gribcode (pour GRIB/NCEP Lambert Conform)
*   - Mode HELP actif si aucun de nom de fichier n'est specifie
*
*  Bernard Dugas, juin 2012 :
*   - Definir les variables logiques GEM2 et GEM3 (hybride normalisee)
*   - Faire la distinction entre les calendriers 'gregorian'/
*     'standard' et 'proleptic_gregorian': On ne peut utiliser
*     les fonctions basees sur UDUNITS avec 'proleptic_gregorian'
*   - Assurer la coherence du contenu des arguments -leap et -calendar
*
*  B. Dugas: mai 2012
*   - Ajouter l'argument '-dtsize', i.e. un interval d'accumulation
*     en heures qui sera utilise lors des calculs des 'time_bnds'
*     mais seulement en mode 'cdf2rpn'
*   - Traiter les arguments optionels '-dt' et '-dateo' dans la
*     direction 'rpncmc' afin de pouvoir definir DATEO,DEET et NPAS,
*     tout en tenant comptes des 'time_bnds'
*   - Definition du postfixe etendue dans la routine lprognam
*   - Activer le mode '-help'. Ignorer alors tous les autres arguments
*
*  B. Dugas: avr 2012
*   Ajouter le calendrier 'proleptic_gregorian' (i.e. 'gregorian')
*
* B. Dugas decembre 2011 :
*   Ajouter la cle '-calendar' ('gregorian', '365_day' ou '360_day')
*
* B. Dugas juillet 2011 :
*   La cle '-leap' est aussi consideree dans le mode cdf2rpn
*
* B. Dugas mai 2009 :
*   Ajouter les arguments -[txyz]coord pour les noms des
*   coordonnees a chercher dans un fichier NetCDF
*
* B. Dugas octobre 2008 :
*   Ajouter l'argument -nongeog (convertir variables non-geographiques)
*
* B. Dugas juillet 2008 :
*   La variable npack est convertie en unites de bits/mots
*
* B. Dugas mai 2008 :
*   La variable ladate est conservee en format date-time-stamp
*   Ajouter l'argument -timdesc (unites temporelles)
*
* B. Dugas ete/automne/hiver 2007:
*   Valeur par defaut de npack = 999 (on utilise alors nf_type)
*   Ajout des cles 'rlonoff', 'hyb_pt', 'hyb_pref', 'hyb_r' et 'phis'
*
* B. Dugas avril 2007:
*   La routine recoit et retourne un jeux de cles supplementaires qui
*   seront traitees par les_arg. On suppose qu'elle est maintenant
*   appellee par jclpnt. Cette derniere routine s'occupe d'ouvrir
*   les unites d'entree/sortie CCCma et/ou RPN/CMC (dont les noms
*   sont retournes dans IONAM) et le fichier attribut_netcdf.dat
*   Le nom (et le mode de fonctionnement) sont determines par
*   un appel a la routine LPROGNAM (dont le source est inclus)
*   On ajoute enfin les cles 'udunits' et 'rpn'
*
* A. Frigon aout 2006:
*   Ajoute 28e cle cle_nhem (valeur de defaut invalide=99) pour generaliser grille
*   polar_stereographic sur hemispheres NORD(1) SUD(2)
*   si direction=cccma
*   (car avant cdf2ccc 1.5 on assignait nhem=1 dans rdlatlon2.f).
*   Cette nouvelle variable cle_nhem a ete introduite dans cdf2ccc.h
*   et ajoutee dans common card.
*   Tout ceci parce que trier.f neglige d'identifier
*   la variable polar_stereographic car definie de type char.
*
* A. Frigon juin 2006 : 
*   Corrige description d60 "vraie a 60N" par "vraie a 60 deg"
*   car vers netcdf tout est general pour PS nord/sud 
*   selon IBUF(7) lu et assigne a nhem
*   dans def_dim_coord.f
*
* G. Bergeron juin 2005:
*   -Introduction de cles pour grille lon/lat regionale et le type derive "project".
*   -Cle "-invj" pour controler l'inversion des indices "j".
*
* A. Frigon Avril 2005:
*    Remplace DATA par assignation correcte
*    pour valeurs par defaut de miss_ccc_def et miss_ccc_def 
*
* C. Desrochers Mar 2005:
*    Ajoute 22e cle fill_ccc facultative pour assigner valeurs 
*    de remplissage dans sortie ccc. Valeur logique fill_ccc_def 
*    pour detecter utilisation cle. Par defaut : fill_ccc_def=.false.
*    et on ne modifie pas les valeurs de remplissage dans la conversion
*    cdf vers ccc.

* Anne Frigon Oct 2003 : 
*    Ajoute 17e cle miss_ccc facultative pour assigner valeurs 
*    manquantes dans sortie ccc. Valeur logique miss_ccc_def 
*    pour detecter utilisation cle Par defaut : miss_ccc_def=.false.
*    On ne modifie pas les valeurs manquantes dans la conversion
*    cdf vers ccc.
*
* Anne Frigon Oct 2003: 
*    Defini file_attr,local comme ch*128 au lieu de ch*80* 
*
*******

      character(512) file_attr,local
      parameter(local='attribut_netcdf.dat')
      parameter(file_attr=
     .'/LOGICIELS/cdf2ccc/etc/attribut_netcdf.dat')

******les_arg(ccard)

      integer, parameter :: ncle = 46 ! nombre de cles

      character*16  cles(ncle)    ! nom de la cle
      character*60  def(ncle)     ! defenition de la cle
      character*512 def1(ncle)    ! sortie et defaut si la cle n'est pas la
      character*512 def2(ncle)    ! defaut si juste la cle est la
      character*512 evalue
      
*     Les NCLE parametres locaux affectent les variables 
*     correspondantes dans cdf2ccc.h

      data
     .  cles(1) /'cdf'    /,  def1(1) /'?'       /,  def2(1) /'?'     /,
     .  cles(2) /'ccc'    /,  def1(2) /'?'       /,  def2(2) /'?'     /,
     .  cles(3) /'dir'    /,  def1(3) /'def'     /,  def2(3) /'netcdf'/,
     .  cles(4) /'leap'   /,  def1(4) /'?'       /,  def2(4) /'no'    /,
     .  cles(5) /'dateo'  /,  def1(5) /'?'       /,  def2(5) /'0'     /,
     .  cles(6) /'dt'     /,  def1(6) /'?'       /,  def2(6) /'0.0'   /,
     .  cles(7) /'tlbl'   /,  def1(7) /'no'      /,  def2(7) /'yes'   /,
     .  cles(8) /'lev'    /,  def1(8) /'?'       /,  def2(8) /'?'     /,
     .  cles(9) /'tm'     /,  def1(9) /'220.0'   /,  def2(9) /'?'     /,
     .  cles(10)/'ht'     /,  def1(10)/'?'       /,  def2(10)/'?'     /,
     .  cles(11)/'grid'   /,  def1(11)/'?'       /,  def2(11)/'?'     /,
     .  cles(12)/'ni'     /,  def1(12)/'?'       /,  def2(12)/'?'     /,
     .  cles(13)/'nj'     /,  def1(13)/'?'       /,  def2(13)/'?'     /,
     .  cles(14)/'pi'     /,  def1(14)/'?'       /,  def2(14)/'?'     /,
     .  cles(15)/'pj'     /,  def1(15)/'?'       /,  def2(15)/'?'     /,
     .  cles(16)/'dgrw'   /,  def1(16)/'?'       /,  def2(16)/'?'     /,
     .  cles(17)/'d60'    /,  def1(17)/'?'       /,  def2(17)/'?'     /,
     .  cles(18)/'0lon'   /,  def1(18)/'GLOBAL'  /,  def2(18)/'?'     /,
     .  cles(19)/'0lat'   /,  def1(19)/'?'       /,  def2(19)/'?'     /,
     .  cles(20)/'dlon'   /,  def1(20)/'?'       /,  def2(20)/'?'     /,
     .  cles(21)/'dlat'   /,  def1(21)/'?'       /,  def2(21)/'?'     /,
     .  cles(22)/'invj'   /,  def1(22)/'yes'     /,  def2(22)/'no'    /
     .  cles(23)/'npack'  /,  def1(23)/'999'     /,  def2(23)/'?'     /,
     .  cles(24)/'lalo'   /,  def1(24)/'no'      /,  def2(24)/'yes'   /,
     .  cles(25)/'attr'   /,  def1(25)/file_attr /,  def2(25)/local   /,
     .  cles(26)/'miss_ccc'/, def1(26)/'?'       /,  def2(26)/'ERR'   /,
     .  cles(27)/'fill_ccc'/, def1(27)/'?'       /,  def2(27)/'ERR'   /,
     .  cles(28)/'cle_nhem'/, def1(28)/'?'       /,  def2(28)/'?'     /,
     .  cles(29)/'udunits'/,  def1(29)/'default' /,  def2(29)/'?'     /,
     .  cles(30)/'rlonoff'/,  def1(30)/'?'       /,  def2(30)/'?'     /,
     .  cles(31)/'hyb_pt' /,  def1(31)/'?'       /,  def2(31)/'?'     /,
     .  cles(32)/'hyb_pref'/, def1(32)/'?'       /,  def2(32)/'?'     /,
     .  cles(33)/'hyb_r'  /,  def1(33)/'?'       /,  def2(33)/'?'     /,
     .  cles(34)/'rpn'    /,  def1(34)/'?'       /,  def2(34)/'?'     /,
     .  cles(35)/'phis'   /,  def1(35)/'?'       /,  def2(35)/'?'     /,
     .  cles(36)/'timdesc'/,  def1(36)/'hours'   /,  def2(36)/'?'     /, 
     .  cles(37)/'nongeog'/,  def1(37)/'oui'     /,  def2(37)/'non'   /,
     .  cles(38)/'xcoord' /,  def1(38)/'!@#$%^&' /, def2(38)/'!@#$%^&'/,
     .  cles(39)/'ycoord' /,  def1(39)/'!@#$%^&' /, def2(39)/'!@#$%^&'/,
     .  cles(40)/'zcoord' /,  def1(40)/'!@#$%^&' /, def2(40)/'!@#$%^&'/,
     .  cles(41)/'tcoord' /,  def1(41)/'!@#$%^&' /, def2(41)/'!@#$%^&'/,
     .  cles(42)/'dtsize' /,  def1(42)/'0.0'     /,  def2(42)/'?'     /,
     .  cles(43)/'calendar'/, def1(43)/'?'       /,  def2(43)/'?'     /,
     .  cles(44)/'gribcode'/, def1(44)/'?'       /,  def2(44)/'?'     /,
     .  cles(45)/'cell_method'/, def1(45)/'?'    /,  def2(45)/'?'     /,
     .  cles(46)/'title'  /,  def1(46)/' '       /,  def2(46)/' '     /

      data
     .  def(1) /'(C) Nom du fichier netCDF'                           /, 
     .  def(2) /'(C) Nom du fichier CCCma'                            /, 
     .  def(3) /'(C) Convertir de/vers (cccma ou rpncmc) a/de netcdf)'/,
     .  def(4) /'(C) Annee bissextile (no / yes)'                     /,
     .  def(5) /'(I) Date de depart de la simulation (AAAAMMJJHH)  '  /, 
     .  def(6) /'(R) Pas de temps (secondes)  '                       /, 
     .  def(7) /"(I) L'etiquette temporel cccma = AAAAMMDDHH"         /,
     .  def(8) /'(C) Type de niveaux'                                 /, 
     .  def(9) /'(R) variable TMOYEN de PARAMETRES  '                 /,
     .  def(10)/'(R) Hauteur du toit du modele en metres'             /,
     .  def(11)/'(C) Type de projection'                              /, 
     .  def(12)/'(I) Nbre de points de grille en X (type f)'          /,
     .  def(13)/'(I) Nbre de points de grille en Y (type f)'          /,
     .  def(14)/'(R) Coordonnee selon x du pole (nbr de dx)'          /,
     .  def(15)/'(R) Coordonnee selon y du pole (nbr de dy)'          /,
     .  def(16)/"(R) Angle entre Greenwich et l'axe X (deg. ouest)"   /,
     .  def(17)/'(R) Longueur de la maille vraie a 60 deg (metres)'   /,
     .  def(18)/"(R) Longitude d'origine (degres ouest)"              /,
     .  def(19)/"(R) Latitude d'origine (degres nord)"                /,
     .  def(20)/"(R) Longueur de la maille selon longitude (degres)"  /,
     .  def(21)/"(R) Longueur de la maille selon latitude (degres)"   /,
     .  def(22)/"(C) Inverse l'ordre de l'indice 'j' dans la sortie"  /,
     .  def(23)/'(I) Densite de compression 0,1,2,4,-64,-32,-16'      /,
     .  def(24)/"(C) Sortir les latitudes et longitudes (no / yes)"   /,
     .  def(25)/'(C) Le fichier attribut_netcdf.dat'                  /,
     .  def(26)/'(R) Argument desuet. Utiliser plutot fill_ccc'       /,
     .  def(27)/'(R) Valeur de remplissage dans fichier CCCma'        /,
     .  def(28)/'(I) Hemisphere 1nord 2sud pour grille PS vers CCC'   /,
     .  def(29)/'(C) Chemin complet du fichier udunits2.xml'          /,
     .  def(30)/"(R) Deplacement des longitudes d'une grille tournee" /,
     .  def(31)/'(R) Pression au toit de la coordonnee hybride (Pa)'  /,
     .  def(32)/'(R) Pression de reference pour la coordonnee hybride'/,
     .  def(33)/'(R) Exposant pour la coordonnee hybride'             /,
     .  def(34)/'(C) Nom du fichier RPN/CMC'                          /,
     .  def(35)/'(C) Nom du fichier PHIS (Option Gal-Chen)'           /,
     .  def(36)/'(C) Unites associes a la variable temporelle'        /,
     .  def(37)/'(C) Convertir les variables non-geographiques'       /,
     .  def(38)/'(C) Nom de la coordoonnee en X du fichier NetCDF'    /,
     .  def(39)/'(C) Nom de la coordoonnee en Y du fichier NetCDF'    /,
     .  def(40)/'(C) Nom de la coordoonnee en Z du fichier NetCDF'    /,
     .  def(41)/'(C) Nom de la coordoonnee en T du fichier NetCDF'    /,
     .  def(42)/"(R) Interval d'accumulation en heures"               /,
     .  def(43)/'(C) Nom du calendrier (gregorian, 365_day, 360-day)' /,
     .  def(44)/'(I) Code GRIB pour une grille Lambert conforme conic'/,
     .  def(45)/'(C) Cell Method utilisee dans les calculs temporels' /,
     .  def(46)/'(C) Optional "title" meta-data'                      /

******

      integer i,nlen,nis,njs,idate
      logical ok

***    Champs de travail locaux pour les_arg

      CHARACTER(4)    DEF_PKTYP
      CHARACTER(16),  DIMENSION(:), ALLOCATABLE :: ACLES
      CHARACTER(512), DIMENSION(:), ALLOCATABLE :: ADEF,ADEF1,ADEF2

      CHARACTER(512)  :: UDUNITS2_DEF

      INTEGER         NBRCLE,NDATE,NHELP,PART1,PART2

      INTEGER         NEWDATE,datchek,L_argenv
      EXTERNAL        NEWDATE

      COMMON         /ZZDEFPK/ DEF_PKTYP

      LOGICAL         IS_ON,IS_OFF
      EXTERNAL        IS_ON,IS_OFF

*-----------------------------------------------------------------------
      udunits2_def = trim( udunits2_def1 ) //
     .               trim( udunits2_def2 ) //
     .               trim( udunits2_def3 )

***    Allocate LES_ARG work fields.

      NBRCLE = NBRGEN + NCLE
      ALLOCATE ( ACLES(NBRCLE),ADEF (NBRCLE),
     .           ADEF1(NBRCLE),ADEF2(NBRCLE) )

      call define_possibilite()

***    Copier les NCLE parametres locaux

      DO  I=1,NCLE
          ACLES(I) = CLES(I)
          ADEF1(I) = DEF1(I)
          ADEF2(I) = DEF2(I)
          ADEF (I) = DEF (I)
      END DO

***    Copier les NBRGEN parametres generaux

      NDATE = -1
      NHELP = -1

      DO  I=1,NBRGEN
          ACLES(NCLE+I) = GCLES(I)
          CALL up2low( ACLES(NCLE+I),acles(ncle+i) )
          ADEF1(NCLE+I) = GDEF1(I)
          ADEF2(NCLE+I) = GDEF2(I)
          ADEF (NCLE+I) = 'Argument general pris en charge par JCLPNT'
          if (GCLES(I).eq.'DATE') NDATE = I
          if (GCLES(I).eq.'HELP') NHELP = I
      END DO

*     Valeurs par defaut :

      miss_ccc_def=.false.
      fill_ccc_def=.false.

*     Lire les parametres d'entrees :

      call les_arg(acles,adef,adef1,adef2,NBRCLE,version)

***    Sauver les NCLE parametres locaux

      DO  I=1,NCLE
          DEF1(I) = ADEF1(I)
      END DO

***    Sauver les NBRGEN parametres generaux

      DO  I=1,NBRGEN
          GDEF1(I) = ADEF1(NCLE+I)
      END DO

***    IPOS(1) contient le nom du programme executant

      IPOS=1
      call lprognam( IONAM(IPOS) )

*     Le mode de fonctionnement est define par le nom
*     du programme appellant, cdf2cc ou bien cdf2rpn.

      cdf2_mode = 'cdf2rpn'
      if (index( IONAM(IPOS),'cdf2rpn' ).eq.0 .and.
     .    index( IONAM(IPOS),'cdf2ccc' ).ne.0)
     .   cdf2_mode = 'cdf2ccc'
         
      if (cdf2_mode.eq.'cdf2ccc') DEF_PKTYP = ' '
      if (cdf2_mode.eq.'cdf2rpn') DEF_PKTYP = 'SQ98'

      if (cdf2_mode.eq.'cdf2ccc') then
         if (NDATE.gt.0) then
            if (GDEF1(NDATE).eq.' ') 
     +         GDEF1(NDATE) = '0'
         endif
      endif

*     Si aucun nom de fichier n'est specifie, mode HELP actif.

      if(def1(1)  == '?' .and. 
     .   def1(2)  == '?' .and. 
     .   def1(34) == '?' ) GDEF1(NHELP) = 'OUI'

      if (GDEF1(NHELP) == 'OUI') return

*     Affecter les variables correspondantes :

      if(def1(1).eq.'?') then
         write(6,6001) ' -cdf "nom du fichier netCDF" ?'
         call                                       xit('lire_arg',  -1)
      else
         netcdf_file = def1(1)
      endif

      if(def1(2).eq.'?' .and. def1(34).eq.'?') then
         if(cdf2_mode.eq.'cdf2ccc') then
            write(6,6001) ' -ccc "nom du fichier CCCma ou RPN/CNMC" ?'
         else
            write(6,6001) ' -rpn "nom du fichier CCCma ou RPN/CNMC" ?'
         endif
         call                                       xit('lire_arg',  -2)
      elseif (def1(2).ne.'?' .and. def1(34).ne.'?') then
         if(cdf2_mode.eq.'cdf2ccc') then
            cccma_file = def1(2)
            I = index(cccma_file,'/')
            IPOS= IPOS+1
            IONAM(IPOS) = def1(2)
            if (i.eq.0) IONAM(IPOS) = './' // IONAM(IPOS)
         else
            cccma_file = def1(34)
            IPOS= IPOS+1
            IONAM(IPOS) = def1(34)
         endif
      elseif (def1(2).ne.'?') then
         cccma_file = def1(2)
         I = index(cccma_file,'/')
         IPOS= IPOS+1
         IONAM(IPOS) = def1(2)
         if (i.eq.0) IONAM(IPOS) = './' // IONAM(IPOS)
      else
         cccma_file = def1(34)
         IPOS= IPOS+1
         IONAM(IPOS) = def1(34)
      endif

      if (def1(3).eq.'def') then
         if(cdf2_mode.eq.'cdf2ccc') then
            def1(3)='cccma'
         else
            def1(3)='rpncmc'
         endif
      endif
      if( def1(3).eq.'netcdf' .or.
     .    def1(3).eq.'rpncmc' .or.
     .    def1(3).eq.'cccma') then
         direction = def1(3)
      else
         write(6,6001)' -dir est mal definie '
         call                                       xit('lire_arg',  -3)
      endif

*      Check for leap-year support.

      call getenvc( 'LEAP_YEAR_CONTROL',EVALUE )

      if (EVALUE /= ' ') then
         if (IS_OFF( EVALUE )) then
            def1(4)='off'
         else if (IS_ON( EVALUE )) then
            def1(4)='on'
         end if
      endif

      EVALUE=' '
      if (def1(4) == '?') then
         leap=.true.
      else if (IS_OFF( def1(4) )) then
         leap=.false.
         EVALUE='off'
      else if (IS_ON( def1(4) )) then
         leap=.true.
         EVALUE='on'
      else
         write(6,6001) ' -leap "no or yes" ?'
         call                                       xit('lire_arg',  -4)
      endif
      if (EVALUE /= ' ') def1(4)=EVALUE


      if(direction.eq.'netcdf') then
         if (len_trim( def1(5) ).ne.10) then
            if (cdf2_mode.eq.'cdf2rpn') then
               ladate = -1
            else
               write(6,6001)' -dateo "yyyymmddhh" '//trim( def1(5) )
               call                                 xit('lire_arg',  -5)
            endif
         else

            read(def1(5),9003) ladate
            PART1 = ladate / 100
            PART2 = ladate - PART1*100
            PART2 = PART2  * 1 000 000

***         convertir ladate en DATE-TIME-STAMP

            datchek = newdate( ladate, PART1,PART2, +3 )
            if (datchek.ne.0) then
               datchek = newdate( ladate, PART1,PART2, +5 )
               if (datchek.ne.0) call               xit('lire_arg',  -5)
            endif

         endif

      else if (cdf2_mode           == 'cdf2rpn' .and.
     .         len_trim( def1(5) ) ==  10)      then
         
         read(def1(5),9003) ladate
         PART1 = ladate / 100
         PART2 = ladate - PART1*100
         PART2 = PART2  * 1 000 000

***      convertir ladate en DATE-TIME-STAMP

         datchek = newdate( ladate, PART1,PART2, +3 )
         if (datchek.ne.0) then
            datchek = newdate( ladate, PART1,PART2, +5 )
            if (datchek.ne.0) call                  xit('lire_arg',  -5)
         endif

      else
         ladate = -1
      endif

      if (IS_OFF( def1(7) )) then
         tlbl=.false.
      else if (IS_ON( def1(7) )) then
         tlbl=.true.
      else
         write(6,6001) ' -tlbl "no or yes" ?'
         call                                       xit('lire_arg',  -7)
      endif

      if(.not.tlbl .and. direction.eq.'netcdf' ) then
         if(def1(6).ne.'?') then
            read(def1(6),9004,err=1000) dt
            goto 150
  100       dt = -1.0
  150       if(dt.le.0.0) then
               write(6,6001)' -dt "pas de temps" '//trim( def1(6) )
               call                                 xit('lire_arg',  -6)
            endif
         else if (cdf2_mode.eq.'cdf2rpn') then
            dt = -1.0
         else
            write(6,6001)' -dt "pas de temps" ? '
            call                                    xit('lire_arg',  -6)
         endif
      elseif(direction == 'rpncmc') then
         dt = -1.0
         if(def1(6) /= '?') then
            read(def1(6),9004,err=200) dt
            goto 250
  200       dt = -1.0
  250       if(dt < 0.0) then
               write(6,6001)' -dt "pas de temps" '//trim( def1(6) )
               call                                 xit('lire_arg',  -6)
            endif
         endif
      endif

      level_desc = ' ' ; gem2=.false. ; gem3 = .false.
      if (def1(8).ne.'?') then
         level_desc= def1(8)
         if (direction.eq.'netcdf') then
            evalue = level_desc ; call up2low( evalue,evalue )
            if(evalue == 'gem2' .or. evalue == 'gem3') then
               level_desc='Hybrid Levels' 
               if (evalue == 'gem2') gem2 =.true.
               if (evalue == 'gem3') gem3 =.true.
            endif
         endif
         ok=.false.
         do i=1,nlvl
            if(level_desc.eq.possible%level(i)) then
               ok=.true. ; exit
            endif
         enddo
         if(.not.ok) then
            write(6,6010) cles(8)//' = '//trim( def1(8) ),
     +                  ' Valeurs possibles de level_desc :'
            do i=1,nlvl
               write(6,6009) trim(possible%level(i))
            enddo
            call                                    xit('lire_arg',  -8)
         endif
      else if (direction.eq.'netcdf') then
         if (cdf2_mode.eq.'cdf2rpn') then
            level_desc = ' '
         else
            write(6,6001) ' -lev "level_desc" ?'
            call                                    xit('lire_arg',  -8)
         endif
      endif

      if(def1(9).eq.'?' ) then
         write(6,6001) ' -tm "TMOYEN" ?'
         call                                       xit('lire_arg',  -9)
      else
         read(def1(9),9004,err=1000) tmoyen
      endif

         IPOS        = IPOS+1
      if(def1(35).ne.'?') then
         IONAM(IPOS) = def1(35)
         if(cdf2_mode.eq.'cdf2ccc') then
            I = index(IONAM(IPOS),'/')
            if (i.eq.0) IONAM(IPOS) = './' // IONAM(IPOS)
         endif
      else
         IONAM(IPOS) =  '+'
      endif

      if(level_desc.eq.'Gal-Chen Levels')then

         if(def1(10).eq.'?') then
            write(6,6001) ' -ht "HTOIT" ?'
            call                                   xit('lire_arg',  -10)
         else
            read(def1(10),9004,err=1000) htoit
         endif

         if(def1(35).eq.'?') then
            write(6,6001) '-phis "PHIS FILE NAME" ?'
            call                                   xit('lire_arg',  -35)
         endif
         
      endif

      if(def1(11).ne.'?' .and. direction.eq.'netcdf') then
         if(.false.)grid_desc = def1(11)             !debug
         project%name = def1(11)
         ok=.false.
         do i=1,ngrd
            if(project%name.eq.possible%grid(i)) then
               ok=.true. ; exit
            endif
         enddo
         if(.not.ok) then
            write(6,6010) cles(11)//' = '//trim( def1(11) ),
     +                  ' Valeurs possibles de grid_desc :'
            do i=1,ngrd
               write(6,6007)trim(possible%grid(i))
            enddo
            call                                   xit('lire_arg',  -11)
         endif
      else if (direction.eq.'netcdf') then
         if (cdf2_mode.eq.'cdf2rpn') then
            project%name = ' '
         else
            write(6,6001) ' -grid "grid_desc" ?'
            call                                   xit('lire_arg',  -11)
         endif
      endif

      if(project%name.eq.'polar-stereographic')then

         project%name='polar_stereographic'
         project%len=7
         
         if(def1(12).eq.'?') then
            write(6,6001) ' -ni "NI" '
            call                                   xit('lire_arg',  -12)
         else
            read(def1(12),9003) nis
            project%nampar(5) = 'nis'
            project%value(5) = float(nis)
         endif

         if(def1(13).eq.'?') then
            write(6,6001) ' -nj "NJ" '
            call                                   xit('lire_arg',  -13)
         else  
            read(def1(13),9003) njs
            project%nampar(6) = 'njs'
            project%value(6) = float(njs)
         endif

         if(def1(14).eq.'?' )then
            write(6,6001) ' -pi "PI" ?'
            call                                   xit('lire_arg',  -14)
         else
            project%nampar(1) = 'pi'
            read(def1(14),9004,err=1000) project%value(1)
         endif

         if(def1(15).eq.'?' )then
            write(6,6001) ' -pj "PJ" ?'
            call                                   xit('lire_arg',  -15)
         else
            project%nampar(2) = 'pj'
            read(def1(15),9004,err=1000) project%value(2)
         endif
         
         if(def1(16).eq.'?' )then
            write(6,6001) ' -dgrw "DGRW" ?'
            call                                   xit('lire_arg',  -16)
         else
            project%nampar(3) = 'dgrw'
            read(def1(16),9004,err=1000) project%value(3)
         endif
         
         if(def1(17).eq.'?' )then
            write(6,6001) ' -d60 "D60" ?'
            call                                   xit('lire_arg',  -17)
         else
            project%nampar(4) = 'd60'
            read(def1(17),9004,err=1000) project%value(4)
         endif

      endif
*** FIN DE      if(project%name.eq.'polar_stereographic')then

      if(def1(23).eq.'?') then
         write(6,6001) ' -npack "NPACK" ?'
         call                                      xit('lire_arg',  -23)
      else

         read(def1(23),9003) npack

         if (npack.ne.999) then
            if (npack < -64 .or. npack > 64) then
               write(6,'(A,I4)') ' Valeur illegale pour npack ',npack
               call                                xit('lire_arg',  -23)
            else if (npack == 0) then
               npack = -64
            else if (npack > 0) then
               npack = -(64/npack)
            endif
         endif

      endif


      if(def1(25).eq.'?') then
         write(6,6001) ' Fichier attribut_netcdf.dat ?'
         call                                      xit('lire_arg',  -25)
      else
         attr_file = def1(25)
         IPOS=IPOS+1
         IONAM(IPOS) = def1(25)
      endif


      if (IS_OFF( def1(24) )) then
         lalo=.false.
      else if (IS_ON( def1(24) )) then
         lalo=.true.
         IPOS=IPOS+1
         IONAM(IPOS)='./latitude'
         IPOS=IPOS+1
         IONAM(IPOS)='./longitude'
      else
         write(6,6001) ' -lalo "no or yes" ?'
         call                                      xit('lire_arg',  -24)
      endif


      if(def1(26).eq.'ERR') then
* dans ce cas, la cle est presente mais il manque la valeur
         write(6,6001) ' -miss_ccc "VALEUR" ?'
         call                                      xit('lire_arg',  -26)
      else if (def1(26).ne.'?') then 
* cle et valeur sont presentes 
CCC     if (direction.eq.'cccma'   .or.
CCC  .      direction.eq.'rpncmc') then
           read(def1(26),9004,err=1000) miss_ccc
           miss_ccc_def=.true.
CCC     else
CCC      write(6,6001) ' -miss_ccc incompatible avec -dir'
CCC      call                                      xit('lire_arg',  -26)
CCC     endif
      endif

* Lecture de cle_nhem  (valeur de defaut=99)
* si direction=cccma necessaire pour grille polaire stereo seulement
* (car avant on assignait nhem=1 dans rdlatlon2.f).
         cle_nhem=99
         if(direction.eq.'cccma') then
            if(def1(28).ne. '?') then
               if (def1(28).eq.'1' .or. def1(28).eq.'2') then
                  read(def1(28),9003) cle_nhem
               else
                  write(6,6001) ' -cle_nhem avec valeur de 1 pour 
     .HEMIS NORD ou 2 pour HEMIS SUD si grille polaire stereographique
     . et direction=cccma'
                  call                             xit('lire_arg', -28)
               endif
            endif
         endif

***    print * , 'dans lire_arg.f CLE 28 def1(28)=======',def1(28) !debug
***    print * , 'dans lire_arg.f CLE 28 cle_nhem=======',cle_nhem !debug

      call get_environment_variable('UDUNITS2_XML_PATH',evalue,L_argenv)

      if (L_argenv > 0) then
         udunits_dat = evalue
      else
         if(def1(29).eq.'?') then
            write(6,6001) ' Fichier udunits2.xml ?'
            call                                   xit('lire_arg',  -29)
         else if(def1(29).eq.'default') then
            udunits_dat = udunits2_def
         else
            udunits_dat = def1(29)
         endif
      end if

      if(cdf2_mode.eq.'cdf2rpn') then
         rlonoff=-1000.0
         if (def1(30).ne.'?') read(def1(30),9004,err=1000) rlonoff
      endif

      if(def1(31)   == '?'             .and. 
     .   cdf2_mode  /= 'cdf2rpn'       .and.
     .   level_desc == 'Hybrid Levels') then
         write(6,6001) ' -hyb_pt "PT" (Pascal) ?'
         call                                      xit('lire_arg',  -31)
      else if(def1(31) /= '?') then
         read(def1(31),9004,err=1000) hyb_pt
      else
         hyb_pt = -1.0
      endif

      if(def1(32)   == '?'             .and. 
     .   cdf2_mode  /= 'cdf2rpn'       .and.
     .   level_desc == 'Hybrid Levels') then
         write(6,6001) ' -hyb_pref "PREF" (Pascal) ?'
         call                                      xit('lire_arg',  -32)
      else if(def1(32) /= '?') then
         read(def1(32),9004,err=1000) hyb_pref
      else
         hyb_pref = -1.0
      endif

      if(def1(33)   == '?'             .and. 
     .   cdf2_mode  /= 'cdf2rpn'       .and.
     .   level_desc == 'Hybrid Levels') then
         write(6,6001) ' -hyb_r "R" ?'
         call                                      xit('lire_arg',  -33)
      else if(def1(33) /= '?') then
         read(def1(33),9004,err=1000) hyb_r
      else
         hyb_r = -1.0
      endif

      if (def1(36).ne.'?' .and. direction.eq.'netcdf') then
         time_desc= def1(36)
         ok=.false.
         do i=1,ntun
            if(time_desc.eq.possible%time(i)) then
               ok=.true. ; exit
            endif
         enddo
         if(.not.ok) then
            write(6,6010) cles(36)//' = '//trim( def1(36) ),
     +                  ' Valeurs possibles de time_desc :'
            do i=1,ntun
               write(6,6009) trim(possible%time(i))
            enddo
            call                                   xit('lire_arg',  -36)
         endif
      else if (direction.eq.'netcdf') then
         write(6,6001) ' -timdesc "time_desc" ?'
         call                                      xit('lire_arg',  -36)
      else
         time_desc = ' '
      endif

      if (IS_ON( def1(37) )) then
         non_geographique = .true.
      else if (IS_OFF( def1(37) )) then
         non_geographique = .false.
      else
         write(6,6010) 'nongeog',trim( def1(37) )
         call                                      xit('lire_arg',  -37)
      endif

* lecture de la cle 'dtsize' seulement en mode 'cdf2rpn'

      if (cdf2_mode == 'cdf2rpn') then
         if (def1(42) /= '?') then
            read(def1(42),9004,err=300) dtsize
            goto 400
  300       dtsize = -1.0
  400       if (dtsize < 0) then
               write(6,6001)' -dtsize "delta accumul" '//trim( def1(42))
               call                                xit('lire_arg',  -42)
            else if (nint( dtsize*3600. ) < 1) then
               if (dtsize > 0.0_8)
     .            write(6,6007) 'Warning: dtsize reset to 0.0'
               dtsize = 0.0
            endif
         else
            write(6,6001)' -dtsize "delta accumul" ? '
            call                                   xit('lire_arg',  -42)
         endif
      endif

* lecture de la cle 'calendar', tout en tenant compte
* d'une possible definition prealable de la cle -leap

      evalue = def1(43) ; call up2low( evalue,evalue )

      if (evalue /= '?' .and. def1(4) == '?') then
         if      (evalue == 'standard'  .or.
     .            evalue == 'gregorian') then
            cccvx = .false. ; leap = .true. ; noUD = .false.
         else if (evalue == 'proleptic_gregorian') then
            cccvx = .false. ; leap = .true. ; noUD = .true.
         else if (evalue == 'noleap' .or.
     .            evalue == '365_day') then
            cccvx = .false. ; leap = .false.
         else if (evalue == '360_day') then
            cccvx = .true. ; leap = .false.
         else
            write(6,6010) 'calendar',trim( def1(43) )
            call                                   xit('lire_arg',  -43)
         endif
      else if (
     .  (def1(4) == 'on'  .and. (evalue == 'noleap'   .or.
     .                           evalue == '365_day'.  or.
     .                           evalue == '360_day'))     .or.
     .  (def1(4) == 'off' .and. (evalue == 'standard' .or.
     .                           evalue == 'gregorian'.or.
     .                           evalue == 'proleptic_gregorian'))
     .        ) then
         write(6,'(/A/)') 'conflit entre arguments -leap et -calendar'
         call                                      xit('lire_arg',  -43)
      else
         cccvx = .false. ; noUD = .true.
      endif
      
      if (leap) then
         call Accept_LeapYear()
      else
         call Ignore_LeapYear()
      end if

* (optionnellement) definir un code de table GRIB NCEP

      gribcode = -999
      if (def1(44) /= '?') read(def1(44),'(BN,I4)',err=1000) gribcode
      
* Toujours definir "cell_method" (voir CF-Metadata),
* meme si tout ce qu'il y a a transfere est un '?'

      cell_method = def1(45)

* (optionnellement) titre inserable dans la liste
* des attributs globaux du fichier NetCDF    

      meta_title = def1(46)

* (optionnellement) definir les noms des coordonnees NetCDF

      xcoord = def1(38)
      ycoord = def1(39)
      zcoord = def1(40)
      tcoord = def1(41)

      if(def1(27).eq.'ERR') then
* dans ce cas, la cle est presente mais il manque la valeur
         write(6,6001) ' -fill_ccc "VALEUR" ?'
         call                                      xit('lire_arg',  -27)
      else if (def1(27).ne.'?') then 
* cle et valeur sont presentes 
CCC     if (direction.eq.'cccma'   .or.
CCC  .      direction.eq.'rpncmc') then
           read(def1(27),9004,err=1000) fill_ccc
           fill_ccc_def=.true.
CCC     else
CCC      write(6,6001) ' -fill_ccc incompatible avec -dir'
CCC      Call                                      xit('lire_arg',  -27)
CCC     endif
      endif

* Assurer la consistance des options "fill" et "miss",
* tout en s'assurant que "fill" ait preseance lorsque
* defini comme argument. Aussi definir fill_toler

      if (fill_ccc_def .neqv. miss_ccc_def) then
         if (miss_ccc_def) then
            fill_ccc_def = .true. ; fill_ccc = miss_ccc
         endif
      endif

      if (fill_ccc_def) fill_toler = abs( fill_ccc ) * 0.001

      if(project%name.eq.'gaussian')then

            project%len =1

      elseif(project%name.eq.'lon/lat')then

         if(def1(18).eq.'?') then
            write(6,6001) ' -0lon "0LON" ?'
            call                                   xit('lire_arg',  -18)
         elseif(def1(18).eq.'GLOBAL') then
            project%len =1
            project%name ="lon/lat global B"
         else
            project%len=5
            project%name="lon/lat regional"
            project%nampar(1)="0lon"
            read(def1(18),9004,err=1000) project%value(1)
         endif

         if (project%name.eq."lon/lat regional")then
            if(def1(19).eq.'?') then
               write(6,6001) ' -0lat "0LAT" ?'
               call                                xit('lire_arg',  -19)
            elseif(project%name.ne."lon/lat global B") then
               project%nampar(2)="0lat"
               read(def1(19),9004,err=1000) project%value(2)
            endif

            if(def1(20).eq.'?') then
               write(6,6001) ' -dlon "DLON" ?'
               call                                xit('lire_arg',  -20)
            elseif(project%name.ne."lon/lat global B") then
               project%nampar(3)="dlon"
               read(def1(20),9004,err=1000) project%value(3)
            endif

            if(def1(21).eq.'?') then
               write(6,6001) ' -dlat "DLAT" ?'
               call                                xit('lire_arg',  -21)
            elseif(project%name.ne."lon/lat global B") then
               project%nampar(4)="dlat"
               read(def1(21),9004,err=1000) project%value(4)
            endif

         endif

      endif


      if (IS_OFF( def1(22) )) then
         invj=.false.
      else if (IS_ON( def1(22) )) then
         invj=.true.
      else
         write(6,6001) ' -invj "no or yes" ?'
         call                                      xit('lire_arg',  -22)
      endif

      if(.false.)then                                                    !debug
      write(6,*) attr_file                                               !debug
      write(6,*) "lire_arg : "                                           !debug
c      write(6,*) "grid_desc :", grid_desc                               !debug
      write(6,*) "porjection name:",project%name                         !debug
      write(6,*) "param. project% nbr:",project%len                      !debug
      do i=1,project%len                                                 !debug
         write(6,*)i,trim(project%nampar(i)),":",project%value(i)        !debug
      enddo                                                              !debug
      endif                                                              !debug

      DEALLOCATE ( ACLES,ADEF,ADEF1,ADEF2 )

      return

 1000 write(6,6099)
      call                                         xit('lire_arg',  -99)

*-----------------------------------------------------------------------
 6001 format(/" A l'appel, definir : ",a/)
 6007 format(33x,a)
 6009 format(34x,a)
 6010 format(/'Mauvaise valeur pour -',a/a)
 6099 format(/"Probleme de lecture d'arguments..."/)
 9003 format(i10)
 9004 format(e20.0)
*-----------------------------------------------------------------------
      end
      subroutine lprognam( nom )

      implicit none

      integer I,J,K,N
      character*(*) nom
      CHARACTER lnom*512

***    Recuperer le nom du programme.

      Call getarg( 0,lnom )

***    Enlever tout prefixe (positionel).

      N = len_trim( lnom )

      Do i=N,1,-1
         If (lnom(i:i).EQ.'/') Then
            Do j=i+1,N
               k         = j-i
               lnom(k:k) = lnom(j:j)
            End Do
            lnom(N-i+1:512) = ' '
            Exit
         End If
      End Do

***    Enlever tout postfixe debutant avec '_','-' ou '.'.

      J = index( lnom,'_' )
      If (j.ne.0) lnom = lnom(1:j-1)

      J = index( lnom,'-' )
      If (j.ne.0) lnom = lnom(1:j-1)

      J = index( lnom,'.' )
      If (j.ne.0) lnom = lnom(1:j-1)

***    Minimiser le nom du programme.

      CALL up2low( lnom, lnom )
      nom = lnom

      return
      end


