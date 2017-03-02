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

      subroutine encodate( DTIME, intime,timeunit )
*
      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'

      real(8) DTIME               ! The output NetCDF time value
      integer intime              ! The input CMC/RPN date-time-stamp
      character(len=128) timeunit ! The unit attribute for time in the netCDF file
 

******
*
*AUTEUR
*  Guy Bergeron  mai 2003
*
*     Definir un nombre d'heures ecoule depuis la date contenue dans 
*     timeunit.
*
*REVISIONS
*
*  Bernard Dugas, fevrier 2017 :
*   Tout le traitement de type udunits (y compris 
*   l'initialisation) se fait maintenant dans
*   la routine udparse3 qui remplace udparse2
*
*  Bernard Dugas, septembre 2014 :
*   Corriger le code d'initiation pour le cas .NOT.leap.
*
*  Bernard Dugas, juillet 2012 :
*   Tenir compte des differentes valeurs de delta%type dans
*   les calculs pour les calendriers '360_day', '365_day' et
*   'proleptic_gregorian'. 
*
*  Bernard Dugas, juin 2012 :
*   Faire la distinction entre les calendriers 'gregorian'/
*   'standard' et 'proleptic_gregorian': On ne peut utiliser
*   les fonctions basees sur UDUNITS avec 'proleptic_gregorian'
*
*  B. Dugas: avr 2012
*   Ajouter le calendrier 'proleptic_gregorian' (i.e. 'gregorian')
*
*  B. Dugas: dec 2011
*   Tenir compte des differents calendriers possibles, a savoir
*   'gregorian' (par defaut), '365_day' et '360_day'. Le dernier
*   est actif lorsque la variable CCCVX=T, et le second lorsque
*   les variables LEAP=F et CCCVX=F.
*
* B. Dugas: aout 2007
* - Renommer le fichier a encodate.ftn
* - L'inclusion de "udunits.inc" se fait via une pre-processeur.
*   Ceci permets d'utiliser le macro udunits UD_POINTER et de
*   de rendre le tout "machine-dependant".
* - Faire un "call xit" en cas d'erreur
*
******

      integer, parameter :: encod=-1, decod=+1

      logical, save :: first_call = .true.
      integer, save :: iyear0,imonth0,iday0,ihour0,imin0,isec0,intime0

      character(len=128) text
      integer iyear,imonth,iday,ihour,imin,isec
      integer status,datchek,part1,part2
      real    retcode
      real(8) hours

      integer, external :: newdate

*-----------------------------------------------------------------------
      if (first_call .and. (cccvx .or. .not.leap .or. noUD)) then

         hours = 0
         call udparse3(timeunit,udunits_dat,decod,hours,
     .                 iyear0,imonth0,iday0,ihour0,imin0,isec0)
         part1   =  (iyear0*100+imonth0 )*100+iday0
         part2   = ((ihour0*100+imin0)*100+isec0)*100
         datchek = newdate( intime0, part1,part2, +3 )
         if (datchek /= 0) call xit('encodate', -2 )

         if (cccvx) then ! Forcer un calendrier de 12 mois * 30 jours 
            if (imonth0 == 3) then ! que newdate supporte pas
               if (iday0 > 1) then
                  iday0 = iday0-1
               else
                  imonth0 = 2 ; iday0 = 30
               endif
            else if (imonth0 == 2) then
               iday0 = iday0+1
            else if (imonth0 == 1 .and. iday0 == 31) then
               imonth0 = 2 ; iday0 = 1
            endif
         else if (.not.leap) then
            call Ignore_LeapYear()
         endif

         first_call = .false.

      endif

      call def_date2( intime, iyear,imonth,iday, ! Decodons intime
     .                        ihour,imin,  isec, 'decode' )

      if (cccvx) then ! Calendrier 360 jours

         if (imonth == 3) then ! Forcer un calendrier de 12 mois * 30 jours 
            if (iday > 1) then
               iday = iday-1
            else
               imonth = 2 ; iday = 30
            endif
         else if (imonth == 2) then
            iday = iday+1
         else if (imonth == 1 .and. iday == 31) then
            imonth = 2 ; iday = 1
         endif

         dtime = (isec   - isec0  )/3600.0_8
     .         + (imin   - imin0  )/60.0_8
     .         + (ihour  - ihour0 )
     .         + (iday   - iday0  )*24.0_8
     .         + (imonth - imonth0)*720.0_8  ! (30*24)
     .         + (iyear  - iyear0 )*8640.0_8 ! (360*24)
         
         if (delta%type == 's') then
            dtime = dtime *  3600.
         else
     .   if (delta%type == 'm') then
            dtime = dtime *  60.
         else
     .   if (delta%type == 'd') then
            dtime = dtime /  24.
         else
     .   if (delta%type == 'M') then
            dtime = dtime / (24.*30.)
         else
     .   if (delta%type == 'Y') then
            dtime = dtime / (24.*30.*12.)
         endif

      else if (.not.leap .or. noUD) then

         ! Calendier 365 jours ou Gregorien Proleptique
         call difdatr( intime,intime0,dtime )

         if (delta%type == 's') then
            dtime = dtime *  3600.
         else
     .   if (delta%type == 'm') then
            dtime = dtime *  60.
         else
     .   if (delta%type == 'd') then
            dtime = dtime /  24.
         else
     .   if (delta%type /= 'h') then
            text="delta%type="//delta%type//
     .           " non valide avec calendrier courant"
            write(6,'(/A/)') trim( text )
            call xit('encodate',-2)
         endif

      else ! Calendrier Gregorien (par defaut)

         call udparse3( timeunit,udunits_dat,encod,DTIME,
     .                  iyear,imonth,iday,ihour,imin,isec )

      endif

*-----------------------------------------------------------------------
 2222 format(i10.10)
 2223 format(i4,3i2)
*-----------------------------------------------------------------------


      end
