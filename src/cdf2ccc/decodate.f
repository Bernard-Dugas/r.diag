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

      subroutine decodate( ncid, DTIME,ouTime )
*
      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'infomem.h'

      integer ncid
      real*8  DTIME
      
      integer ouTime

* AUTEUR Guy Bergeron   mai 2003
*
*     Definir la date dans un format date-time-stamp CMC/RPN (ouTime) a
*     partir de l'attribut "units" de la variable coordonnee "time" ssi
*     "units" est different de "unitless". Si non, ouTime est affectee
*     par DTIME.
*
*     Les valeurs de la variable de coordonnee "time" correspondent a des 
*     nombres de secondes, d'heures, ou de jours et ce en fonction de 
*     l'attribut "units" de la variable "time". Des exemples de l'attribut 
*     units de la variable "time" :
*     
*                   "seconds since 1800-1-1 00:00:0.0"   
*                   "hours since 1800-1-1 00:00:0.0"   
*                   "days since 1800-1-1 00:00:0.0"   
*                   "months since 1800-1-1 00:00:0.0"   
*     
*     NOTA: unitstring contient la date de reference pour le nombre (de 
*           seconde, d'heure, de jour ou de mois) ecoules.
*
*REVISION:
*
*     Aout 2018 - Bernard Dugas :
*           Passer directement coord(tid)%name a la fonction
*           nf_inq_varid plutot qu'une variable en minuscule
*     Fevrier 2017 - Bernard Dugas :
*           tout le traitement de type udunits (y compris 
*           l'initialisation) se fait maintenant dans
*           la routine udparse3 qui remplace udparse2
*     Juillet 2012 - Bernard Dugas :
*           le traitement de l'option 'unitstring = "day as %y%m%d.%f"'
*           est corrige et rendu independant du calendrier specifie
*
*     Juin  2012 - Bernard Dugas :
*         - faire la distinction entre les calendriers
*           'gregorian'/'standard' et 'proleptic_gregorian',
*           i.e. on ne peut utiliser les fonctions basees sur
*           UDUNITS avec un calendrier 'proleptic_gregorian'
*         - supporter unitstring = "day as %y%m%d.%f" (pour ECMWF ?)
*
*     Avr.  2012 - Bernard Dugas :
*           ajouter le calendrier 'proleptic_gregorian' (i.e. 'gregorian')
*
*     Dec.  2011 - Bernard Dugas :
*           correction des 29 et 30 fevrier au mode '360_day'
*
*     Jul.  2011 - Bernard Dugas :
*           tenir compte de l'attribut calendar, qui peut
*           servir a identifier la longueur d'une annee...
*           supporte les valeurs 'gregorian' (par defaut),
*           '365_day' et '360_day'
*
*     Mar.  2010 - Bernard Dugas :
*           nf_get_att_text utilise coord(tid)%name plutot
*           que 'time' lorsqu'on cherche les unites de la
*           variable temporelle
*
*     Nov.  2009 - Bernard Dugas :
*           ouTime est maintenant calcule par def_date2
*           puisqu'on le garde format DateTimeStamp
*
*     Fev.  2009 - Bernard Dugas :
*           ajouter la routine DECODATE2 (identique a 
*           DECODATE sauf retour ouTime en INTEGER*8)
*
*     Mai   2008 - Bernard Dugas :
*           appel a UDPARSE2 generalise
*           ajout de datec pour le cas "julian day"
*
*     Aout  2007 - B. Dugas :
*         - Renommer le fichier a decodate.ftn
*         - L'inclusion de "udunits.inc" se fait via un pre-processeur.
*           Ceci permets d'utiliser le macro udunits UD_POINTER et de
*           de rendre le tout "machine-dependant".
*         - Supporter le type "julian day"
*
*     Mai   2005 - Guy Bergeron : 
*           ajout de cal_date pour le cas "months since AAAA-MM-DD"
*           ajout de cal_date pour le cas "Season DJF MAM JJA SON"
*
*     Juin  2004 - Guy Bergeron : 
*           agrument d'appel en real*8 DTIME au lieu de real*4 xtime
*
*     Juil. 2004 - Guy Bergeron : 
*           appel de udparse conditionnel a la valeur de "units"
* 
******

      integer, parameter :: sense=+1

      integer timeid,iyear,imonth,iday,ihour,iminute,isecond
      integer part1,part2,timestamp2,ier,nlen,status
      integer, save :: timestamp1=0
      integer iii(32),i,j
      real retcode

      integer, external :: newdate
      logical, save :: CheckLeapYears=.true.
      logical  NoLeapYear

      real(8)       hours,fraction
      integer(8)    elapsed,nombre
      character*128 unitstring,calendarstring,string

      data  unitstring /''/      ! initialisation a vide 
*-----------------------------------------------------------------------
      calendarstring=' '

*     Extraire l'attribut "units" de la coordonnee "time" :

      status=nf_inq_varid(ncid, coord(tid)%name, timeid)
      if (status == nf_noerr) then
         status=nf_get_att_text(ncid,timeid,'units',unitstring)
         call handle_err2(status,'decodate')
         status=nf_get_att_text(ncid,timeid,'calendar',calendarstring)
         if (status /= nf_noerr) calendarstring = 'gregorian'
      else
         unitstring = 'unitless'
         calendarstring = 'gregorian'
      endif

      call clean_char( calendarstring,string,nlen )
      calendarstring = string(1:nlen)

*     decode la date (si necessaire) :

      call up2low( calendarstring,calendarstring )
      call up2low( unitstring,unitstring )

      if(unitstring.eq."unitless") then
         ouTime=DTIME

      else if (unitstring(1:06)=="day as") then

         if (unitstring(8:16) == "%y%m%d.%f") then
            if (DTIME==0) then
               ouTime = -1 
               return
            else
               nombre = DTIME ; fraction = DTIME-nombre
               write(string,'(I8.8)') nombre
               read(string,'(i4.4,i2.2,i2.2)') iyear,imonth,iday
               fraction=int( 86400*fraction )
               ihour   = fraction/3600
               iminute =(fraction-3600*ihour)/60
               isecond = fraction-3600*ihour -60*iminute
            endif
         endif

         call def_date2( ouTime, iyear,imonth, iday,
     .                           ihour,iminute,isecond, 'encode' )

      else

         ihour=0 ; iminute = 0 ; isecond = 0
         
         if (calendarstring == 'standard'   .or.
     .       calendarstring == 'gregorian') then

            if (CheckLeapYears) then
               call Get_LeapYear_Status( NoLeapYear )
               if (NoLeapYear) then
                  write(6,6101)
                  call Accept_LeapYear()
               endif
               CheckLeapYears = .false.
            endif

            if (unitstring(1:6).eq."season") then
               call cal_date(unitstring,DTIME,iyear,imonth,iday,ihour)

            else if (unitstring(1:10).eq."julian day") then
               call datec( nint( DTIME ), iyear,imonth,iday )

            else if (unitstring(1:12).eq."months since" .or.
     .               unitstring(1:11).eq."years since") then
               call cal_date(unitstring,DTIME,iyear,imonth,iday,ihour)

            else if (unitstring(1:10).eq."days since"     .or.
     .               unitstring(1:11).eq."hours since"    .or.
     .               unitstring(1:13).eq."minutes since"  .or.
     .               unitstring(1:13).eq."seconds since") then
               call udparse3(unitstring,udunits_dat,sense,DTIME,
     .                       iyear,imonth,iday,ihour,iminute,isecond)

            endif
               
         else if (calendarstring == 'noleap'  .or.
     .            calendarstring == '365_day' .or.
     .            calendarstring == 'proleptic_gregorian') then

*           Une annee ayant toujours 365 jours (pas de jours bissextiles)
*           ou le calendrier Greogorien Proleptique (avec jours bissextiles)

            if (CheckLeapYears) then
               call Get_LeapYear_Status( NoLeapYear )
               if (calendarstring == 'noleap'   .or.
     .             calendarstring == '365_day') then
                  if (.not.NoLeapYear) then
                     write(6,6100)
                     call Ignore_LeapYear()
                  endif
               else if (calendarstring == 'proleptic_gregorian') then
                  if (NoLeapYear) then
                     write(6,6101)
                     call Accept_LeapYear()
                  endif
               endif
               CheckLeapYears = .false.
            endif

            if (unitstring(1:10) /= 'days since'     .and.
     .          unitstring(1:11) /= 'hours since'    .and.
     .          unitstring(1:13) /= 'minutes since'  .and.
     .          unitstring(1:13) /= 'seconds since') then
               write(6,6001) trim( unitstring )
               call                                xit('decodate', -1 )
            endif

            if (timestamp1 == 0) then

*              decoder l'origine de la coordonnee temporelle

               hours = 0
               call udparse3(unitstring,udunits_dat,sense,hours,
     .                       iyear,imonth,iday,ihour,iminute,isecond)
            
               part1 =  (iyear*100+imonth )*100+iday
               part2 = ((ihour*100+iminute)*100+isecond)*100
               ier   = newdate( timestamp1, part1,part2, +3 )

            endif

            if      (unitstring(1:10) ==    'days since') then
               hours = DTIME*24.
            else if (unitstring(1:11) ==   'hours since') then
               hours = DTIME
            else if (unitstring(1:13) == 'minutes since') then
               hours = DTIME/60.
            else if (unitstring(1:13) == 'seconds since') then
               hours = DTIME/3600.
            endif

            call incdatr( timestamp2,timestamp1, hours )
            ier   = newdate( timestamp2, part1,part2, -3 )

            iyear   =      part1/10000
            imonth  = mod( part1/100,100 )
            iday    = mod( part1    ,100 )
            ihour   =      part2/1000000
            iminute = mod( part2/10000,100 )
            isecond = mod( part2/100,100 )

         else if (calendarstring == '360_day') then

*           Au depart, on a une annee de 12 mois ayant chacun 30 jours.

*           On modifie ca pour se conformer aux routines de dates
*           calendrier de la facon suivante: janvier et mars auront
*           31 jours et fevrier aura toujours 28 jours. Tous les
*           autres mois auront 30 jours.

            if (CheckLeapYears) then
               call Get_LeapYear_Status( NoLeapYear )
               if (.not.NoLeapYear) then
                  call Ignore_LeapYear()
               endif
               CheckLeapYears = .false.
            endif

            if (unitstring(1:11) /=   'years since'  .and. 
     .          unitstring(1:12) /=  'months since'  .and.
     .          unitstring(1:10) /=    'days since'  .and.
     .          unitstring(1:11) /=   'hours since'  .and.
     .          unitstring(1:13) /= 'minutes since'  .and.
     .          unitstring(1:13) /= 'seconds since') then
               write(6,6001) trim( unitstring )
               call                                xit('decodate', -1 )
            endif

*           decoder l'origine de la coordonnee temporelle

            hours = 0
            call udparse3(unitstring,udunits_dat,sense,hours,
     .                    iyear,imonth,iday,ihour,iminute,isecond)
            
            if      (unitstring(1:11) ==   'years since') then
               elapsed  = nint( DTIME*360*86400,8 )
            else if (unitstring(1:12) ==  'months since') then
               elapsed  = nint( DTIME*30*86400,8 )
            else if (unitstring(1:10) ==    'days since') then
               elapsed  = nint( DTIME*86400,8 )
            else if (unitstring(1:11) ==   'hours since') then
               elapsed  = nint( DTIME*3600,8 )
            else if (unitstring(1:13) == 'minutes since') then
               elapsed  = nint( DTIME*60,8 )
            else if (unitstring(1:13) == 'seconds since') then
               elapsed  = nint( DTIME,8 )
            endif

            isecond     = isecond + mod( elapsed, 60_8 )
            elapsed     = elapsed/60

            if (isecond > 59) then
               isecond  = isecond-60
               iminute  = iminute+1
            endif

            iminute     = iminute + mod( elapsed, 60_8 )
            elapsed     = elapsed/60

            if (iminute > 59) then
               iminute  = iminute-60
               ihour    = ihour+1
            endif

            ihour       = ihour +  mod( elapsed, 24_8 )
            elapsed     = elapsed/24

            if (ihour   > 23) then
               ihour    = ihour-24
               iday     = iday+1
            endif

            iday        = iday + mod( elapsed, 30_8 )
            elapsed     = elapsed/30

            if (iday    > 30) then
               iday     = iday-30
               imonth   = imonth+1
            endif

            imonth      = imonth + mod( elapsed,12_8 )
            elapsed     = elapsed/12

            if (imonth  > 12) then
               imonth   = imonth-12
               iyear    = iyear+1
            endif

            iyear       = iyear + elapsed

*           Correction pour les 29 et 30 fevrier:
*           Deux journees sont enlevees de fevrier et
*           on ajoute les 31 janvier et 31 mars.

            if (imonth == 2) then
               if (iday == 1) then
                  imonth = 1 ; iday = 31
               else if (iday < 30) then
                  iday = iday-1
               else
                  imonth = 3 ; iday = 1
               endif
            else if (imonth == 3) then
               iday = iday+1
            endif

         endif

         call def_date2( ouTime, iyear,imonth, iday,
     .                           ihour,iminute,isecond, 'encode' )

      endif

      if (.false.)write(6,*)' decodate :',trim(unitstring),': ',      !debug
     .                        DTIME,' : ',ouTime                      !debug
*-----------------------------------------------------------------------
 6001 format(/'DECODAT: Un-supported time units with',
     .        ' non-gregorian calendar... ',A/)
 6100 format(/' In DECODATE: Leap Years will be ignored...'//)
 6101 format(/' In DECODATE: Leap Years will be accounted for...'//)

      end
      subroutine decodate2( ncid, DTIME,ouTime )

      implicit none

*     Definir la date dans le format AAAAMMDDHH (outime) a partir de 
*     l'attribut "units" de la variable coordonnee "time" ssi "units" est 
*     different de "unitless". Si non, outime est affectee par dtime.

      integer, parameter :: sense=1

      integer   ncid
      real*8    DTIME
      integer*8 ouTime
 
      include 'netcdf.inc'
      include 'cdf2ccc.h'

      integer timeid,iyear,imonth,iday,ihour,iminute,isecond,status
      integer*8, parameter :: cent = 100
      real retcode
      
      character*128 unitstring,string

      data  unitstring /''/      ! initialisation a vide 
*-----------------------------------------------------------------------

*     Extraire l'attribut "units" de la coordonnee "time" :

      status=nf_inq_varid(ncid,'time', timeid)
      if (status == nf_noerr) then
         status=nf_get_att_text(ncid,timeid,'units',unitstring)
         call handle_err2(status,'decodate2')
      else
         unitstring = 'unitless'
      endif

*     decode la date (si necessaire) :

      call up2low( unitstring,unitstring )

      if(unitstring.eq."unitless") then
         ouTime=DTIME

      else

         ihour=0 ; iminute = 0 ; isecond = 0
         
         if (unitstring(1:6).eq."season") then
            call cal_date(unitstring,DTIME,iyear,imonth,iday,ihour)

         else if (unitstring(1:10).eq."julian day") then
            call datec( nint( DTIME ), iyear,imonth,iday )

         else if (unitstring(1:12).eq."months since" .or.
     .            unitstring(1:11).eq."years since") then
            call cal_date(unitstring,DTIME,iyear,imonth,iday,ihour)

         else if (unitstring(1:10).eq."days since"     .or.
     .            unitstring(1:11).eq."hours since"    .or.
     .            unitstring(1:13).eq."minutes since"  .or.
     .            unitstring(1:13).eq."seconds since") then
            call udparse3(unitstring,udunits_dat,sense,DTIME,
     .                    iyear,imonth,iday,ihour,iminute,isecond)

         else
            call udparse3(unitstring,udunits_dat,sense,DTIME,
     .                   iyear,imonth,iday,ihour,iminute,isecond)

         endif

         ouTime =          isecond 
     .          + cent * ( iminute
     .          + cent * ( ihour
     .          + cent * ( iday 
     .          + cent * ( imonth 
     .          + cent *   iyear ) ) ) )

      endif

      if (.false.)write(6,*)" decodate2 :",trim(unitstring),": ",     !debug
     .                        DTIME," : ",ouTime                      !debug
*-----------------------------------------------------------------------
      end
