!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     This subroutine is part of the netopen Fortran callable CDC
!     group of routines for reading generic netCDF files in the new
!     cooperative standard as presented in:
!	http://ferret.wrc.noaa.gov/noaa_coop/coop_cdf_profile.html
!
!     This routine takes a time value pulled from the netCDF file,
!     determines if it is in the new or old time format (by using the
!     udunits function) and then returns the representative year,
!     month, day and hour represented by the time value.
!
!     Written by Cathy Smith of CDC on ???
!     Modified by Tom Baltzer of CDC on Feb. 8, 1995
!	- To determine if time is old or new format and parse
!		accordingly
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine udparse3( TimeUnit,udunits_dat,sense,ioTime, &
                           oyear,omonth,oday,ohour,ominute,osecond )

      use f_udunits_2
      implicit none

!*****
!
!AUTEUR
! G. Bergeron  juin 2005
!
!    Modification de la valeur de intime si months. Patch pour regler
!    un bug ou il manquait une seconde apres conversion de intime en 
!    date. J'ai ajoute la plus petite valeur pour que sa fonctionne.
!
!REVISIONS 
!
!   Bernard Dugas  fev  2017 :
!   - Ajouter les arguments "udunits_dat" et "sense". Ce
!     dernier determine dans quel sens se fait le travail
!   - Convertir a l'interface f_udunits_2 de MFV qui
!     gere tous les appels FORTRAN a udunits v2
!   - Coder avec FORTRAN 90/95
!   - Renommer a UDPARSE3
!   Bernard Dugas  sept 2014 :
!   - Corriger la declaration de TimeUnit
!   Bernard Dugas  juin 2012 :
!   - Corriger les arguments a XIT
!   Bernard Dugas  avril 2011 :
!   - Corriger l'appel a utcaltime
!   Bernard Dugas  mai  2008 :
!   - Renommer udparse2
!   - Ajouter les arguments ominute,osecond
!   Bernard Dugas  aout 2007 :
!   - Renommer le fichier a udparse.ftn
!   - L'inclusion de "udunits.inc" se fait via un pre-processeur.
!     Ceci permets d'utiliser le macro udunits UD_POINTER et de
!     de rendre le tout "machine-dependant".
!   - Faire un "call xit" en cas d'erreur
!
!*****

      ! Arguments I/O

      ! The unit attribute for time in the netCDF file
      character(len=*), intent(IN)    :: TimeUnit
      ! The UDUNITS2_XML_PATH pointing to the udunits2.xml file
      character(len=*), intent(IN)    :: udunits_dat
      ! The input/output time value (depending on the value of sense)
      real(8),            intent(INOUT) :: ioTime
      ! The input/output year,month,day and hour
      integer,            intent(INOUT) :: oyear,omonth,oday,ohour,ominute,osecond
      ! sense=+1 (decode ioTime), sense=-1 (encode ioTime)
      integer,            intent(IN)    :: sense

      ! Variables locales

      type(UT_SYSTEM_PTR),    save :: unitSystem
      type(CV_CONVERTER_PTR), save :: converter,iconverter
      type(UT_UNIT_PTR),      save :: secondu,unit
      character(len=512),     save :: TimeUnit0, udunits_dat0
      logical,                save :: boot=.true.,oldCDC=.false.

      integer                      :: charset,status
      real(8)                      :: xx,second,resolution,udTime

!-----------------------------------------------------------------------
      if (boot) then

         ! Initialisation pour UDUNITS v2

         unitSystem   = UT_SYSTEM_PTR_NULL
         secondu      = UT_UNIT_PTR_NULL
         unit         = UT_UNIT_PTR_NULL
         secondu      = UT_UNIT_PTR_NULL

         charset      = UT_ASCII
         TimeUnit0    = TimeUnit
         udunits_dat0 = udunits_dat
         
         unitSystem = f_ut_read_xml( trim( udunits_dat0 ) )

         if (UD_is_null( unitSystem ) ) then
            call ud_parse_error( status )
            if (status /= UT_SUCCESS) call xit('udparse3',-1 )
         end if

         secondu = f_ut_get_unit_by_name( unitSystem, "second" )
         if (UD_is_null( secondu ) ) then
            call ud_parse_error( status )
            if (status /= UT_SUCCESS) call xit('udparse3',-2 )
         end if

         unit = f_ut_parse( unitsystem,trim( TimeUnit0 ),charset )
         if (UD_is_null( unit ) ) then
            call ud_parse_error( status )

            if (sense  == +1         .and. &
                status == UT_UNKNOWN .and. &
               (TimeUnit(1:1) == 'y' .or. TimeUnit(1:1) == 'Y')) then

               ! Assume old CDC standard
               ! if time unit is unknown
               oldCDC = .true.

            else if (status /= UT_SUCCESS) then
               call xit('udparse3',-3 )
            end if
         end if

         if (.not.oldCDC) then

            converter = f_ut_get_converter( unit,secondu )
            if (UD_is_null( converter ) ) then
               call ud_parse_error( status )
               if (status /= UT_SUCCESS) call xit('udparse3',-4 )
            end if

            iconverter = f_ut_get_converter( secondu,unit )
            if (UD_is_null( iconverter ) ) then
               call ud_parse_error( status )
               if (status /= UT_SUCCESS) call xit('udparse3',-5 )
            end if

         end if

         boot = .false. ! Pour ne plus repasser ici

      else if (TimeUnit /= TimeUnit0) then

         print *,'More than one NetCDF time unit detected...'
         print *,'initial... ',trim( TimeUnit0 )
         print *,'current... ',trim( TimeUnit )
         call xit('udparse3',-6 )

      end if

      if (oldCDC) then

         ominute=0 ; osecond=0 ; xx = 0.

         oyear  = nint(  ioTime    /10000000000._8 ) ; xx =    oyear*10000000000._8
         omonth = nint( (ioTime-xx)/100000000._8 )   ; xx = xx+omonth*100000000._8
         oday   = nint( (ioTime-xx)/1000000._8 )     ; xx = xx+oday*1000000._8
         ohour  = nint( (ioTime-xx)/10000._8 )

      else if (sense == +1) then

         ! On decode ioTime dans oyear,omonth,oday,ohour,ominute,osecond

         udTime = f_cv_convert_double( converter, ioTime )
         call ud_parse_error( status )
         if (status /= UT_SUCCESS) call xit('udparse3',-7 )

         call f_ut_decode_time( udTime, &
              oyear,omonth,oday,ohour,ominute,second,resolution )
         call ud_parse_error( status )
         if (status /= UT_SUCCESS) call xit('udparse3',-8 )

         osecond = nint( second )

      else if (sense == -1) then

         ! On encode oyear,omonth,oday,ohour,ominute,osecond dans ioTime

         second = osecond
         udTime = f_ut_encode_time( oyear,omonth,oday,ohour,ominute,second )
         call ud_parse_error( status )
         if (status /= UT_SUCCESS) call xit('udparse3',-9 )

         ioTime = f_cv_convert_double( iconverter, udTime )
         call ud_parse_error( status )
         if (status /= UT_SUCCESS) call xit('udparse3',-10)

      else
         print *,'Unknown value for argument sense =',sense
         call xit('udparse3',-11)
      end if

      if (.false.)write(6,6001)oyear,omonth,oday,ohour, &        !debug
                               ominute,osecond                   !debug
 6001 format("UDPARSE3 :", &                                     !debug
             i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",f5.2)  !debug
      
      return
!-----------------------------------------------------------------------

   end subroutine udparse3
   subroutine ud_parse_error( status )

      use f_udunits_2

      integer status

!-----------------------------------------------------------------------
      status = f_ut_get_status()
 
      if (status == UT_SUCCESS)         return
      if (status == UT_BAD_ARG)         print *,"UT_BAD_ARG: An argument violates the function's contract"
      if (status == UT_EXISTS)          print *,'UT_EXISTS: Unit, prefix, or identifier already exists'
      if (status == UT_NO_UNIT)         print *,'UT_NO_UNIT: No such unit exists'
      if (status == UT_OS)              print *,'UT_OS: Operating-system error.  See "errno".'
      if (status == UT_NOT_SAME_SYSTEM) print *,'UT_NOT_SAME_SYSTEM: The units belong to different unit-systems'
      if (status == UT_MEANINGLESS)     print *,'UT_MEANINGLESS: The operation on the unit(s) is meaningless'
      if (status == UT_NO_SECOND)       print *,"UT_NO_SECOND: The unit-system doesn't have a unit named "//'"second"'
      if (status == UT_VISIT_ERROR)     print *,'UT_VISIT_ERROR: An error occurred while visiting a unit'
      if (status == UT_CANT_FORMAT)     print *,"UT_CANT_FORMAT: A unit can't be formatted in the desired manner"
      if (status == UT_SYNTAX)          print *,'UT_SYNTAX: String unit representation contains syntax error'
      if (status == UT_UNKNOWN)         print *,'UT_UNKNOWN: String unit representation contains unknown word'
      if (status == UT_OPEN_ARG)        print *,"UT_OPEN_ARG: Can't open argument-specified unit database"
      if (status == UT_OPEN_ENV)        print *,"UT_OPEN_ENV: Can't open environment-specified unit database"
      if (status == UT_OPEN_DEFAULT)    print *,"UT_OPEN_DEFAULT: Can't open installed, default, unit database"
      if (status == UT_PARSE)           print *,'UT_PARSE: Error parsing unit specification'

      return
!-----------------------------------------------------------------------

   end subroutine ud_parse_error



