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
      subroutine scanfile(FUNIT,IBUF,ntime,nlev,lablvl)

      implicit none

      include 'cdf2ccc.h'
      include 'workmem.h'
      include 'varmem.h'

      integer FUNIT,IBUF(*)
      integer ntime,nlev,lablvl(maxlev)
      
******
*
*AUTEUR Guy Bergeron   juin 2003
*
*     Evaluer le nombre d'enregistrement temporelle et leurs valeurs. 
*     Evaluer le nombre de nivaux et leurs valeurs. Definir les extrema 
*     des valeurs de la variable contenue dans funit :
*
*REVISIONS
*
*  Bernard Dugas juin 2013 :
*  - Tenir compte des valeurs manquantes lors du calcul de infvar(:)%range
*  Bernard Dugas aout 2012 :
*  - Tenir compte du jour de depart dans les calculs des time_bnds.
*    Ceci corrige un probleme avec les echantillons quotidiens
*  - Lorsque nsample > 1 et que la valeur de nhour=0 au premier echantillon
*    (i=1), on suppose que le dateo du moment corresponds a time_bnds(1,1)
*  Bernard Dugas juin 2012 :
*  - Ameliorer les calculs des time_bnds
*  - Utiliser les champs infvar(:)%niv, infvar(:)%pas,
*    infvar(:)%ip3 et infvar(:)%dateo pour
*    definir rtime, time_bnds et lablvl
*  - Utiliser GETSET2 plutot que RECGET/RECUP2
*  Bernard Dugas mai 2012 :
*  - Calculer les limites temporelles associees a des
*    moyennes de plusieurs echantillons (time_bnds)
*  Bernard Dugas novembre 2009 :
*  - On ne calcule plus stamp lorsque tlbl est vrai
*    puisque ibuf(2) est alors deja dans le bon format
*  Bernard Dugas fevrier 2009 :
*  - Tenir compte des minutes et secondes dans le
*    calcul de rtime pour les fichiers CMC/RPN
*  Bernard Dugas octobre 2008 :
*  - Limiter les messages notant la presence de plus d'un delta-T
*  Bernard Dugas juillet 2008 :
*  - Tenir compte de delta%type dans le calcul de rtime (si tlbl est faux)
*  Bernard Dugas mai 2008 :
*  - Supporter delta%type={s,m,h,d,M,y}
*  Bernard Dugas automne/hiver 2007 :
*  - Enlever la mod de juillet 2004, le champs de montagne
*    etant maintenant lu dans un autre fichier
*  - Les mix/max sont sauves dans infvar(:)%range(1:2)
*  Guy Bergeron juillet 2004 : Hybrid height coordinate
*
******

      real*8 xmin   ! valeur minimale
      real*8 xmax   ! valeur maximale

      integer i,cvar,nwds,itime

      character*4 type,name
      integer  iname,iiname,initl, nsample,ier,nk,nt,iset
      integer  year, month, day, hour, minute, sec,idelta
      integer nyear,nmonth,nday,nhour,nminute,nsec

      integer dateo,dateom, part1,part2,   reste,count, istart
      integer datev,datevp, lev(max_levs), nlev0,ndays
      real*8  ddelta,date8,nhour8,pdelta
      integer(8) :: dixala8=100000000_8

      character*128 string,dummy
      character(16) date16

      logical ok,deja_annonce,zero_initial

      integer, external :: newdate,gethigh
*-----------------------------------------------------------------------
      deja_annonce = .false. ; zero_initial = .false.

      ntime=0
      nlev=0

      delta%dval=0

                            dummy='hours since '              
      if (time_desc.ne.' ') dummy=trim( time_desc )//' since '              

      if (dummy.eq.'seconds since')  delta%type='s'
      if (dummy.eq.'minutes since')  delta%type='m'
      if (dummy.eq.'hours since'  )  delta%type='h'
      if (dummy.eq.'days since'   )  delta%type='d'
      if (dummy.eq.'months since' )  delta%type='M'         
      if (dummy.eq.'years since'  )  delta%type='y'         
         
      call def_date2( ladate, year,month,day,hour,minute,sec, 'decode' )

      if (sec == 0) then
         write(string,7779) trim( dummy ),year,month,day,hour,minute
      else
         write(string,7780) trim( dummy ),year,month,day,hour,minute,sec
      endif

      do cvar=1,maxvar

         if (.not.infvar(cvar)%var_ok) cycle

         nt    = infvar(cvar)%tid
         ntime = max( ntime, infvar(cvar)%len(nt) )

         nk    = infvar(cvar)%zid
         nlev  = max( nlev, infvar(cvar)%len(nk) )

      enddo

      ! Definir rtime, ddelta et (optionellement) time_bnds
 
      infvar(1:maxvar)%time_bnds_L = .false.

      do cvar=maxvar,1,-1

         if (.not.infvar(cvar)%var_ok) cycle

         nt = infvar(cvar)%tid
         
         if (infvar(cvar)%len(nt) /= ntime) cycle

         datev = -1

         do i=1,ntime

            if (tlbl) then        

*              Convertir infvar(cvar)%pas(i) en terme
*              de nbre d'heures ecoules depuis "LADATE".

               write(date16,'(I16.16)') infvar(cvar)%pas(i)
               read(date16,'(I4.4,2I2.2,3I2.2)') nyear,nmonth, nday,
     .                                           nhour,nminute,nsec

               if (delta%type == 'y') then

                  if (month  == nmonth  .and.
     .                day    == nday    .and.
     .                hour   == nhour   .and.
     .                minute == nminute .and.
     .                sec    == nsec   ) then
                     rtime(i) = nyear-year
                  else
                     write(6,6001)
                     call                          xit('scanfile',-1)
                  endif

               else if (delta%type == 'M') then

                  if (day    == nday    .and.
     .                hour   == nhour   .and.
     .                minute == nminute .and.
     .                sec    == nsec   ) then
                     rtime(i) = (nyear-year)*12 + (nmonth-month)
                  else
                     write(6,6001)
                     call                          xit('scanfile',-1)
                  endif

               else

                  part1 = infvar(cvar)%pas(i) / dixala8
                  part2 = infvar(cvar)%pas(i) - part1 * dixala8

                  datevp = datev ! Sauver la valeur precedente de datev

                  ier = newdate( datev, part1,part2, +3 )
                  if (ier /= 0) call               xit('scanfile',-2)

                  call encodate(rtime(i),datev,string) ! selon UDUNITS

                  if (i > 1 .and. .not.infvar(cvar)%time_bnds_L) cycle

                  if (ccc_pktyp(1:2) == 'SQ' .or. dtsize > 0.00001) then

                     nsample = infvar(cvar)%ip3(i)

                     if (dtsize > 0.00001 .or. nsample > 1) then

                        ! Calculer les limites temporelles associees
                        ! a des moyennes de plusieurs echantillons.

                        if (i == 1) infvar(cvar)%time_bnds_L = .true.

                        if (dtsize > 0.00001) then

                           ddelta = -dtsize
                           call incdatr( dateo,datev, ddelta )

                           ier   = newdate( dateo,part1,part2, -3 )
                           if (ier /= 0) call      xit('scanfile',-4)

                           call encodate( date8,dateo,string )

                        else

                           dateo = infvar(cvar)%dateo(i)
                           ier   = newdate( dateo,part1,part2, -3 )
                           if (ier /= 0) call      xit('scanfile',-5)

                           ! Calculons le nombre d'heures entre le moment
                           ! courant et le jour 01 et l'heure 00 du mois. On
                           ! suppose que tous les ensembles d'echantillons 
                           ! commencent de telle sorte

                           ndays = mod( part1 , 1 00 )
                           nhour = (ndays - 1) * 24 + part2 / 1 00 00 00

                           call difdatr( datev,dateo, ddelta )

                           ! Definir l'interval de sauvetage

                           nhour8 = ddelta/(nsample-1)

                           if (nhour > 0 .and. i == 1) then

                              count = nhour

                              do while (count >= 24 .and. ndays > 1)
                                 ndays = ndays - 1
                                 count = count - 24
                              enddo

                              if (count >= 24)
     .                        infvar(cvar)%time_bnds_L = .false.

                              part1 =(part1 / 1 00) * 1 00 +  ndays
                              reste = part2 - count * 1 00 00 00

                              ! Enlever nhour a dateo si nhour est egal
                              ! a la taille de l'interval de sauvetage
                              ! pour definir time_bnds(1,i)

                              if (nhour == nint( nhour8 )) then

                                 ier = newdate( dateom,part1,reste, +3 )
                                 if (ier /= 0) call xit('scanfile',-6)
                                 call encodate( date8,dateom,string )

                              else
                                 infvar(cvar)%time_bnds_L = .false.
                              endif

                           else if (i == 1) then

                              ! Denoter qu'on ne devra pas tenir compte
                              ! d'un offset de nhours8 par la suite

                              zero_initial = .true.

                              ! Pas besoin de verifier pourquoi nhour /= 0
                              ! et on suppose que dateo corresponds deja au
                              ! debut du premier interval ("time_bnds")

                              call encodate( date8,dateo,string )

                           else

                              if (zero_initial) then ! Pas de calcul d'offset
                                 dateom = dateo
                              else
                                 ! si dateom (= dateo - l'interval de sauvetage)
                                 ! est egal au datev precedent, dateom definira
                                 ! le time_bnds(1,i)
                                 call incdatr( dateom, dateo,-nhour8 )
                              endif

                              if (dateom == datevp) then
                                 call encodate( date8,dateom,string )
                              else
                                 infvar(cvar)%time_bnds_L = .false.
                              endif

                           endif

                        endif

                        if (infvar(cvar)%time_bnds_L) then
                           time_bnds(1,i) = date8
                           time_bnds(2,i) = rtime(i)
                        endif

                     else
                        infvar(cvar)%time_bnds_L = .false.
                     endif

                  endif

               endif

            else

*        Si etiquette temporel du fichier cccma = nbre de dt

               if      (delta%type == 's') then
                  rtime(i)=dble(infvar(cvar)%pas(i))*dt
               else if (delta%type == 'm') then
                  rtime(i)=dble(infvar(cvar)%pas(i))*dt/60.
               else if (delta%type == 'h') then
                  rtime(i)=dble(infvar(cvar)%pas(i))*dt/3600.
               else if (delta%type == 'd') then
                  rtime(i)=dble(infvar(cvar)%pas(i))*dt/86400.
               else
               
                  write(6,6002)  trim( dummy )
                  call                             xit('scanfile',-2)

               endif

            endif

            if (i > 1) then
               delta%dval=rtime(i)-rtime(i-1)
               if (i > 2                                          .and. 
     .            .not.deja_annonce                               .and.
     .            abs(pdelta-delta%dval) > 0.0000001*abs(pdelta) ) then
                  deja_annonce = .true.
                  print *,'More than one delta...',pdelta,delta%dval
               endif
               pdelta = delta%dval ! Sauver pour comparaison
            endif

         enddo

         exit

      enddo

      time_bnds_L = .false.
      do i=1,maxvar

         if (.not.infvar(i)%var_ok) cycle
         time_bnds_L = (time_bnds_L .or.infvar(i)%time_bnds_L)

      enddo

      ! Definir lablvl(1:nlev)

      do cvar=maxvar,1,-1

         if (.not.infvar(cvar)%var_ok) cycle

         nk = infvar(cvar)%zid
         
         if (infvar(cvar)%len(nk) /= nlev) cycle

         lablvl(1:nlev)=infvar(cvar)%niv(1:nlev)

         exit
            
      enddo

      ! Definir infvar(:)%range(1:2)

      iset=0

  100 call getset2( funit, dval,lev,nlev0, ibuf,maxpk,ok )

      if (ok) then

         iset=iset+1
         iiname=ibuf(3)
         write(name,'(a)') iiname

         nwds=ibuf(5)*ibuf(6)
         if(spec) nwds=nwds*2
         nwds=nwds*nlev0

         if (.not.(spec .or. fill_ccc_def)) then
            xmin = minval( dval(1:nwds) )
            xmax = maxval( dval(1:nwds) )
         else if (fill_ccc_def) then
                                ! Tenir compte des valeurs manquantes
            do istart=1,nwds
               if (abs( dval(istart)-fill_ccc ) >= fill_toler) exit
            enddo
            if (istart > nwds) then
               xmin = fill_ccc ; xmax = xmin
            else
               xmin = dval(istart) ; xmax = xmin
               do i=istart+1,nwds
                  if (abs( dval(i)-fill_ccc ) >= fill_toler) then
                     xmin = min( xmin,dval(i) )
                     xmax = max( xmax,dval(i) )
                  endif
               enddo
            endif
         endif

         do cvar=1,maxvar
            if(infvar(cvar)%name == name) then
               infvar(cvar)%range(1) = min( infvar(cvar)%range(1),xmin )
               infvar(cvar)%range(2) = max( infvar(cvar)%range(2),xmax )
            endif
         enddo

         goto 100

      else
         if (iset == 0) call                       xit('scanfile',-3)
      end if

      call precede(funit,-1)
*-----------------------------------------------------------------------
      include 'format.h'
 6001 format('Timedesc est "months since" et delta_t pas entier.')
 6002 format('Timedesc ',A," n'est pas encore supporte.")
*-----------------------------------------------------------------------
      end
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
      subroutine minmax (x,xmin,xmax,nwds)

      implicit none
      integer i,nwds
      real*8 x(nwds),xmin,xmax

      do  i=1,nwds,1
        xmin = min(x(i),xmin)
        xmax = max(x(i),xmax)
      enddo
*-----------------------------------------------------------------------
      end
