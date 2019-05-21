# if defined (RDIAG_LICENCE)
!/* 
! * Copyright (C) 2017-2019  UQAM centre ESCER
! *
! * This software is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This software is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this software; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
# endif
 module averages_common   ! tables and table management routines
    use stats_signatures, only: time_mean_signature, variance_signature
    implicit none
    public
    
    ! Author Michel Valin, Centre ESCER/UQAM, 2017-2019

    ! Revision history :
    ! May 2019 by B.Dugas (ESCER/UQAM) - Adapted to the r.diag framework

    type :: field
      integer*8 :: date_lo, date_hi                ! earliest date, latest date collected
      real *8, dimension(:,:,:), pointer :: stats  ! (:,:,1) = sum, (:,:,2) = sum of squares
      integer :: ni, nj                            ! field dimensions
      integer :: nsamples, sample                  ! number of samples collected, duration of sample
 !    integer :: npas_min, npas_max                ! lowest and highest time step collected
      integer :: ip1, ip2, ip3, dateo, deet        ! from field's standard file parameters
      integer :: npas
      integer :: ig1, ig2, ig3, ig4                ! grid
      integer :: level2                            ! second level if 2 level data
      character(len=12) :: etiket                  ! from field's standard file parameters
      character(len=4)  :: nomvar                  ! from field's standard file parameters
      character(len=2)  :: typvar                  ! from field's standard file parameters
      character(len=1)  :: grtyp                   ! from field's standard file parameters
      logical :: special                           ! special record, ignore stats
    end type

    integer, parameter :: PAGE_SIZE  = 256   ! 256 fields per page in table (MUST BE A POWER OF 2)
    integer, parameter :: ENTRY_MASK = 255   ! PAGE_SIZE-1, mask to get slot number from index slot=and(ix,ENTRY_MASK)
    type :: page                             ! used to implement an array of pointers
      type(field), dimension(:), pointer :: p
    end type

    integer, parameter :: MAX_PAGES = 128   ! 128 pages at most (32768 records)
    integer, parameter :: PAGE_SHIFT = -8   ! right shift count to get page number from index
    type(page), dimension(:), pointer, save :: ptab => NULL()  ! statistics table, MAX_PAGES pages of size PAGE_SIZE

    integer, save :: next = -1             ! index of last valid record in statistics table
    integer, save :: verbose = 2           ! ERROR + WARNING
    logical, save :: variance = .true.     ! variance or standard deviation required [true unless -novar is used]
    integer, save :: fstdmean = 0          ! unit number for standard file into which averages are written
    integer, save :: fstdvar = 0           ! unit number for standard file into which variances/std deviatoins are written
    logical, save :: std_dev = .false.     ! output standard deviation rather than variance
    logical, save :: newtags = .false.     ! use new ip1/2/3 taggging style (not implemented yet)
    logical, save :: check_dateo = .false. ! check that all samples have the same date of origin (single model run)
    logical, save :: skip_npas0 = .true.   ! skip record if npas == 0 (default)
    logical, save :: weight_ip3 = .false.  ! use ip3 as a weight
    logical, save :: weight_time = .false. ! use time as a weight
    logical, save :: weight_abs = .false.  ! use a specific constant weight
    integer, save :: time_weight = 24      ! weight is in days, set to 1 for weight in hours
    logical, save :: strict = .false.      ! non strict mode by default
    logical, save :: ensemble = .false.    ! ensemble mode
    logical, save :: select_etiket = .false. ! etiket 1s a significant item if .true.

    character(len=4), dimension(1024), save :: specials
    integer, save :: nspecials=0
  contains

    subroutine create_page(n)  ! create and initialize page n that will contain PAGE_SIZE entries
      implicit none
      integer, intent(IN) :: n  ! page number
      integer :: i, status

      if(n < 0 .or. n >= MAX_PAGES) then
        print *,"FATAL: page number",n,' is invalid. it must be between 0 and',MAX_PAGES-1
        call quit('Averages_common - Create_page', 1 )
      endif
      if(associated(ptab(n)%p)) then
        print *,"FATAL: page number",n,' is already allocated'
        call quit('Averages_common - Create_page', 2 )
      endif
      allocate(ptab(n)%p(0:PAGE_SIZE-1),STAT=status)  ! allocate page
      if(status .ne. 0) then
        print *,"FATAL: page number",n,' cannot be allocated'
        call quit('Averages_common - Create_page', 3 )
      endif
      do i = 0, PAGE_SIZE-1
        nullify(ptab(n)%p(i)%stats)         ! nullify entries in page
      enddo
    end subroutine

    function create_page_entry(pg,slot,ni,nj) result(ix) ! create and initialize entry slot in page pg
      implicit none
      integer, intent(IN) :: pg, slot    ! page number, slot number
      integer, intent(IN) :: ni, nj      ! field dimensions (neede for stats)
      integer :: ix                      ! return value is index in table (ix < 0 indicates failure)
      integer :: status

      ix = -1
      if(slot < 0 .or. slot >= PAGE_SIZE ) then
        if(verbose > 0) print *,"ERROR: slot number",slot," not valid"
        return
      endif
      if(pg < 0 .or. pg >= MAX_PAGES ) then
        if(verbose > 0) print *,"ERROR: page number",pg," not valid"
        return
      endif

      if(.not. associated(ptab(pg)%p)) then  ! page does not exist, try to create it
        call create_page(pg)
      endif

      if(associated(ptab(pg)%p(slot)%stats)) return ! stats array exists, return with failure code

      ix = slot + PAGE_SIZE * pg                    ! build proper index number
      if(variance) then ! need sum and sum of squares
        allocate(ptab(pg)%p(slot)%stats(ni,nj,2),STAT=status)
      else              ! only need sum
        allocate(ptab(pg)%p(slot)%stats(ni,nj,1),STAT=status)
      endif

      if(verbose > 4) print *,'DEBUG: allocated stats - slot, ni,nj =',ix,ni,nj
      if(status .ne. 0) then
        print *,"FATAL: entry number",slot,'in page',pg,' cannot be allocated'
        call quit('Averages_common - create_page_entry', 1 )
      endif
!     initialize entry to EMPTY but READY
      ptab(pg)%p(slot)%stats = 0.0   ! set sum and sum of squares to zero
      ptab(pg)%p(slot)%ni = ni       ! dimensions of field
      ptab(pg)%p(slot)%nj = nj
      ptab(pg)%p(slot)%nsamples = 0  ! no samples yet
!       ptab(pg)%p(slot)%npas_min = 999999999
!       ptab(pg)%p(slot)%npas_max = -1
      ptab(pg)%p(slot)%dateo = -1
      ptab(pg)%p(slot)%date_lo = 0
      ptab(pg)%p(slot)%date_hi = 0
      ptab(pg)%p(slot)%sample = 0
      ptab(pg)%p(slot)%deet = -1
      ptab(pg)%p(slot)%npas = -1
      ptab(pg)%p(slot)%ip1 = -1
      ptab(pg)%p(slot)%level2 = -1
      ptab(pg)%p(slot)%ip2 = -1
      ptab(pg)%p(slot)%ip3 = -1
      ptab(pg)%p(slot)%ig1 = -1
      ptab(pg)%p(slot)%ig2 = -1
      ptab(pg)%p(slot)%ig3 = -1
      ptab(pg)%p(slot)%ig4 = -1
      ptab(pg)%p(slot)%etiket = ""
      ptab(pg)%p(slot)%nomvar = ""
      ptab(pg)%p(slot)%typvar = ""
      ptab(pg)%p(slot)%grtyp = ""
      ptab(pg)%p(slot)%special = .false.   ! by default not a special record
      return
    end function create_page_entry

    function new_page_entry(ni,nj) result(ix) ! create a new entry in tables
      implicit none
      integer, intent(IN) :: ni, nj
      integer :: ix
      integer :: pg, slot

      next = next + 1             ! next valid index
      pg = next / PAGE_SIZE       ! page number
      slot = mod(next,PAGE_SIZE)  ! slot number
      ix = create_page_entry(pg,slot,ni,nj)
      if(ix < 0) then
        print *,"FATAL: error creating entry",slot,' in page',pg
        call quit('Averages_common - New_page_entry', 1 )
      endif
      return
    end function new_page_entry

    function date_stamp_64(date_stamp) result(date_64)  ! compute date_64 from CMC date_stamp
      implicit none
      integer, intent(IN) :: date_stamp
      integer*8 :: date_64  ! in seconds

      integer :: t1, t2 
      integer year, month, day, hour, minute, second, julian

      call newdate(date_stamp,t1,t2,-3)  ! date stamp to t1(YYYYMMDD), t2(HHMMSShh)

      year = t1 / 10000
      t1 = t1 - (10000*year)
      month = t1 / 100
      day = t1 - (month*100)

      hour = t2 / 1000000
      t2 = t2 - (1000000*hour)
      minute = t2 / 10000
      t2 = t2 - (10000*minute)
      second = t2 / 100

      call jdatec(julian,year,month,day)
      date_64 = julian
!       date_64 = (date_64*86400) + (hour*3600) + (minute*60) + second     ! 86400 seconds in a day
      date_64 = date_64*86400
      date_64 = date_64 + hour*3600
      date_64 = date_64 + minute*60
      date_64 = date_64 + second
      return
    end function date_stamp_64

    function date_stamp_32(date_64) result(date_stamp)  ! convert date_64 to CMC datestamp
      implicit none
      integer*8, intent(IN) :: date_64  ! in seconds
      integer :: date_stamp

      integer :: t1, t2, t0  ! used later when implementing this function
      integer year, month, day, hour, minute, second, hms, ymd
      integer*8 :: ymd8

      ymd8 = date_64 / 86400    ! julian day
      ymd = ymd8
      call datec(ymd,year,month,day)
      t1 = year*10000 + month*100 + day   ! YYYYMMDD
      hms = date_64 -(ymd*86400)
      hour = hms / 3600
      hms = hms - (hour*3600)
      minute = hms / 60
      second = hms - (minute*60)
      t2 = hour*1000000 + minute*10000 + second*100  ! HHMMSS00
      call newdate(t0,t1,t2,3)  ! t1(YYYYMMDD), t2(HHMMSShh) to  date stamp
      date_stamp = t0
      return
    end function date_stamp_32
!
!   process record read from one of the input standard files
!
    function process_entry(z,ni,nj,ip1,ip2,ip3,dateo,deet,npas,etiket,nomvar,typvar,grtyp,ig1,ig2,ig3,ig4) result(ix)
      implicit none
      integer, intent(IN) :: ni, nj, ip1, ip2, ip3, npas, dateo, deet
      integer, intent(IN) :: ig1, ig2, ig3, ig4
      character(len=*), intent(IN) :: etiket, nomvar, typvar, grtyp
      real, dimension(ni,nj), intent(IN) :: z
      integer :: ix

      integer :: i, j, pg, slot, sample, it1, it2, it3, level2
      integer*8 :: date_lo, date_hi, dnp
      type(field), pointer :: p
      real :: weight, r1, r2, r3
      logical :: is_special
      character(len=128) :: string

      ix = -1
      level2 = -1  ! a priori, one level data
      call diag_convip_plus( ip1, r1, it1, -1, string, .false. )    ! convert ip1
      if(ip3 == 0 .and. ip2 < 1000000) then    ! special "old style coding" for time in hours
        it2 = 10
        r2 = ip2
      else
         call diag_convip_plus( ip2, r2, it2, -1, string, .false. )  ! convert ip2
      endif
      call diag_convip_plus( ip3, r3, it3, -1, string, .false. )    ! convert ip3

      if( (it1 == it3) .and. (ip1 > 0) .and. (ip3 > 0) ) level2 = ip3   ! 2 level data
      if( (it1 == it2) .and. (ip1 > 0) .and. (ip2 > 0) ) level2 = ip2   ! 2 level data

      dnp = npas
      dnp = dnp * deet
      is_special = any(nomvar == specials(1:nspecials))
      weight = 1.0
      sample = 0  ! "instantaneous" sample
      if((it2 == 10) .and. (it3 == 10))then  ! both ip2 and ip3 are time tags
        sample = 3600*ABS(r2 - r3)           ! sample interval in seconds
        if(verbose > 4) print *,'DEBUG: sample =',sample
      endif
      if(is_special .or. ensemble) then   !  special names or ensemble mode
        date_lo = 0
        date_hi = 0
      else                  ! regular data / averages
        if( weight_ip3 .or. typvar(2:2) .eq. time_mean_signature ) then ! weight is IP3 or number of samples in IP3
          i = ip3
          if(ishft(ip3,-24) == 15) i = iand(ip3,Z'000FFFFF')   ! keep lower 20 bits (type 15)
          weight = max(1,i)     ! number of samples
        endif
        if(weight_time) then                                 ! time weight
          weight = (dnp) / (3600.0 * time_weight)
        endif
        if(weight_abs) then
          weight = time_weight                  ! explicit weight
        endif
        date_lo = date_stamp_64(dateo)      ! compute 64 bit date_lo from dateo
        date_hi = date_lo + dnp             ! date of validity of sample
        if(weight == 1.0) then
          date_lo = date_hi                 ! one date kept (date of validity)
        else
          sample = (date_hi - date_lo)      ! interval in seconds between the 2 dates
          if(weight_ip3) sample = nint( (date_hi - date_lo) / (weight-1.0) )
        endif
      endif
      do i = 0 , next               ! do we have an entry that matches this record's parameters
        slot = iand(i,ENTRY_MASK)   ! slot from index
        pg = ishft(i,PAGE_SHIFT)    ! page from index
        p => ptab(pg)%p(slot)       ! pointer to data
        if(p%ip1 .ne. ip1 .or. p%level2 .ne. level2) cycle    ! not same level(s)
        if(p%ni .ne. ni .or. p%nj .ne. nj) cycle              ! not same dimensions
        if(trim(p%nomvar) .ne. trim(nomvar)) cycle            ! not same name
        if(trim(p%grtyp) .ne. trim(grtyp)) cycle              ! not same grid type
        if((p%ig1 .ne. ig1) .or. (p%ig2 .ne. ig2) .or. (p%ig3 .ne. ig3) .or. (p%ig4 .ne. ig4) ) cycle  ! not same grid
        if(trim(p%etiket) .ne. trim(etiket) .and. select_etiket) cycle            ! not same experiment
        if((p%dateo .ne. dateo) .and. check_dateo) cycle      ! dateo verification is optional
        if(is_special .or. ensemble)then                      ! ensemble mode : records must have same ip1/ip2/ip3
          if(p%ip2 .ne. ip2 .or. p%ip3 .ne. ip3) cycle        ! special records must have same ip1/ip2/ip3
        endif
        if(p%typvar(2:2) .eq. time_mean_signature) p%typvar = trim(typvar)  ! force typvar into p%typvar if average
        if(p%typvar(1:1) .ne. typvar(1:1)) cycle               ! check first character of typvar
        if(sample .ne. p%sample .and. weight == 1.0) then
           if(verbose > 1) print *,'WARNING: sample interval mismatch, got',sample,' expected',p%sample
           if(strict) call quit('Averages_common - Process_entry', 1 )
        endif
        ix = i                      ! a matching entry has been found
        if(verbose > 4) print *,'DEBUG: found entry, previous samples',ix,p%nsamples
        exit
      enddo

      if(ix == -1) then              ! not found, make a new entry
        ix = new_page_entry(ni,nj)
        slot = iand(ix,ENTRY_MASK)   ! slot from index
        pg = ishft(ix,PAGE_SHIFT)    ! page from index
        p => ptab(pg)%p(slot)        ! pointer to data
        p%etiket = trim(etiket)      ! initialize entry using parameters from record
        p%nomvar = trim(nomvar)
        p%typvar = trim(typvar)
        p%grtyp = trim(grtyp)
        p%dateo = dateo
        p%date_lo = date_lo
        p%date_hi = date_hi
        p%sample = sample
        p%deet = deet
        p%npas = npas
        p%ip1 = ip1
        p%level2 = level2
        p%ip2 = ip2
        p%ip3 = ip3
        p%ig1 = ig1
        p%ig2 = ig2
        p%ig3 = ig3
        p%ig4 = ig4
        p%ni = ni
        p%nj = nj
        p%nsamples = 0
        p%special = is_special
!         print *,'DEBUG: dateo, deet,npas =',dateo, deet,npas,date_lo
      endif
      p%nsamples = p%nsamples + weight ! add weight to number of samples
!       if(p%nsamples == 48) print *,'DEBUG: dateo, deet,npas =',dateo, deet,npas,p%date_lo,date_hi
      if(is_special) return  ! no stats for special records, just add one to sample count
!       p%npas_max = max(p%npas_max,npas)
!       p%npas_min = min(p%npas_min,npas)
      p%date_lo = min(p%date_lo , date_lo) ! update earliest/latest date
      p%date_hi = max(p%date_hi , date_hi)
      do j = 1 , nj
      do i = 1 , ni
         p%stats(i,j,1) = p%stats(i,j,1) + z(i,j)*weight               ! update sum
         if(variance) p%stats(i,j,2) = p%stats(i,j,2) + z(i,j)*dble(z(i,j))*weight  ! update sum of squares if necessary
      enddo
      enddo
    end function process_entry
 end module averages_common
