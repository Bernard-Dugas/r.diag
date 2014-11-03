*netCDF Dimensions:

      integer xdid          ! id de la dimension x
      integer ydid          ! id de la dimension y
      integer zdid          ! id de la dimension z  
      integer timedid       ! id de la dimension t
      integer bndsdid       ! id de la dimension t
      integer numdid        ! id de la dimension spectral
      integer level2did
      integer unlimdimid    ! id de la dimension unlimited

      common /did_com/ xdid,ydid,zdid,timedid,numdid,level2did,
     .                   unlimdimid,bndsdid

      TYPE dimension
        sequence
        character*80 name           ! nom de la dimension
        integer      len            ! longueur de la dimension
        integer      duplic         ! repetition en bout de linge (shifted grid)
      END TYPE dimension

      TYPE (dimension), dimension (:), pointer :: dim

      common /dim_com/ dim
