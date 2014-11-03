!
!netCDF Variables et variables coordonnees :

      TYPE information
       sequence
       real(8)       mult            ! fact. mult. conversion d''unite netCDF->CCCma
       real(8)       add             ! fact. add. conversion d''unite netCDF->CCCma
       character(80) name            ! nom  de la variable
       integer(4)    type            ! type de la variable
       integer(4)    ndim            ! nbre de dimension
       integer(4)    dimid(max_dims) ! id de la dimension correspondante
       integer(4)    nattr           ! nombre d''attribut
       character(52) padding         ! taille totale = 240 octets
      END TYPE information

      TYPE (information), dimension (:), pointer :: var,coord

!
!netCDF ID des Variables coordonnees :
 
      integer(4) xid      ! ID de la coordonnee x
      integer(4) yid      ! ID de la coordonnee y
      integer(4) zid      ! ID de la coordonnee z
      integer(4) tid      ! ID de la coordonnee t

      logical(4) no_time  ! switch pour cas intemporel

      character(80) lon   ! nom la variable de longitude dans fichier netCDF
      character(80) lat   ! nom la variable de latitude dans fichier netCDF

      common /id_coord_comm/ xid,yid,zid,tid,no_time,lon,lat


!
!netCDF file list :


      TYPE liste
       sequence
       character(80)   name            ! nom  de la variable
       integer(4)      type            ! type de la variable
       integer(4)      ndim            ! nbre de dimension
       integer(4)      dimid(max_dims) ! id de la dimension correspondante
       integer(4)      nattr           ! nombre d''attribut
       logical(4)      var_ok          ! variable logique
       character(64)   padding         ! taille totale = 240 octets
      END TYPE liste

      TYPE (liste), dimension (:), pointer :: list

      common /list_com/ var,coord,list

