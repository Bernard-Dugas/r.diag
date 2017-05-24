!
!netCDF Variables et variables coordonnees :
!
      real*8, dimension(:,:), pointer :: variable
      real*8, dimension(:,:), pointer :: dcoordonne,time_bnds
!
      common /var_com/ dcoordonne,variable,time_bnds
