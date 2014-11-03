*
*netCDF spectral :

      real*8, dimension(:), pointer :: mean          ! premier coefficient sepctral
      real*8, dimension(:), pointer :: add_offset    ! fact. add. de compression du coef. spec.
      real*8, dimension(:), pointer :: scale_fact    ! fact. mult. de compression du coef. spec.

      common /spec_com/ mean,add_offset,scale_fact
