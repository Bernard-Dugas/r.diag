*
*workmem.inc
*
      integer*1, dimension(:), pointer :: i1val
      integer*2, dimension(:), pointer :: i2val
      integer*4, dimension(:), pointer :: ival
      real*4,    dimension(:), pointer :: rval
      real*8,    dimension(:), pointer :: dval,rtime
*
      common /val_com/ i1val,i2val,ival,rval,dval,rtime
*

