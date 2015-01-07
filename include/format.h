
c	les formats des sous-programmes attribut_coord.f et enleve_bissextile.f
c	doivent concorder.

C7777 format(a,i5.4,2('-',i2.2),x,2(i2.2,':'),f3.1,'}')           !units
 7777 format(a,i5.4,2('-',i2.2),x,2(i2.2,':'),i2.2,'}')           !units
 7778 format(a,i5.4,2('-',i2.2),x,2(i2.2,':'),i2.2,'.0}')         !units
C7779 format(a,i5.4,2('-',i2.2),x,2(i2.2,':'),'0.0')              !units
 7779 format(a,i5.4,2('-',i2.2),x,2(i2.2,':'),'00')               !units
 7780 format(a,i5.4,2('-',i2.2),x,2(i2.2,':'),i2.2,'.0')          !units

 8888 format(a,i4.4,2('-',i2.2),x,2(i2.2,':'),i2.2,'}')           !delta_t
 8887 format(a12,i4.4,2(x,i2.2),x,2(i2.2,x))      

