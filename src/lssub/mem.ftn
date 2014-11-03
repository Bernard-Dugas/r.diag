#     if !defined (valeur_nmax)
#         define   valeur_nmax 100000
#     endif
#     if !defined (valeur_mmax)
#         define   valeur_mmax 5000
#     endif
#     if defined (RDIAG_LICENCE)
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
#     endif
! code extrait de:
! maximum entropy method from numerical recipes.
! Chapter 13, Section 6 (Linear prediction,
! Linear predictive coding)
!
      REAL FUNCTION evlmem(fdt,d,m,xms)

      IMPLICIT none

      INTEGER  m
      REAL     fdt,xms,d(m)

      INTEGER  i
      REAL     sumi,sumr
      REAL*8   theta,wi,wpi,wpr,wr,wtemp

      theta=6.28318530717959d0*fdt
      wpr=cos(theta)
      wpi=sin(theta)
      wr=1.d0
      wi=0.d0
      sumr=1.
      sumi=0.
      do 11 i=1,m
        wtemp=wr
        wr=wr*wpr-wi*wpi
        wi=wi*wpr+wtemp*wpi
        sumr=sumr-d(i)*sngl(wr)
        sumi=sumi-d(i)*sngl(wi)
11    continue
      evlmem=xms/(sumr**2+sumi**2)

      return
      END

      SUBROUTINE fixrts(d,m)

      IMPLICIT    none

      INTEGER     MMAX
      PARAMETER ( MMAX = valeur_mmax )

      INTEGER     m
      REAL        d(m)
!U    USES zroots
      INTEGER     i,j
      LOGICAL     polish
      COMPLEX     a(MMAX),roots(MMAX)

      EXTERNAL    zroots

      a(m+1)=cmplx(1.,0.)
      do 11 j=m,1,-1
        a(j)=cmplx(-d(m+1-j),0.)
11    continue
      polish=.true.
      call zroots(a,m,roots,polish)
      do 12 j=1,m
        if(abs(roots(j)).gt.1.)then
          roots(j)=1./conjg(roots(j))
        endif
12    continue
      a(1)=-roots(1)
      a(2)=cmplx(1.,0.)
      do 14 j=2,m
        a(j+1)=cmplx(1.,0.)
        do 13 i=j,2,-1
          a(i)=a(i-1)-roots(j)*a(i)
13      continue
        a(1)=-roots(j)*a(1)
14    continue
      do 15 j=1,m
        d(m+1-j)=-real(a(j))
15    continue

      return
      END

      SUBROUTINE laguer(a,m,x,its)

      IMPLICIT none

      REAL       EPSS
      INTEGER    MAXIT,MR,MT
      PARAMETER (EPSS=2.e-7,MR=8,MT=10,MAXIT=MT*MR)

      INTEGER    m,its
      COMPLEX    a(m+1),x

      INTEGER    iter,j
      REAL       abx,abp,abm,err,frac(MR)
      COMPLEX    dx,x1,b,d,f,g,h,sq,gp,gm,g2

      SAVE       frac
      DATA       frac /.5,.25,.75,.13,.38,.62,.88,1./

      do 12 iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=abs(b)
        d=cmplx(0.,0.)
        f=cmplx(0.,0.)
        abx=abs(x)
        do 11 j=m,1,-1
          f=x*f+d
          d=x*d+b
          b=x*b+a(j)
          err=abs(b)+abx*err
11      continue
        err=EPSS*err
        if(abs(b).le.err) then
          return
        else
          g=d/b
          g2=g*g
          h=g2-2.*f/b
          sq=sqrt((m-1)*(m*h-g2))
          gp=g+sq
          gm=g-sq
          abp=abs(gp)
          abm=abs(gm)
          if(abp.lt.abm) gp=gm
          if (max(abp,abm).gt.0.) then
            dx=m/gp
          else
            dx=exp(cmplx(log(1.+abx),real(iter)))
          endif
        endif
        x1=x-dx
        if(x.eq.x1)return
        if (mod(iter,MT).ne.0) then
          x=x1
        else
          x=x-dx*frac(iter/MT)
        endif
12    continue
      
      print *, 'too many iterations in laguer'

      return
      END

      SUBROUTINE memcof(data,n,m,xms,d)

      IMPLICIT    none

      INTEGER     MMAX
      PARAMETER ( MMAX = valeur_mmax )
      INTEGER     NMAX
      PARAMETER ( NMAX = valeur_nmax )

      INTEGER     m,n
      REAL        xms,d(m),data(n)

      INTEGER     i,j,k
      REAL        denom,p,pneum,wk1(NMAX),wk2(NMAX),wkm(MMAX)

      if (m.gt.MMAX.or.n.gt.NMAX) then
         print *, 'workspace NMAX too small in memcof/mem.ftn'
         call xit(' memcof ', -1 )
      endif

      p=0.
      do 11 j=1,n
        p=p+data(j)**2
11    continue
      xms=p/n
      wk1(1)=data(1)
      wk2(n-1)=data(n)
      do 12 j=2,n-1
        wk1(j)=data(j)
        wk2(j-1)=data(j)
12    continue
      do 17 k=1,m
        pneum=0.
        denom=0.
        do 13 j=1,n-k
          pneum=pneum+wk1(j)*wk2(j)
          denom=denom+wk1(j)**2+wk2(j)**2
13      continue
        d(k)=2.*pneum/denom
        xms=xms*(1.-d(k)**2)
        do 14 i=1,k-1
          d(i)=wkm(i)-d(k)*wkm(k-i)
14      continue
        if(k.eq.m)return
        do 15 i=1,k
          wkm(i)=d(i)
15      continue
        do 16 j=1,n-k-1
          wk1(j)=wk1(j)-wkm(k)*wk2(j)
          wk2(j)=wk2(j+1)-wkm(k)*wk1(j+1)
16      continue
17    continue
      print *, 'never get here in memcof'

      END

      SUBROUTINE predic(data,ndata,d,m,future,nfut)

      IMPLICIT    none

      INTEGER     MMAX
      PARAMETER ( MMAX = valeur_mmax )

      INTEGER     ndata,nfut,m
      REAL        d(m),data(ndata),future(nfut)

      INTEGER     j,k
      REAL        discrp,sum,reg(MMAX)

      do 11 j=1,m
        reg(j)=data(ndata+1-j)
11    continue
      do 14 j=1,nfut
        discrp=0.
        sum=discrp
        do 12 k=1,m
          sum=sum+d(k)*reg(k)
12      continue
        do 13 k=m,2,-1
          reg(k)=reg(k-1)
13      continue
        reg(1)=sum
        future(j)=sum
14    continue

      return
      END

      SUBROUTINE zroots(a,m,roots,polish)

      IMPLICIT   none

      REAL        EPS
      INTEGER     MAXM
      INTEGER     MMAX

      PARAMETER ( EPS  = 1.e-6 )
      PARAMETER ( MMAX = valeur_mmax )
      PARAMETER ( MAXM = MMAX + 1 )

      INTEGER    m
      COMPLEX    a(m+1),roots(m)
      LOGICAL    polish

CU    USES laguer
      INTEGER    i,j,jj,its
      COMPLEX    ad(MAXM),x,b,c

      EXTERNAL   laguer

      do 11 j=1,m+1
        ad(j)=a(j)
11    continue
      do 13 j=m,1,-1
        x=cmplx(0.,0.)
        call laguer(ad,j,x,its)
        if(abs(aimag(x)).le.2.*EPS**2*abs(real(x))) x=cmplx(real(x),0.)
        roots(j)=x
        b=ad(j+1)
        do 12 jj=j,1,-1
          c=ad(jj)
          ad(jj)=b
          b=x*b+c
12      continue
13    continue
      if (polish) then
        do 14 j=1,m
          call laguer(a,m,roots(j),its)
14      continue
      endif
      do 16 j=2,m
        x=roots(j)
        do 15 i=j-1,1,-1
          if(real(roots(i)).le.real(x))goto 10
          roots(i+1)=roots(i)
15      continue
        i=0
10      roots(i+1)=x
16    continue

      return
      END
