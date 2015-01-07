# if defined (RDIAG_LICENCE)
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
# endif
  SUBROUTINE FWGHTS(L0,K,W,LOPASS,SUMW) 
  
      IMPLICIT  none

      INTEGER   K
      LOGICAL   LOPASS
      REAL      L0,W(-K:K),SUMW

!     * JUILLET 18/85 - F. ZWIERS, N.E. SARGENT, B.DUGAS. 
  
!     * COMPUTE WEIGHTS FOR:  
  
!     * 1) FWGHTS ====================================== 
!     * IDEAL LOW PASS FILTER WITH CUTOFF FREQUENCY L0
!     * USING HANNING TAPER WITH K WEIGHTS, 
  
!     * 2) FWGHTSS ====================================== 
!     * PRODUCT OF TWO SEQUENTIALLY APPLIED SYMMETRIC AVERAGE 
!     * FILTERS WITH BASE CUTOFF FREQUENCY L0.

!     Local variables

      INTEGER   I,J,JR,JRP
      REAL(8)   W8(2*K+1),SUMW8,PI,PII,FK,FI,AK,HK,BK,F
      
!---------------------------------------------------------------------------- 
  
      PI  = 4.0*ATAN(1.D0)
      PII = 1./PI

!     * IDEAL FILTER. 
  
      FK      = K
      W8(K+1) = L0/PI
      SUMW8   = W8(K+1) 
      DO  I=1,K
          FI        = I
          AK        = SIN(FI*L0)*PII/FI
          HK        = 0.5*(1.0+COS(PI*FI/FK))
          BK        = HK*AK
          W8(K+1+I) = BK 
          W8(K+1-I) = BK 
          SUMW8     = SUMW8+BK+BK 
      END DO
  
      IF (.NOT.LOPASS) SUMW8 = -SUMW8

      SUMW = SUMW8

      DO  I=-K,K
          W(I) = W8(K+1+I)/SUMW8
      END DO

      IF (.NOT.LOPASS) W(0) = W(0)+1

      RETURN

      ENTRY FWGHTSS(L0,K,W,LOPASS,SUMW) 
  
!     * SEQUENTIAL MOVING AVERAGE FILTER. 
  
      PI  = 4.0*ATAN(1.D0)

      F   = L0/(2.0*PI) 
      JR  = INT(1.0/(2.0*F))
      JRP = 2*JR/3
      K   = JR+JRP

      W = 0.0
  
      DO  I=-JR,JR 
         DO  J=-JRP,JRP 
            W(I+J) = W(I+J) + 1.0/FLOAT((2*JR+1)*(2*JRP+1))
         END DO
      END DO
  
      IF (.NOT.LOPASS)                                     THEN 

         W = -W ; W(0) = W(0)+1 
  
      ENDIF 
  
      RETURN

  END SUBROUTINE FWGHTS
  SUBROUTINE GAINSQ(L0,K,LAMBDA,B,NORMAL) 
  
!     * AOUT 02/85 - F. ZWIERS, B.DUGAS.
  
!     * COMPUTE AND DISPLAY THE SQUARED GAIN OF THE IDEAL LOW PASS
!     * FILTER WITH WEIGHTS TAPERED BY HANNING TAPER
!     * 
!     * L0 IS THE CUTOFF FREQUENCY OF THE IDAL LOW-PASS FILTER
!     * K IS THE TRUNCATION POINT OF THE FILTER WEIGHTS 
!     * 
!     * LAMBDA RETURNS THE FREQUENCIES AT WHICH THE SQUARED GAIN
!     * IS EVALUATED
!     * B RETURNS THE CORRESPONDING VALUES OF THE SQUARED GAIN
!     * 
!     * NOTE
!     *------------------------------------------------------------------ 
!     * 
!     * ROUTINE CADRE IS USED TO INTEGRATE THE TRANSFER FUNCTION
!     * OF THE IDEAL LOW PASS FILTER TRUNCATED AFTER K TERMS.
!     * 
!     *-----------------------------------------------------------------
  
      IMPLICIT none

      INTEGER  K
      REAL     L0,LAMBDA(*),B(*),NORMAL

!     LOCAL VARIABLES.

      INTEGER, PARAMETER :: NFREQ = 199

      INTEGER  KAY,NF2,I,IER(3),MAXIER
      REAL(8)  AERR,RERR,ERR,PIOV2,PI,TWOPI,D,PIBK
      REAL(8)  F0,F1,G0,G1,H0,H1
      REAL     RANGE(4)

      REAL(8), EXTERNAL :: CADREF

      MAXIER = 0
  
      AERR   = 1.E-6
      RERR   = 1.E-6

      PIOV2  = 2.0*ATAN(1.D0)
      PI     = 2.0*PIOV2
      TWOPI  = 2.0*PI

      IF (L0 < PIOV2)                                          THEN 
          D = 2.0*L0/DBLE(NFREQ-1) 
      ELSE
          D = 2.0*(PI-L0)/DBLE(NFREQ-1)
      ENDIF 
  
      NF2 = (NFREQ-1)/2 
      DO  I=-NF2,NF2
          LAMBDA(I+NF2+1) = L0+DBLE(I)*D 
      END DO
  
      PIBK = PI/DBLE(K)
      DO  I=1,NFREQ 
          F0        = LAMBDA(I)-L0
          F1        = LAMBDA(I)+L0
          G0        = F0-PIBK 
          G1        = F1-PIBK 
          H0        = F0+PIBK 
          H1        = F1+PIBK 
          B(I)      = (0.50*CADREF(F0,F1,K,AERR,RERR,ERR,IER(1)) &
                    +  0.25*CADREF(G0,G1,K,AERR,RERR,ERR,IER(2)) &
                    +  0.25*CADREF(H0,H1,K,AERR,RERR,ERR,IER(3)))/TWOPI 
          B(I)      = (B(I)/NORMAL)**2
          LAMBDA(I) = LAMBDA(I)/TWOPI 
          MAXIER = MAX( MAXIER, IER(1),IER(2),IER(3) )
      END DO
      IF (MAXIER > 2) CALL            XIT(' Gainsq ',-1 )
  
      WRITE(6,6030) 
      WRITE(6,6040) (LAMBDA(I),B(I),I=1,NFREQ)
  
      RANGE(1) = LAMBDA(1)
      RANGE(2) = LAMBDA(NFREQ)
      RANGE(3) = -0.4 
      RANGE(4) =  1.6 
  
      RETURN
!---------------------------------------------------------------------- 
 6030 FORMAT('0SQUARED GAIN OF THE LOW PASS FILTER:'/  &
             5(7X,'F',5X,'ABS(B(F))**2')/5('  -----------------------'))
 6040 FORMAT(5(1X,2E12.4)) 

  END SUBROUTINE GAINSQ 
  SUBROUTINE  LISSXY()

      WRITE(0,'(A)')'SUBROUTINE LISSXY2 REPLACES LISSXY' 
      STOP 'In outdated call to LISSXY...'

  END SUBROUTINE LISSXY
  SUBROUTINE LISSXY2( G, LISSE,NLON,NLAT )

!**    September 2004 - B.Dugas

!**    SEQUENTIALLY APPLY DIGITAL FILTERS OF ORDER LISSE
!**    IN THE X AND Y DIRECTIONS ON THE FIELD G(NLON,NLAT).
!**    W IS A WORK ARRAY THAT CAN STILL HOLD PREVIOUSLY
!**    CALCULATED FILTER WEIGHTS.

      IMPLICIT none

!**    ARGUMENTS.

      INTEGER  LISSE,NLON,NLAT
      REAL     G(NLON,NLAT)

!**    LOCAL DECLARATIONS.

      REAL     L0,PI,K0,SUMW
      INTEGER  I,J,K,MAXSIZ,NL
      REAL,    DIMENSION (:), ALLOCATABLE :: HOLDX,HOLDY
      REAL,    DIMENSION(:,:),ALLOCATABLE,SAVE :: W

      INTEGER  OLD_LISSE
      SAVE     OLD_LISSE
      DATA     OLD_LISSE / -1 /

!----------------------------------------------------------------------------
      MAXSIZ = MIN( NLON,NLAT )
      IF (LISSE >= MAXSIZ) RETURN

      ALLOCATE( HOLDX(NLON) , HOLDY(NLAT) )

!**    CALCULATE THE FILTER WEIGHTS.

      IF (LISSE /= OLD_LISSE) THEN

         IF (OLD_LISSE.NE.-1) DEALLOCATE( W )

         ALLOCATE( W(-LISSE:LISSE,LISSE) )

         W = 0.

         PI = 4.0*ATAN(1.0D0)
         L0 = PI/LISSE

         DO K=2,LISSE
            CALL FWGHTS( L0,K, W(-K,K),.TRUE.,SUMW )
         END DO

         OLD_LISSE = LISSE

      END IF

!**    SMOOTH IN X.

      DO J=1,NLAT
         HOLDX      = 0.0
         HOLDX(1) = G(1,J)
         DO I=2,LISSE-1
            DO K=-(I-1),+(I-1)
               HOLDX(I) = HOLDX(I)+W(K,I)*G(I+K,J)
            ENDDO
         ENDDO
         DO I=LISSE,NLON-LISSE+1
            DO K=-LISSE+1,LISSE-1
               HOLDX(I) = HOLDX(I)+W(K,LISSE)*G(I+K,J)
            ENDDO
         ENDDO
         DO I=NLON-LISSE+2,NLON-1
            DO K=-(NLON-I),+(NLON-I)
               HOLDX(I) = HOLDX(I)+W(K,NLON-I+1)*G(I+K,J)
            ENDDO
         ENDDO
         HOLDX(NLON) = G(NLON,J)
         G(:,J) = HOLDX
      ENDDO

!**    SMOOTH IN Y.

      DO I=1,NLON
         HOLDY    =   0.0
         HOLDY(1) = G(I,1)
         DO J=2,LISSE-1
            DO K=-(J-1),+(J-1)
               HOLDY(J) = HOLDY(J)+W(K,J)*G(I,J+K)
            ENDDO
         ENDDO
         DO J=LISSE,NLAT-LISSE+1
            DO K=-LISSE+1,LISSE-1
               HOLDY(J) = HOLDY(J)+W(K,LISSE)*G(I,J+K)
            ENDDO
         ENDDO
         DO J=NLAT-LISSE+2,NLAT-1
            DO K=-(NLAT-J),+(NLAT-J)
               HOLDY(J) = HOLDY(J)+W(K,NLAT-J+1)*G(I,J+K)
            ENDDO
         ENDDO
         HOLDY(NLAT) = G(I,NLAT)
         G(I,:) = HOLDY
      ENDDO

      DEALLOCATE( HOLDX , HOLDY )
!----------------------------------------------------------------------------

      RETURN

  END SUBROUTINE LISSXY2
  REAL(8) FUNCTION CADREF ( a, b, k, abserr, relerr, error, ind )
!
!***********************************************************************
!
!! CADREF estimates the integral of F(X) from A to B.
!  On Output, returns the approximate value of the integral.
!
!  Discussion:
!
!    CADREF is the Cautious Adaptive Romberg Extrapolator.
!
!    Obtained from http://orion.math.iastate.edu/burkardt/,
!    a page created by John Burkardt, from 2000 to 2002
!    at the Math department of Iowa State University
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    Carl DeBoor and J R Rice,
!    CADRE: An algorithm for numerical quadrature,
!    Mathematic Software, pages 417-449,
!    Academic Press, New York, 1971.
!
!  Modified:
!
!    30 October 2000
!     6 February 2014 - B. Dugas (UQAM) :
!            - real => real(8) ; subroutine => function
!            - Integrates F, the transfer function of the
!              ideal low pass filter truncated at K terms
!            - Renamed to CADREF
!
!  Parameters:
!
!    Input, real A, the lower limit of integration.
!
!    Input, real B, the upper limit of integration.
!
!    Input, integer K, the truncation point of the filter weights 
!
!    Input, real ABSERR, the absolute error tolerance.
!
!    Input, real RELERR, the relative error tolerance.
!
!    Output, real ERROR, an estimate of the absolute error.
!
!    Output, integer IND, reliability indicator.
!    If IND <= 2, RESULT is very reliable.  Higher values of
!    IND indicate less reliable values of RESULT.
!
  implicit none
!
  integer, parameter :: mxstge = 30
  integer, parameter :: maxtbl = 10
  integer, parameter :: maxts = 2049
!
  real(8) a
  real(8) abserr
  real(8) ait(maxtbl)
  logical aitken
  real(8), parameter :: aitlow = 1.1_8
  real(8), parameter :: aittol = 0.1_8
  real(8) astep
  real(8) b
  real(8) beg
  real(8) begin(mxstge)
  real(8) bma
  real(8) curest
  real(8) dif(maxtbl)
  real(8) diff
  real(8) end
  real(8) ergoal
  real(8) erra
  real(8) errer
  real(8) error
  real(8) errr
  real(8) est(mxstge)
  real(8) fbeg
  real(8) fbeg2
  real(8) fend
  real(8) fextm1
  real(8) fextrp
  real(8) finis(mxstge)
  real(8) fn
  real(8) fnsize
  logical h2conv
  real(8) h2next
  real(8) h2tfex
  real(8), parameter :: h2tol = 0.15_8
  real(8) hovn
  integer i
  integer ibeg
  integer ibegs(mxstge)
  integer iend
  integer ii
  integer iii
  integer ind
  integer istage
  integer istep
  integer istep2
  integer it
  integer k
  integer l
  integer lm1
  integer n
  integer n2
  integer nnleft
  real(8) prever
  real(8) r(maxtbl)
  logical reglar
  logical reglsv(mxstge)
  real(8) relerr
  real(8) result
  logical right
  real(8) rn(4)
  real(8) rnderr
  real(8) sing
  real(8) singnx
  real(8) slope
  real(8) stage
  real(8) step
  real(8) stepmn
  real(8) sum1
  real(8) sumabs
  real(8) t(maxtbl,maxtbl)
  real(8) tabs
  real(8) tabtlm
  real(8), parameter :: tljump = 0.01_8
  real(8) ts(2049)
  real(8) vint

! FUNC, THE TRANSFER FUNCTION OF THE IDEAL
! LOW PASS FILTER TRUNCATED AT K TERMS 
  
  real(8) FUNC,LAMBDA
  FUNC(LAMBDA) = MERGE( 2.0*DBLE( K )+1.0,                    &
                        SIN( (2.0*DBLE( K )+1.0)*LAMBDA/2.0 ) &
               /        SIN( LAMBDA/2.0 ),                    &
                        LAMBDA == 0 )

  if ( a == b ) then
    cadref = 0.0
    return
  end if
 
  begin(1:mxstge) = 0.0
  est(1:mxstge) = 0.0
  finis(1:mxstge) = 0.0
  ibegs(1:mxstge) = 0
  reglsv(1:mxstge) = .false.
 
  vint = 0.0
 
  rn(1:4) = (/ 0.7142005_8, 0.3466282_8, 0.8437510_8, 0.1263305_8 /)
 
  rnderr = epsilon ( rnderr )
  cadref = 0.0
  error = 0.0
  ind = 1
  bma = abs ( b - a )
  errr = min ( 0.1_8, max ( abs ( relerr ), 10.0*rnderr) )
  erra = abs ( abserr )
  stepmn = max ( bma / 2**mxstge, max ( bma, abs ( a ), abs ( b ) ) * rnderr )
  stage = 0.5
  istage = 1
  curest = 0.0
  fnsize = 0.0
  prever = 0.0
  reglar = .false.
  beg = a
  fbeg = func(beg) / 2.0
  ts(1) = fbeg
  ibeg = 1
  end = b
  fend = func(end) / 2.0
  ts(2) = fend
  iend = 2
 
10 continue
 
  right = .false.
 
20 continue

  step = end - beg
  astep = abs ( step )
 
  if ( astep < stepmn ) then
    ind = 5
    cadref = curest + vint
    return
  end if
 
  t(1,1) = fbeg+fend
  tabs = abs ( fbeg ) + abs ( fend )
  l = 1
  n = 1
  h2conv = .false.
  aitken = .false.
  go to 40
 
30 continue
 
40 continue
 
  lm1 = l
  l = l+1
  n2 = n*2
  fn = n2
  istep = (iend-ibeg)/n

  if ( istep > 1 ) then
    go to 60
  end if

  ii = iend
  iend = iend + n

  if ( iend > maxts ) then
    go to 440
  end if

  hovn = step / fn
 
  iii = iend
  do i = 1, n2, 2
    ts(iii) = ts(ii)
    ts(iii-1) = func(end-i*hovn)
    iii = iii-2
    ii = ii-1
  end do
 
  istep = 2
 
60 continue
 
  istep2 = ibeg+istep/2
 
  sum1 = 0.0
  sumabs = 0.0
  do i = istep2, iend, istep
    sum1 = sum1 + ts(i)
    sumabs = sumabs + abs ( ts(i) )
  end do
 
  t(l,1) = t(l-1,1) / 2.0 + sum1 / fn
  tabs = tabs / 2.0 + sumabs / fn
 
  n = n2
  it = 1
  vint = step * t(l,1)
  tabtlm = tabs * rnderr
  fnsize = max ( fnsize, abs ( t(l,1) ) )
  ergoal = max ( astep * rnderr * fnsize, &
    stage * max ( erra , errr * abs ( curest+vint ) ) )
  fextrp = 1.0
  do i = 1, lm1
    fextrp = fextrp * 4.0
    t(i,l) = t(l,i) - t(l-1,i)
    t(l,i+1) = t(l,i) + t(i,l) / ( fextrp - 1.0 )
  end do
 
  errer = astep * abs ( t(1,l) )
  if ( l > 2 ) go to 90
  if ( abs ( t(1,2) ) <= tabtlm ) go to 290
  go to 40
 
90 continue
 
  do i = 2, lm1

    if ( abs ( t(i-1,l) ) > tabtlm ) then
      diff = t(i-1,lm1) / t(i-1,l)
    else
      diff = 0.0
    end if

    t(i-1,lm1) = diff

  end do
 
  if ( abs ( 4.0 - t(1,lm1) ) <= h2tol ) go to 130
  if ( t(1,lm1) == 0.0 ) go to 120
  if ( abs ( 2.0 - abs ( t(1,lm1) ) ) < tljump ) go to 280
  if (l==3) go to 30
  h2conv = .false.
  if ( abs ( ( t(1,lm1) - t(1,l-2) ) / t(1,lm1) ) <= aittol ) go to 160
 
  if ( .not. reglar .and. l == 4 ) go to 30
 
120 continue
 
  if ( errer <= ergoal ) go to 310
  go to 380

130 continue

  if ( .not. h2conv ) then
    aitken = .false.
    h2conv = .true.
  end if

140 continue

  fextrp = 4.0

150 continue

  it = it+1
  vint = step * t(l,it)
  errer = abs ( step / (fextrp-1.0) * t(it-1,l))
  if ( errer <= ergoal ) go to 340
  if ( it == lm1 ) go to 270
  if ( t(it,lm1) == 0.0 ) go to 150
  if ( t(it,lm1) <= fextrp ) go to 270

  if ( abs ( t(it,lm1) / 4.0 - fextrp ) / fextrp < aittol ) then
    fextrp = fextrp*4.0
  end if

  go to 150
 
160 continue

  if ( t(1,lm1) < aitlow ) then
    go to 380
  end if
 
  if ( .not. aitken ) then
    h2conv = .false.
    aitken = .true.
  end if
 
170 continue

  fextrp = t(l-2,lm1)
  if ( fextrp > 4.5 ) go to 140
  if ( fextrp < aitlow ) go to 380

  if ( abs ( fextrp - t(l-3,lm1) ) / t(1,lm1) > h2tol ) then
    go to 380
  end if

  sing = fextrp
  fextm1 = fextrp - 1.0

  ait(1) = 0.0
  do i = 2, l
    ait(i) = t(i,1) + (t(i,1)-t(i-1,1)) / fextm1
    r(i) = t(1,i-1)
    dif(i) = ait(i) - ait(i-1)
  end do

  it = 2

190 continue

  vint = step*ait(l)

200 continue

  errer = errer / fextm1
 
  if ( errer <= ergoal ) then
    ind = max ( ind, 2 )
    go to 340
  end if
 
210 continue

  it = it+1
  if ( it == lm1 ) go to 270

  if ( it <= 3 ) then
    h2next = 4.0
    singnx = 2.0 * sing
  end if

  if ( h2next < singnx ) go to 230
  fextrp = singnx
  singnx = 2.0 * singnx
  go to 240

230 continue

  fextrp = h2next
  h2next = 4.0 * h2next

240 continue
 
  do i = it, lm1
    if ( abs ( dif(i+1) ) > tabtlm ) then
      r(i+1) = dif(i) / dif(i+1)
    else
      r(i+1) = 0.0
    end if
  end do
 
  h2tfex = -h2tol*fextrp
  if ( r(l) - fextrp < h2tfex ) go to 270
  if ( r(l-1) - fextrp < h2tfex ) go to 270
  errer = astep * abs ( dif(l) )
  fextm1 = fextrp - 1.0
  do i = it, l
    ait(i) = ait(i)+dif(i) / fextm1
    dif(i) = ait(i)-ait(i-1)
  end do
 
  go to 190
 
270 continue

  fextrp = max(prever/errer,aitlow)
  prever = errer
  if (l<5) go to 40
  if (l-it>2.and.istage<mxstge) go to 370
  if (errer/fextrp**(maxtbl-l)<ergoal) go to 40
  go to 370
 
280 continue

  if ( errer > ergoal ) go to 370
  diff = abs ( t(1,l) ) * 2.0 * fn
  go to 340
 
290 continue

  slope = (fend-fbeg) * 2.0
  fbeg2 = fbeg * 2.0
 
  do i = 1, 4
    diff = abs ( func ( beg + rn(i) * step ) - fbeg2 - rn(i) * slope )
    if ( diff > tabtlm ) go to 330
  end do
 
  go to 340
 
310 continue

  slope = (fend-fbeg)*2.0
  fbeg2 = fbeg*2.0
  i = 1
 
320 continue

  diff = abs ( func(beg+rn(i)*step) - fbeg2 - rn(i) * slope )
 
330 continue

  errer = max ( errer, astep * diff )
  if (errer > ergoal) go to 380
  i = i+1
  if ( i <= 4 ) go to 320
  ind = 3
 
340 continue

  cadref = cadref+vint
  error = error+errer
 
350 continue

  if (right) go to 360
  istage = istage-1
  if (istage==0) return
  reglar = reglsv(istage)
  beg = begin(istage)
  end = finis(istage)
  curest = curest-est(istage+1)+vint
  iend = ibeg-1
  fend = ts(iend)
  ibeg = ibegs(istage)
  go to 400
 
360 continue

  curest = curest+vint
  stage = stage*2.0
  iend = ibeg
  ibeg = ibegs(istage)
  end = beg
  beg = begin(istage)
  fend = fbeg
  fbeg = ts(ibeg)
  go to 10
 
370 continue

  reglar = .true.
 
380 continue
 
  if ( istage == mxstge ) then
    ind = 5
    cadref = curest+vint
    return
  end if
 
390 continue

  if (right) go to 410
  reglsv(istage+1) = reglar
  begin(istage) = beg
  ibegs(istage) = ibeg
  stage = stage/2.0

400 continue

  right = .true.
  beg = (beg+end)/2.0
  ibeg = (ibeg+iend)/2
  ts(ibeg) = ts(ibeg) / 2.0
  fbeg = ts(ibeg)
  go to 20

410 continue

  nnleft = ibeg-ibegs(istage)
  if (end+nnleft>=maxts) go to 440
  iii = ibegs(istage)
  ii = iend
  do i = iii, ibeg
    ii = ii+1
    ts(ii) = ts(i)
  end do
 
  do i = ibeg, ii
    ts(iii) = ts(i)
    iii = iii+1
  end do
 
  iend = iend+1
  ibeg = iend-nnleft
  fend = fbeg
  fbeg = ts(ibeg)
  finis(istage) = end
  end = beg
  beg = begin(istage)
  begin(istage) = end
  reglsv(istage) = reglar
  istage = istage+1
  reglar = reglsv(istage)
  est(istage) = vint
  curest = curest+est(istage)
  go to 10

440 continue

  ind = 4

460 continue

  cadref = curest+vint

  return
  end function cadref

