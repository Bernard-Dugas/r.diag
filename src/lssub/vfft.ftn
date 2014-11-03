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
C     $Log: vfft.ftn,v $
C     Revision 3.12  2014/10/16 12:00:46  dugas
C     Modifications commandees par la compilation avec GFORTRAN et IFORT.
C
C     Revision 3.11  2014/09/25 18:42:04  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.10  2011/11/22 19:16:56  dugas
C     Corriger deux declarations du parametre MAXFAC.
C
C     Revision 3.9  2010/03/27 21:15:50  dugas
C     - Toutes les references a Temperton ont ete
C       remplacees par des references a ECMWF.
C     - De meme, le TYPE='TEMP' devient TYPE='ECMW'.
C     - Appel a FFT772 remplace par un appel a FFT_RD.
C
C     Revision 3.8  2008/04/30 19:47:55  dugas
C     Corriger l'usage des macros pour r.gppf (passe 3).
C
C     Revision 3.7  2008/04/28 21:38:53  dugas
C     Corriger l'usage des macros pour r.gppf (passe 2).
C
C     Revision 3.6  2007/08/29 20:00:29  dugas
C     Utiliser ALLOCATE/DEALLOCATE plutot que HPALLOC/HPDEALLC et ajouter FFFT8_B/RSFFT4.
C
C     Revision 3.5  2006/07/04 16:29:02  dugas
C     - Faire passer lot_maximum de 4 a 256.
C     - Utiliser RSFFT3 dans RSFFT au lieu de RSFFT2.
C     - Nouveau RSFFT3 qui fait des transformees complexes de
C       longueur N/2, ou N est le nombre de points reels.
C     - ALLOCATE/DEALLOCATE des champs de travail dans RSFFT et CSFFT.
C
C     Revision 3.4  2005/07/28 17:24:12  dugas
C     Modifier le code pour enlever les messages d'avertissement de F90.
C
C     Revision 3.3  2003/09/15 16:19:46  dugas
C     Re-definir le macro lot_maximum de 1024 a 4 suite a
C      une modification equivalente dans les versions libpriv.a
C      de RPASSM8 et QPASSM8 sur nos frontaux survenues recamment.
C
C     Revision 3.2  1995/11/01 16:56:54  armnrbd
C     De-allouer les appels multi-taches.
C
C     Revision 3.1  1994/11/17  14:14:16  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:56:21  13:56:21  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.0  93/10/13  13:32:19  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.5  93/08/19  16:20:43  16:20:43  armnrbd (Bernard Dugas)
C     Modifications cosmetiques.
C     
C     Revision 1.4  93/03/05  23:22:51  armnrbd
C     Corriger les declarations des parametres d'appel de CSFFT2.
C     
C     Revision 1.3  93/02/16  14:38:43  armnrbd
C     Corriger erreur dans boucle parrallele.
C     
C     Revision 1.2  93/02/16  11:37:26  armnrbd
C     Mettre a jour les routines rsfft,rsfft2,csfft et csfft2.
C     Finaliser le driver (stand-alone) vfft2.
C     
C     Revision 1.1  92/12/21  11:11:10  armnrbd
C     Utiliser la derniere version de SFFT (niveau FACTOR5a).
C     
C     Revision 1.0  92/02/21  11:34:53  armnrbd
C     Initial revision
C     
#     if !defined (nombre_de_facteurs)
#         define   nombre_de_facteurs 19
#     endif
#     if !defined (lot_maximum)
#         define   lot_maximum 256
#     endif
#     if !defined (nombre_de_taches)
#         define   nombre_de_taches 1
#     endif
C     SUBROUTINE 'FTSETUP2' COMPUTES FACTORS OF N & TRIGONOMETRIC
C     FUNCTIONS REQUIRED BY THE VFFT2 DRIVER.

      SUBROUTINE FTSETUP2( TRIG,IFAX,N,INC, TYPE,IER )

C     * AUTEUR: B.Dugas - 6 novembre 1991.

      IMPLICIT none

C     * PARAMETRES:
C     * TRIG (OUT) - VECTEUR DE FACTEURS TRIGONOMETRIQUES DANS LES
C     *              CAS STOCKHAM ET ECMWF, ILIST DANS LE CAS PFA.
C     * IFAX (OUT) - VECTEUR DE FACTEURS DE DECOMPOSITION DE N.
C     * N    (IN)  - NOMBRE DE POINTS DANS LA TRANSFORMEE.
C     * INC  (IN)  - INCREMENT ENTRE DEUX POINTS DE DONNEES (POUR PFA).
C     * TYPE (IN)  - PEUT FORCER LE TYPE DE TRANSFORMEE.
C     *      (OUT) - LE TYPE UTILISE.
C     * IER  (OUT) - DRAPEAU D'ERREUR. SI TOUT VA BIEN IL CONTIENT 0.
C     *              S'IL Y A UN PROBLEME, IL CONTIENT 1.

C     * NOTE:  1) LES VALEURS ACCEPTEES DE TYPE EN I/O SONT 
C     *           TYPE =  " PFA"  ET
C     *           TYPE =   "STOC"  ET
C     *           TYPE =   "ECMW" .

      INTEGER    MAXFAC
      PARAMETER( MAXFAC = nombre_de_facteurs )

      CHARACTER  TYPE*4
      INTEGER    N,IER,INC, IFAX(MAXFAC+1)
      REAL       TRIG(*)

      EXTERNAL   SETPFA2,SET77,FFTini2
C----------------------------------------------------------------------
      IER   = 1

C     * FIRST CASE: ECMWF FFT's FACTORS OF N.

      IF (TYPE /= 'STOC' .AND. TYPE /= ' PFA' )
     +    CALL SETfft_RD( N,IER )

      IF (IER /= 0 .OR. TYPE /= 'ECMW')                        THEN

C         * SECOND CASE: USE (R/Q)PFA FACTOR ROUTINE. NOTE THAT
C         * IN THIS CASE, TRIG WILL CONTAIN THE ILIST FIELD RATHER
C         * THAN THE TRIG VALUES.

          IF (TYPE.NE.'STOC') 
     +        CALL SETPFA2( IFAX,TRIG,N,INC,INC,+1, IER )

          IF (IER.NE.0 .AND. TYPE.NE.' PFA')                   THEN

C             * LAST CASE: USE STOCKHAM FACTOR ROUTINE.

              CALL FFTini2( TRIG,IFAX,N/2, IER )

C             * SIGNAL THAT THE STOCKHAM TRANSFORMS WERE CONSIDERED.

              IF (IER.EQ.0) TYPE = 'STOC'

          ELSE

C             * SIGNAL THAT THE (R/Q)PFA TRANSFORMS WERE CONSIDERED.

              TYPE = ' PFA'

          ENDIF

      ELSE

C         * SIGNAL THAT THE FFT_M8 TRANSFORMS WERE CONSIDERED.

          TYPE = 'ECMW'

      END IF

      RETURN
C-------------------------------------------------------------------

      END 

C     SUBROUTINE VFFT2 IS A MULTIPLE FAST REAL PERIODIC TRANSFORM
C     DRIVER. REAL TRANSFORMS OF LENGTH N PERFORMED  BY REMOVING
C     REDUNDANT OPERATIONS  FROM  COMPLEX  TRANSFORM OF LENGTH N.
C     ECMWF FFT_RD AND (R/Q)PFA AS WELL AS THE STOCKHAM FFT'S
C     ARE SUPPORTED, EITHER EXPLICITELY (THROUGH 'IER') OR 
C     IMPLICITELY (AS DETERMINED BY THE FACTOR DECOMPOSITIONS).

      SUBROUTINE VFFT2( A, INC,JUMP,N,LOT,ISIGN, TYPE,IER ) 

C     * AUTEUR: B.Dugas - 6 novembre 1991.

      IMPLICIT    none

C     * PARAMETERS...

C     * A IS THE ARRAY CONTAINING INPUT/OUTPUT DATA.                 (IN/OUT)
C     * INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'.              (IN)
C     *            (E.G. INC=1 FOR CONSECUTIVELY STORED DATA) 
C     * JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR. (IN)
C     *            (SHOULD BE AT LEAST (N+2))
C     * N IS THE LENGTH OF THE DATA VECTORS                          (IN)
C     * LOT IS THE NUMBER OF DATA VECTORS.                           (IN)
C     *            (SMALLER THAN MAXLOT)
C     * ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT.         (IN)
C     *       = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL.
C     * TYPE HOLDS ONE OF THE ACCEPTED TRANSFORM TYPE INDICATORS.    (IN)
C     *       THE ACCEPTED I/O TRANSFORM VALUES FOR TYPE ARE:
C     *          TYPE =  "    "  AND
C     *               =  "ECMW"  AND
C     *               =  " PFA"  AND
C     *               =  "STOC"
C     *       BLANK INPUT MEANS TRY FOR ECMW/PFA/STOC IN THAT ORDER.
C     * IER IS THE OUTPUT ERROR CONDITION.                           (OUT)
C     *       THE ERROR CONDITIONS OUTPUT ARE:
C     *           IER = 0, ALL OK
C     *               = 1, UNABLE TO FACTOR N
C     *               = 2, fft_RD unable to allocate work space
C     *               = 3, ILLEGAL FACTOR
C     *               = 4, PFA FACTORS ARE NOT MUTUALLY PRIME
C     *               = 5, STOCKHAM DOES NOT SUPPORT INC.NE.1
C     *               = 6, LOT IS TOO LARGE ("PFA" or "STOC")

C     * NOTE:
C     *  TRIG (A LIST OF TRIG FUNCTION VALUES) AND IFAX (A LIST OF
C     *  FACTORS OF N) ARE NO LONGER READ-IN. THE OLD WORK AREA IS 
C     *  NOW DYNAMICALLY REQUESTED.

C     * ORDERING OF COEFFICIENTS: 
C     *     A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2) 
C     *     WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED 

C     * ORDERING OF DATA: 
C     *     X(0),X(1),X(2),...,X(N-1), 0 , 0 ; (N+2) LOCATIONS REQUIRED 

C     * VECTORIZATION IS ACHIEVED BY DOING THE TEANSFORMS IN PARALLEL.

C     * N MUST BE COMPOSED OF FACTORS RECOGNIZED BY SET77 AND AS SUCH
C     * DOES NOT HAVE TO BE EVEN.

C     * DEFINITION OF TRANSFORMS: 
C     * ------------------------- 

C     * ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N)) 
C     *     WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K) 
C     * 
C     * ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N)) 
C     *           B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N)) 

      CHARACTER*4 TYPE
      INTEGER     INC,JUMP,N,LOT,ISIGN,IER,IERR(0:nombre_de_taches)
      INTEGER     INCLOT,NPLOT,DEBUT,FIN,M, IA, I,L
      REAL        A(N)

      INTEGER     MAXLOT
      PARAMETER ( MAXLOT = lot_maximum )

      INTEGER     N0,INC0, IFAX(20)
      SAVE        N0,INC0, IFAX

      real(8),    save, pointer :: TRIG8(:)
      real,       allocatable   :: W(:)

      real        trig
      pointer    (itrig,trig(1))
      
      CHARACTER   FTSCAL*8,ACTIVE*4
      SAVE        FTSCAL  ,ACTIVE

      INTEGER     NP
      SAVE        NP

      INTEGER     mpserv,Bidon
      EXTERNAL    mpserv

      LOGICAL               INFO
      COMMON      /ZZVERBO/ INFO

      EXTERNAL    FTSETUP2,RSFFT,FFT_RD,RPFA2,QPFA2

      DATA        NP   / 0 /
      DATA        ACTIVE / 'NONE' /

C------------------------------------------------------------------- 
C     * CHECK THAT FFT SETUP HAS BEEN DONE.

      IF (     FTSCAL.NE.'FFTINIOK' 
     +    .OR.   N   .NE.   N0    .OR.    NP .EQ.  0
     +    .OR. (TYPE .NE. ACTIVE  .AND. TYPE .NE. ' ')
     +    .OR. (INC  .NE.  INC0   .AND.' PFA'.EQ. ACTIVE) )    THEN

C         * RELEASE OLD TRIG ARRAY.

          IF (FTSCAL.EQ.'FFTINIOK' .and. 
     +        associated( TRIG8 ) ) deallocate( TRIG8 )

C         * SAVE ACTIVE TRANSFORM TYPE.

          ACTIVE = TYPE

C         * ASSIGN NEW TRIG ARRAY. DETERMINE FACTORS/TRIG CONSTANTS.

          allocate( TRIG8(3*N) ) ; itrig = loc( trig8(1) )
          CALL FTSETUP2( TRIG,IFAX,N,INC, ACTIVE,IER )

C         * QUIT IF UNSUCCESSFULL.

          IF (IER.EQ.1)                                        THEN
              IF (INFO) WRITE(6,6001) N 
              deallocate( TRIG8 )
              RETURN
          END IF

          FTSCAL = 'FFTINIOK'
          N0     =  N
          INC0   =  INC

C         * CHECK THE NUMBER OF PARALLEL THREADS AVAILABLE.

          NP     = nombre_de_taches
          NP     = mpserv('THREADS',LOT)

      ELSE
          itrig = loc( trig8(1) )
      END IF

      IF (LOT/NP.GT.MAXLOT .and. ACTIVE /= 'ECMW')             THEN

C         * LOT IS TOO LARGE.

          IER = 6

          IF (INFO) WRITE(6,6002) LOT,MAXLOT
          FTSCAL = 'NOTFTINI'
          deallocate( TRIG8 )

          RETURN

      ELSE IF (ACTIVE.EQ.'STOC'.AND. INC.NE.1)                 THEN

C         * INC.NE.1 IS NOT SUPPORTED WITH THE STOCKHAM CASE.

          IER = 5

          IF (INFO) WRITE(6,6005) N
          FTSCAL = 'NOTFTINI'
          deallocate( TRIG8 )

          RETURN

      END IF

      INCLOT = (LOT-1)/NP+1

C     * ALLOCATE WORK ARRAY FOR (R/Q)PFA.

      IF (ACTIVE.EQ.' PFA') allocate( W(INC*JUMP*LOT) )
      IF (ACTIVE.EQ.' PFA') W(1) = 0.0

C     * SEPARATE INTO PARRALLEL TASKS.

CC$doacross local( I,L,M,DEBUT,FIN,NPLOT,IA )

      DO 100 M=0,NP-1

          DEBUT  =            M  *INCLOT
          FIN    = MIN( LOT,DEBUT+INCLOT )

          NPLOT  = FIN-DEBUT

          IA     = DEBUT*JUMP + 1

C         *CALCULATE THE TRANSFORMS.

          IERR(M) = 1

          IF (ACTIVE.EQ.'ECMW')                                THEN

              CALL FFT_RD( A(IA), INC,JUMP,N,NPLOT,ISIGN, IERR(M) )

          ELSE IF (ACTIVE.EQ.' PFA')                           THEN

              IF (ISIGN.EQ.-1)                                 THEN
                  CALL RPFA2( A(IA),W(IA), 
     +                             IFAX,TRIG, INC,JUMP,INC,JUMP, 
     +                             N,NPLOT,+1, IERR(M) )
              ELSE IF (ISIGN.EQ.+1)                            THEN
                  CALL QPFA2( A(IA),W(IA),
     +                             IFAX,TRIG, INC,JUMP,INC,JUMP, 
     +                             N,NPLOT,+1, IERR(M) )
              END IF

              
              IF (IERR(M).EQ.2) IERR(M) = 3
              IF (IERR(M).EQ.1) IERR(M) = 4

              DEBUT = IA
              DO  L=1,NPLOT
                  DO  I=DEBUT,DEBUT+(N+1)*INC,INC
                      A(I) = W(I)
                  END DO
                  IF (ISIGN.EQ.+1)                             THEN
                      A(DEBUT+ N   *INC) = 0.
                      A(DEBUT+(N+1)*INC) = 0.
                  END IF
                  DEBUT = DEBUT+JUMP
              END DO

          ELSE IF (ACTIVE.EQ.'STOC')                           THEN

              IERR(M) = 0
              CALL RSFFT(  A(IA),  TRIG8,IFAX,
     +                             JUMP,N,NPLOT, ISIGN )

          END IF
#         if !defined (sgi)
          IF (IERR(M).NE.0) GOTO 200
#         endif
  100 CONTINUE
  200 CONTINUE

C     * SWITCH OFF PARRALLEL TASKS.

CC$    Bidon = mpserv('BLOCK',Bidon)

C     * SAVE ERROR EXIT CODE FROM MULTI-PROCESS LOOP.

      DO  M=0,NP-1
          IF (IERR(M).NE.0)                                    THEN
              IER = IERR(M) 
              GOTO 300
          END IF
      END DO

C     * DE-ALLOCATE (R/Q)PFA WORK ARRAY.

  300 IF (ACTIVE.EQ.' PFA') deallocate( W )

      RETURN 
C-------------------------------------------------------------------
 6001 FORMAT(' N =',I4,' - CONTAINS ILLEGAL FACTORS.')
 6002 FORMAT(' LOT=',I5,' - TOO LARGE IN VFFT. MAXLOT=',I4)
 6005 FORMAT(' N =',I4,' - REQUESTED STOCKHAM TRANSFORM, WHILE INC IS',
     +       ' NOT EQUAL TO ONE.')

      END 

      subroutine ffft8_b( a, inc, jump, lot, isign )

      implicit none

!     arguments

      integer inc, jump, lot, isign
      real*8  a(*)

!     STKFFT8 common variables

      INTEGER     MAXFAC
      PARAMETER ( MAXFAC = nombre_de_facteurs )

      CHARACTER   STRING*8
      INTEGER     N, IFAC
      REAL *8,    POINTER :: TRIG(:,:)

      COMMON    / FTSCAL  / STRING
      COMMON    / STKFFT8 / TRIG, IFAC(MAXFAC+1), N

!     local variables

      INTEGER  I,J,IJ,NP2
      REAL*8,  DIMENSION (:,:), ALLOCATABLE :: B

      EXTERNAL RSFFT4
!---------------------------------------------------------------------

      if (inc.eq.1) then

         call rsfft4( a, JUMP,LOT, ISIGN )

      else

         NP2 = 2*N+2

         allocate( B(NP2,LOT) )

         do J=1,LOT
            IJ = (J-1)*JUMP+1
            do I=1,NP2
               B(I,J) = a(IJ)
               IJ = IJ+inc
            enddo
         enddo

         call rsfft4( B, NP2,LOT, ISIGN )

         do J=1,LOT
            IJ = (J-1)*JUMP+1
            do I=1,NP2
               a(IJ) = B(I,J)
               IJ = IJ+inc
            enddo
         enddo

         deallocate( B )

      endif

      return

!---------------------------------------------------------------------
      end

      SUBROUTINE  RSFFT( A, TRIG,IFAX,JUMP,N,LOT, ISIGN )

      IMPLICIT none

      INTEGER  N,LOT,ISIGN,JUMP,IFAX(*)
      REAL     A(JUMP,LOT)
      REAL*8   TRIG(*)

      INTEGER  R1,I1,R2,I2
      REAL*8,  DIMENSION (:), ALLOCATABLE :: CCW

      LOGICAL           INFO
      COMMON  /ZZVERBO/ INFO

      EXTERNAL RSFFT3

C---------------------------------------------------------------------
      if (mod( n,2 ).ne.0)                                     THEN
          IF (INFO) print *,' RSFFT: Illegal uneven N = ',n
          RETURN
      end if

      ALLOCATE ( CCW((N+2)*LOT*2) )

      R1 = 1
      I1 = R1+(N+2)*LOT/2
      R2 = I1+(N+2)*LOT/2
      I2 = R2+(N+2)*LOT/2

      CALL rsfft3( a, CCW(R1),CCW(I1),CCW(R2),CCW(I2),
     +             trig,ifax, JUMP,N/2, LOT, ISIGN )

      DEALLOCATE ( CCW )

      RETURN
C---------------------------------------------------------------------

      END
      SUBROUTINE rsfft2( a, CCr1,CCi1,CCr2,CCi2,trigs,ifax, 
     +                      JUMP,N,LOT,LOTH, ISIGN )

c     * Auteur: B.Dugas, RPN.

c     * Historique des modifications...
c     *  ...  Fonctionnalite INC,JUMP,LOT le 6 novembre 1991.
c     *  ...  Version originale le 03 juin 1991.

c     * Version (reelle a/de complexe) de la transformee FFT complexe
c     * contenue dans SFFT. On complexifie les donnees et deux trans-
c     * formees reelles sont executees a la foi.  Un nombre impair de
c     * transformees est tout de meme acceptable.

c     * Parametres...
C     * A est le champs contenant les donnees d'entree/sortie.       (IN/OUT)
C     * CCs SONT les champs de travails de dimension LOTH*N.         (IN)
C     * JUMP est l'increment entre le debut de deux vecteurs.        (IN)
C     *            (doit etre au moins "N+2")
C     * N est la longueur de vecteurs de donnees.                    (IN)
C     * LOT est le nombre de vecteurs de donnees                     (IN)
C     *            (doit etre inferieur a MAXLOT)
C     * ISIGN = +1 pour transformee spectral a grille.               (IN)
C     *       = -1 pour transformee grille a spectral.

      IMPLICIT    none

      INTEGER     MAXFAC
      PARAMETER ( MAXFAC = nombre_de_facteurs )

      logical     ODD
      integer     JUMP,LOT
      real        a(JUMP,LOT)
      integer     II,L0,L1,L2,    IFAX(MAXFAC+1),
     +            M,M2,N,ISIGN,   ILH,LOTH,ILOTH
      real*8      CCr1(LOTH,N),   CCr2(LOTH,0:N-1),  
     +            CCi1(LOTH,N),   CCi2(LOTH,0:N-1),   
     +            trigs(3*N),     CCrM,CCrL,CCiM,CCiL,
     +            ON,DON,         a11,a21,a12,a22

      EXTERNAL    sfft

c--------------------------------------------------------------------------
      ILOTH =  LOT /2
      ILH   = (N+1)/2

                         odd = .false.
      if (ILOTH.NE.LOTH) odd = .true.

      ON    = 1.0/DBLE( 2*N )
      DON   = 2.0*ON

      if (ISIGN.eq.-1)                                         then
      
c        * Forward transform (real to spectral).
c        * Transfer grid data from a to CC1.

         IF (N.GT.ILOTH)                                       THEN

            L1 = -1

            DO L0=1,ILOTH

               L1    = L1+2
               DO II =1,N

                  CCr1(L0,II) = a(II,L1  )
                  CCi1(L0,II) = a(II,L1+1)

               END DO

            END DO

         ELSE

            DO II =1,N

               L1 = -1

               DO L0=1,ILOTH

                  L1    = L1+2

                  CCr1(L0,II) = a(II,L1  )
                  CCi1(L0,II) = a(II,L1+1)

               END DO

            END DO

         END IF

         IF (odd)                                              THEN
            DO II =1,N

               CCr1(LOTH,II) = a(II,LOT)
               CCi1(LOTH,II) = 0.0

            END DO
         END IF

c        * Call actual complex transforms. Results in CC2.

         call sfft( -1, N,LOTH,ifax,CCr1,CCi1,CCr2,CCi2,trigs )

         L1   = -1
#        if defined (NEC)
*vdir nodep
#        endif
         DO L0=1,ILOTH

            L1 =  L1+2

c           * Set b wave number 0 components.

            a(1,L1)   = CCr2(L0,0)*DON
            a(1,L1+1) = CCi2(L0,0)*DON
            a(2,L1)   = 0.0
            a(2,L1+1) = 0.0

         END DO

         IF (odd)                                              THEN
            a(1,LOT)   = CCr2(LOTH,0)*DON
            a(2,LOT)   = 0.0
         END IF

         IF (ILH.GT.ILOTH)                                     THEN

            L1 =  -1
            DO L0=1,ILOTH

               L1 =  L1+2

               DO M=1,ILH

                  L2 = N-M
                  M2 = M+M+1

                  CCrM = CCr2(L0,M)
                  CCrL = CCr2(L0,L2)
                  CCiM = CCi2(L0,M)
                  CCiL = CCi2(L0,L2)

c                 * Transfer wave numbers (1 --> ILH) to b.

                  a(M2  ,L1  ) = ( +CCrM + CCrL ) * ON
                  a(M2+1,L1  ) = ( +CCiM - CCiL ) * ON
                  a(M2  ,L1+1) = ( +CCiM + CCiL ) * ON
                  a(M2+1,L1+1) = ( -CCrM + CCrL ) * ON

               END DO

            END DO

         ELSE

            DO M=1,ILH

               L1 =  -1
               L2 = N-M
               M2 = M+M+1

               DO L0=1,ILOTH

                  L1 =  L1+2

                  CCrM = CCr2(L0,M)
                  CCrL = CCr2(L0,L2)
                  CCiM = CCi2(L0,M)
                  CCiL = CCi2(L0,L2)

c                 * Transfer wave numbers (1 --> ILH) to b.

                  a(M2  ,L1  ) = ( +CCrM + CCrL ) * ON
                  a(M2+1,L1  ) = ( +CCiM - CCiL ) * ON
                  a(M2  ,L1+1) = ( +CCiM + CCiL ) * ON
                  a(M2+1,L1+1) = ( -CCrM + CCrL ) * ON

               END DO

            END DO

         END IF

         IF (odd)                                              THEN
            DO M=1,ILH

               L2 = N-M
               M2 = M+M+1

               CCrM = CCr2(LOTH,M)
               CCrL = CCr2(LOTH,L2)
               CCiM = CCi2(LOTH,M)
               CCiL = CCi2(LOTH,L2)

c                 * Transfer wave numbers (1 --> ILH) to b.

               a(M2  ,LOT) = ( +CCrM + CCrL ) * ON
               a(M2+1,LOT) = ( +CCiM - CCiL ) * ON

            END DO
         END IF

      else

c        * Backward transform (spectral to real).
c        * Transfer spectral data from b to CC2.

         L1 = -1
         DO L0=1,ILOTH

            L1 =  L1+2

c           * Set wave number 0 components.

            a11        = a(1,L1)
            a12        = a(1,L1+1)

            CCr2(L0,0) = a11
            CCi2(L0,0) = a12

         END DO

         IF (odd)                                              THEN
            CCr2(LOTH,0) = a(1,LOT)
            CCi2(LOTH,0) = 0.0
         END IF

         IF (ILH.GT.ILOTH)                                     THEN

            L1 =  -1
            DO L0=1,ILOTH

               L1 =  L1+2

               DO M=1,ILH

                  L2 = N-M
                  M2 = M+M+1

                  a11 = a(M2  ,L1  )
                  a21 = a(M2+1,L1  )
                  a12 = a(M2  ,L1+1)
                  a22 = a(M2+1,L1+1)

c                 * Transfer wave numbers (1 --> ILH) to CC2.

                  CCr2(L0,M ) = +a11 - a22
                  CCi2(L0,M ) = +a21 + a12
                  CCr2(L0,L2) = +a11 + a22
                  CCi2(L0,L2) = -a21 + a12

               END DO

            END DO

         ELSE

            DO M=1,ILH

               L2 = N-M
               M2 = M+M+1

               L1 =  -1
               DO L0=1,ILOTH

                  L1 =  L1+2

                  a11 = a(M2  ,L1  )
                  a21 = a(M2+1,L1  )
                  a12 = a(M2  ,L1+1)
                  a22 = a(M2+1,L1+1)

c                 * Transfer wave numbers (1 --> ILH) to CC2.

                  CCr2(L0,M ) = +a11 - a22
                  CCi2(L0,M ) = +a21 + a12
                  CCr2(L0,L2) = +a11 + a22
                  CCi2(L0,L2) = -a21 + a12

               END DO

            END DO

         END IF

         IF (odd)                                              THEN
            DO M=1,ILH

               L2 = N-M
               M2 = M+M+1

               a11 = a(M2  ,LOT)
               a21 = a(M2+1,LOT)

c              * Transfer wave numbers (1 --> ILH) to CC2.

               CCr2(LOTH,M ) = +a11
               CCi2(LOTH,M ) = +a21
               CCr2(LOTH,L2) = +a11
               CCi2(LOTH,L2) = -a21

            END DO
         END IF

         call sfft( +1, N,LOTH,ifax,CCr2,CCi2,CCr1,CCi1,trigs )

c        * Transfer grid data from CC1 to a.

         IF (N.GT.ILOTH)                                       THEN

            L1 = -1
            DO L0=1,ILOTH

               L1 = L1+2
#if defined (NEC)
*vdir nodep
#endif
               DO II =1,N

                  a(II,L1  ) = CCr1(L0,II)
                  a(II,L1+1) = CCi1(L0,II)

               END DO

            END DO

         ELSE

            DO II =1,N

               L1 = -1
               DO L0=1,ILOTH

                  L1         = L1+2
#if defined (NEC)
*vdir nodep
#endif
                  a(II,L1  ) = CCr1(L0,II)
                  a(II,L1+1) = CCi1(L0,II)

               END DO

            END DO

         END IF

         IF (odd)                                              THEN
            DO II =1,N

               a(II,LOT) = CCr1(LOTH,II)

            END DO
         END IF

      end if

      return
c-----------------------------------------------------------------------

      end
      SUBROUTINE rsfft3( a, CCr1,CCi1,CCr2,CCi2,trigs,ifax, 
     +                      JUMP,N,LOT, ISIGN )
c
c     * Auteur: B.Dugas, RPN.
c
c     * Version (reelle a/de complexe) de la transformee FFT complexe
c     * contenue dans SFFT. On complexifie les 2*n donnees reelles et
c     * une transforme complexe de longueur n est executee.
c
*     ****************************************************************
*     * L'algorithme utilise dans cette routine est presente
*     * dans "Numerical recipes, The Art of Scientific Computing",
*     * Cambridge University Press, 1986. Voir pages 390-400.
*     ****************************************************************
C
c     * Parametres...
C     * A est le champs contenant les donnees d'entree/sortie        (IN/OUT)
C     * CCs SONT les champs de travails de dimension LOT*(N+1)       (IN)
C     * JUMP est l'increment entre le debut de deux vecteurs         (IN)
C     *            (doit etre .GE. 2*N+2)
C     * N est la longueur de vecteurs de donnees complexes           (IN)
C     * LOT est le nombre de vecteurs de donnees                     (IN)
C     * ISIGN = +1 pour transformee spectral a grille                (IN)
C     *       = -1 pour transformee grille a spectral 

      IMPLICIT    none

      INTEGER           MAXFAC
      PARAMETER       ( MAXFAC = nombre_de_facteurs )
      integer     IFAX( MAXFAC+1 )

      integer     JUMP,LOT,       N,ISIGN
      real*8      CCr1(LOT,N),    CCr2(LOT,0:N),  
     +            CCi1(LOT,N),    CCi2(LOT,0:N),   
     +            trigs(3*N)
      real        a(JUMP,LOT)

C     * Variables locales.

      logical     ODD
      integer     II,L,L0,L1,L2, 
     +            ILOTH,M,M2,     NS2,NP2
      real*8      CCrM,CCrL,      CCiM,CCiL,
     +            ON,DON,         a11,a12,
     +            Wr,Wi,          Wtemp,Wpi,Wpr,pis,
     +            H1r,H1i,        H2r,H2i, C1,C2

      LOGICAL              INFO
      COMMON     /ZZVERBO/ INFO

      EXTERNAL    sfft

c--------------------------------------------------------------------------
      NS2 = N/2
      NP2 = 2*N+2

      IF (JUMP.LT.NP2)                                         THEN
          IF (INFO) PRINT *,' RSFFT3: JUMP < NP2; JUMP,NP2= ',JUMP,NP2
          RETURN
      END IF
 
      pis = -2.0*ASIN( 1.0D0 )*ISIGN
      Wpr = -2.*SIN( 0.5*pis/N )**2
      Wpi =     SIN(     pis/N )
 
      ON  = 1.0/DBLE( N )
      DON = 2.0*ON
 
      if (ISIGN.eq.-1)                                         then

c        * Forward transform (real to spectral).
c        * Transfer grid data from a to CC1.

         DO L0=1,LOT
            DO II =1,N
               L1  = 2*II-1
               CCr1(L0,II) = a(L1  ,L0)
               CCi1(L0,II) = a(L1+1,L0)
            END DO
         END DO

c        * Call actual complex transforms. Results in CC2.

         call sfft( -1, N, LOT,ifax,CCr1,CCi1,CCr2,CCi2,trigs )

         C1  = 0.25
         DO L0=1,LOT

c           * Set b wave number 0 and N components.

            CCrm      = CCr2(L0,0)
            CCim      = CCi2(L0,0)

            a(1    ,L0) = ( CCrm + CCim ) * DON * C1
            a(2    ,L0) =   0.0
            a(NP2-1,L0) = ( CCrm - CCim ) * DON * C1
            a(NP2  ,L0) =   0.0

         END DO

         C1 =  C1   * ON
         C2 = -0.25 * ON

         DO L0=1,LOT

            Wr  = 1.+Wpr
            Wi  = Wpi

            DO M=1,NS2

               L  = N-M
               M2 = M+M+1
               L2 = NP2-M2

               CCrM = CCr2(L0,M)
               CCrL = CCr2(L0,L)
               CCiM = CCi2(L0,M)
               CCiL = CCi2(L0,L)

c              * Transfer wave numbers (1 --> N-1) to b.

               H1r        = ( CCrM + CCrL ) * C1
               H1i        = ( CCiM - CCiL ) * C1
               H2r        = ( CCiM + CCiL ) * C2
               H2i        = ( CCrM - CCrL ) * C2

               a(M2  ,L0) = H1r - ( Wr*H2r - Wi*H2i )
               a(M2+1,L0) = H1i + ( Wr*H2i + Wi*H2r )
               a(L2  ,L0) = H1r + ( Wr*H2r - Wi*H2i )
               a(L2+1,L0) =-H1i + ( Wr*H2i + Wi*H2r )

               Wtemp = Wr
               Wr = Wr*Wpr-Wi   *Wpi+Wr
               Wi = Wi*Wpr+Wtemp*Wpi+Wi

            END DO

         END DO

      else

c        * Backward transform (spectral to real).
c        * Transfer spectral data from b to CC2.

         DO L0=1,LOT

c           * Set wave number 0 components.

            a11        = a(1    ,L0)
            a12        = a(NP2-1,L0)
            CCr2(L0,0) = a11+a12
            CCi2(L0,0) = a11-a12

         END DO

         DO L0=1,LOT

            Wr  = Wpr+1.
            Wi  = Wpi

            DO M=1,NS2

               L  = N-M
               M2 = M+M+1
               L2 = NP2-M2

               CCrM = a(M2  ,L0)
               CCrL = a(L2  ,L0)
               CCiM = a(M2+1,L0)
               CCiL = a(L2+1,L0)

c              * Transfer wave numbers (1 --> N-1) to CC2.

               H1r        = ( CCrM + CCrL )
               H1i        = ( CCiM - CCiL )
               H2r        = ( CCiM + CCiL )
               H2i        = ( CCrM - CCrL )

               CCr2(L0,M) = H1r - ( Wr*H2r - Wi*H2i )
               CCi2(L0,M) = H1i + ( Wr*H2i + Wi*H2r )
               CCr2(L0,L) = H1r + ( Wr*H2r - Wi*H2i )
               CCi2(L0,L) =-H1i + ( Wr*H2i + Wi*H2r )

               Wtemp = Wr
               Wr = Wr*Wpr-Wi   *Wpi+Wr
               Wi = Wi*Wpr+Wtemp*Wpi+Wi

            END DO

         END DO

         call sfft( +1, N, LOT,ifax,CCr2,CCi2,CCr1,CCi1,trigs )

c        * Transfer grid data from CC1 to a.

         DO L0=1,LOT
            DO II =1,N
               L1 = 2*II-1
               a(L1  ,L0) = CCr1(L0,II)
               a(L1+1,L0) = CCi1(L0,II)
            END DO
         END DO

      end if

      return
c-----------------------------------------------------------------------

      end

      SUBROUTINE rsfft4( a, JUMP,LOT, ISIGN )

      IMPLICIT    none

      integer     JUMP,LOT,ISIGN
      real*8      a(JUMP,LOT)

!     * Auteur: B.Dugas, RPN.

!     * Version (reelle a/de complexe) de la transformee FFT complexe
!     * contenue dans SFFT2. On complexifie les 2*n donnees reelles et
!     * une transforme complexe de longueur n est executee.

!     ****************************************************************
!     * L'algorithme utilise dans cette routine est presente
!     * dans "Numerical recipes, The Art of Scientific Computing",
!     * Cambridge University Press, 1986. Voir pages 390-400.
!     ****************************************************************

!     * Parametres...
!     * A est le champs contenant les donnees d'entree/sortie        (IN/OUT)
!     * JUMP est l'increment entre le debut de deux vecteurs         (IN)
!     *            (doit etre .GE. 2*N+2)
!     * LOT est le nombre de vecteurs de donnees                     (IN)
!     * ISIGN = +1 pour transformee spectral a grille                (IN)
!     *       = -1 pour transformee grille a spectral 

!     STKFFT8 common variables

      INTEGER     MAXFAC
      PARAMETER ( MAXFAC = nombre_de_facteurs )

      CHARACTER   STRING*8
      INTEGER     N, IFAC
      REAL *8,    POINTER :: TRIG(:,:)

      COMMON    / FTSCAL  / STRING
      COMMON    / STKFFT8 / TRIG, IFAC(MAXFAC+1), N

      LOGICAL               INFO
      COMMON    / ZZVERBO / INFO

!     local variables.

      logical     ODD

      integer     II,ILOTH,       L,L0,L1,L2,  
     +            M,M2,NS2,NP2
      real*8      CCrM,CCrL,      CCiM,CCiL,  
     +            ON,DON,         a11,a12,    
     +            Wr,Wi,          Wtemp,Wpi,Wpr,pis, 
     +            H1r,H1i,        H2r,H2i, C1,C2     

      REAL*8,  DIMENSION (:,:), ALLOCATABLE :: CCr1,CCi1,CCr2,CCi2

      EXTERNAL    sfft2
!--------------------------------------------------------------------------
      NS2 = N/2
      NP2 = 2*N+2

      IF (JUMP.LT.NP2)                                         THEN
          IF (INFO) PRINT *,' RSFFT4: JUMP < NP2; JUMP,NP2= ',JUMP,NP2
          RETURN
      END IF

      pis = -2.0*ASIN( 1.0D0 )*ISIGN
      Wpr = -2.*SIN( 0.5*pis/N )**2
      Wpi =     SIN(     pis/N )

      ON  = 1.0/DBLE( N )
      DON = 2.0*ON

      ALLOCATE( CCr1(LOT,  N),CCi1(LOT,  N),
     +          CCr2(LOT,0:N),CCi2(LOT,0:N) )

      if (ISIGN.eq.-1)                                         then

!        * Forward transform (real to spectral).
!        * Transfer grid data from a to CC1.

         DO L0=1,LOT
            DO II =1,N
               L1  = 2*II-1
               CCr1(L0,II) = a(L1  ,L0)
               CCi1(L0,II) = a(L1+1,L0)
            END DO
         END DO

!        * Call actual complex transforms. Results in CC2.

         call sfft2( -1, LOT, CCr1,CCi1,CCr2,CCi2 )

         C1  = 0.25

         DO L0=1,LOT

!           * Set b wave number 0 and N components.

            CCrm      = CCr2(L0,0)
            CCim      = CCi2(L0,0)

            a(1    ,L0) = ( CCrm + CCim ) * DON * C1
            a(2    ,L0) =   0.0
            a(NP2-1,L0) = ( CCrm - CCim ) * DON * C1
            a(NP2  ,L0) =   0.0

         END DO

         C1 =  C1   * ON
         C2 = -0.25 * ON

         DO L0=1,LOT

            Wr  = 1.+Wpr
            Wi  = Wpi

            DO M=1,NS2

               L  = N-M
               M2 = M+M+1
               L2 = NP2-M2

               CCrM = CCr2(L0,M)
               CCrL = CCr2(L0,L)
               CCiM = CCi2(L0,M)
               CCiL = CCi2(L0,L)

!              * Transfer wave numbers (1 --> N-1) to b.

               H1r        = ( CCrM + CCrL ) * C1
               H1i        = ( CCiM - CCiL ) * C1
               H2r        = ( CCiM + CCiL ) * C2
               H2i        = ( CCrM - CCrL ) * C2

               a(M2  ,L0) = H1r - ( Wr*H2r - Wi*H2i )
               a(M2+1,L0) = H1i + ( Wr*H2i + Wi*H2r )
               a(L2  ,L0) = H1r + ( Wr*H2r - Wi*H2i )
               a(L2+1,L0) =-H1i + ( Wr*H2i + Wi*H2r )

               Wtemp = Wr
               Wr = Wr*Wpr-Wi   *Wpi+Wr
               Wi = Wi*Wpr+Wtemp*Wpi+Wi

            END DO

         END DO

      else

!        * Backward transform (spectral to real).
!        * Transfer spectral data from b to CC2.

         DO L0=1,LOT

!           * Set wave number 0 components.

            a11        = a(1    ,L0)
            a12        = a(NP2-1,L0)
            CCr2(L0,0) = a11+a12
            CCi2(L0,0) = a11-a12

         END DO

         DO L0=1,LOT

            Wr  = Wpr+1.
            Wi  = Wpi

            DO M=1,NS2

               L  = N-M
               M2 = M+M+1
               L2 = NP2-M2

               CCrM = a(M2  ,L0)
               CCrL = a(L2  ,L0)
               CCiM = a(M2+1,L0)
               CCiL = a(L2+1,L0)

!              * Transfer wave numbers (1 --> N-1) to CC2.

               H1r        = ( CCrM + CCrL )
               H1i        = ( CCiM - CCiL )
               H2r        = ( CCiM + CCiL )
               H2i        = ( CCrM - CCrL )

               CCr2(L0,M) = H1r - ( Wr*H2r - Wi*H2i )
               CCi2(L0,M) = H1i + ( Wr*H2i + Wi*H2r )
               CCr2(L0,L) = H1r + ( Wr*H2r - Wi*H2i )
               CCi2(L0,L) =-H1i + ( Wr*H2i + Wi*H2r )

               Wtemp = Wr
               Wr = Wr*Wpr-Wi   *Wpi+Wr
               Wi = Wi*Wpr+Wtemp*Wpi+Wi

            END DO

         END DO

         call sfft2( +1, LOT, CCr2,CCi2,CCr1,CCi1 )

!        * Transfer grid data from CC1 to a.

         DO L0=1,LOT
            DO II =1,N
               L1 = 2*II-1
               a(L1  ,L0) = CCr1(L0,II)
               a(L1+1,L0) = CCi1(L0,II)
            END DO
         END DO

      end if

      DEALLOCATE( CCr1,CCi1,CCr2,CCi2 )

      return

!-----------------------------------------------------------------------
      end

      SUBROUTINE  CSFFT( A, JUMP,N,LOT, ISIGN )

C     * AUTEUR: B.DUGAS - AUTOMNE 1992.

      IMPLICIT    none

      INTEGER     MAXFAC
      PARAMETER(  MAXFAC = nombre_de_facteurs + 1 )

      REAL        A(2,*)
      INTEGER     N,LOT,ISIGN,JUMP

      REAL        TRIG
      REAL*8      TRIG8
      REAL*8,     DIMENSION (:), ALLOCATABLE :: CC
      INTEGER     LOT2, R1,I1,R2,I2, ERRCOD, IFAX(MAXFAC)

      POINTER   ( PTR8 , TRIG8(1) ),( PTR , TRIG(1) )
      SAVE        PTR8

      CHARACTER   FTSCAL*8,ACTIVE*4
      SAVE        FTSCAL

      INTEGER     N0,SCC0
      SAVE        N0,SCC0

      LOGICAL              INFO
      COMMON     /ZZVERBO/ INFO

      EXTERNAL    HPALLOC,HPDEALLC,CSFFT2,FTSETUP2

      DATA        FTSCAL / ' ' /
      DATA        N0     /  0  /
      DATA        SCC0   /  0  /

c-------------------------------------------------------------------
      LOT2 = LOT
      IF (MOD(LOT,2).EQ.0) LOT2 = LOT + 1

      IF (FTSCAL.NE.'FFTINIOK' .OR. N      .NE. N0
     +                         .OR. N*LOT2 .GT. SCC0)          THEN

          SCC0 = N*LOT2

C         * RELEASE OLD TRIG ARRAY.

          IF (FTSCAL.EQ.'FFTINIOK') CALL HPDEALLC( PTR8,ERRCOD,0 )

C         * SET ACTIVE TRANSFORM TYPE.

          ACTIVE = 'STOC'

C         * ASSIGN NEW TRIG ARRAY. DETERMINE FACTORS/TRIG CONSTANTS.

          CALL HPALLOC( PTR8,3*N,ERRCOD,8 )
          PTR = LOC( TRIG8(1) )

          CALL FTSETUP2( TRIG,IFAX,N,1, ACTIVE,ERRCOD )

C         * QUIT IF UNSUCCESSFULL.

          IF (ERRCOD.EQ.1)                                     THEN
              IF (INFO) WRITE(6,6000) N 
              FTSCAL = ' '
              RETURN
          END IF

          N0     =  N
          FTSCAL = 'FFTINIOK'

      END IF

      ALLOCATE ( CC(SCC0*4) )

      r1 = 1
      i1 = r1+N*LOT2
      r2 = i1+N*LOT2
      i2 = r2+N*LOT2

      CALL csfft2( a, CC(r1),CC(i1),CC(r2),CC(i2),
     +               trig8,ifax, JUMP,N,LOT,LOT2,ISIGN )

      DEALLOCATE ( CC )

      RETURN
c-------------------------------------------------------------------

 6000 FORMAT(' N =',I4,' - CONTAINS ILLEGAL FACTORS.')

      END
      SUBROUTINE csfft2( a, CCr1,CCi1,CCr2,CCi2,trigs,ifax, 
     +                      JUMP,N,LOT,LOT2,ISIGN )

c     * Auteur: B.Dugas, RPN - AUTOMNE 1992.

c     * Version (complexe a complexe) de la transformee FFT complexe
c     * contenue dans SFFT.

c     * Parametres...
C     * A est le champs complexe d'entree/sorties.                   (IN/OUT)
C     * CC est le champs de travail de dimension LOT2*N*4.           (IN)
C     * JUMP est l'increment entre le debut de deux vecteurs.        (IN)
C     *            (doit etre au moins N)
C     * N est la longueur de vecteurs de donnees.                    (IN)
C     * LOT est le nombre de vecteurs de donnees                     (IN)
C     *            (doit etre inferieur a MAXLOT)
C     * LOT2 est le eal a LOT+mod(LOT+1,2).                          (IN)
C     * ISIGN = +1 pour transformee spectral a grille.               (IN)
C     *       = -1 pour transformee grille a spectral.

      IMPLICIT    none

      INTEGER     MAXFAC
      PARAMETER ( MAXFAC = nombre_de_facteurs + 1 )

      integer     II,JJ,          IFAX(MAXFAC), JUMP
      integer     N,ISIGN,        LOT,LOT2
      real*8      CCr1(LOT2,N),   CCr2(LOT2,N), ON
      real*8      CCi1(LOT2,N),   CCi2(LOT2,N), trigs(3*N)
      real        a(2,JUMP,LOT)

      EXTERNAL    sfft

c--------------------------------------------------------------------------
      if (ISIGN.eq.-1)                                         then
          ON = 1.0/FLOAT( N )
      else
          ON = 1.0
      end if

      IF (N.GT.LOT)                                            THEN

         DO JJ=1,LOT

            DO II =1,N

               CCr1(JJ,II) = a(1,II,JJ)
               CCi1(JJ,II) = a(2,II,JJ)

            END DO

         END DO

      ELSE

         DO II =1,N

            DO JJ=1,LOT

               CCr1(JJ,II) = a(1,II,JJ)
               CCi1(JJ,II) = a(2,II,JJ)

            END DO

         END DO

      END IF

      IF (LOT.NE.LOT2)                                         THEN

          DO II =1,N
             CCr1(LOT2,II) = 0.0
             CCi1(LOT2,II) = 0.0
          END DO

      END IF

c     * Call actual complex transforms. Results in CC2.

      call sfft( ISIGN, N,LOT2,ifax,CCr1,CCi1,CCr2,CCi2,trigs )

      IF (N.GT.LOT)                                            THEN

         DO JJ=1,LOT

            DO II=1,N

               a(1,II,JJ) = CCr2(JJ,II) * ON
               a(2,II,JJ) = CCi2(JJ,II) * ON

            END DO

         END DO

      ELSE


         DO II=1,N

            DO JJ=1,LOT

               a(1,II,JJ) = CCr2(JJ,II) * ON
               a(2,II,JJ) = CCi2(JJ,II) * ON

            END DO

         END DO

      END IF

      return
c-----------------------------------------------------------------------

      end
