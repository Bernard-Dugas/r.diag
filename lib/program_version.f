      SUBROUTINE PROGRAM_VERSION ( mode )

      IMPLICIT      none

      CHARACTER*(*) mode

***    Auteur: B.Dugas

***    Objet: Imprimer de l'information sur la version courante.

      CHARACTER *80 RDIAG,VERSION,AMODE*3
      CHARACTER     DAT*20,REV*6

      EXTERNAL      RMNLIB_version
*---------------------------------------------------------------------

      AMODE = mode
      CALL LOW2UP( AMODE,AMODE )

      REV   = 'x.y.z'
      DAT   = 'XXX YY ZZZZ'
      RDIAG = 'Version ' // trim( REV ) // ' ' // trim( DAT )

      CALL RMNLIB_version( VERSION, .false.)

      IF (AMODE.EQ.'REV')                                      THEN

          write(6,'(A)') REV
          stop

      ELSE IF (AMODE.EQ.'DAT')                                 THEN

          write(6,'(A)') DAT
          stop

      ELSE IF (AMODE.EQ.'ALL' .OR.AMODE.EQ.'OUI')              THEN

          write(6,6000) trim( RDIAG ),trim( VERSION )

      END IF

      RETURN

*---------------------------------------------------------------------
 6000 format(/' The current R.DIAG is : ',A/
     +        ' And it is linked with ',A/)

      END
