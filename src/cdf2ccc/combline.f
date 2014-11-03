!
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
!
      subroutine combline5( VALUE,gg,INDICE,K,NI,NJ,NK,SCALE,OFFSET,
     .                      MULT,ADD,FILL_OK,FILL_CDF,INVJJ,LASLON,
     .                      FILL_ALL,FILL_CDF_NAN )

      implicit none
      
      include 'cdf2ccc.h'
      
      integer indice,k,ni,nj,nk,laslon

      logical fill_ok,invjj,fill_all,fill_cdf_nan
      real*8  scale,offset,mult,add,fill_cdf
      real*8  value(ni*nj*nk)
      real*8  gg(ni*nj)

******
*
* AUTEUR Guy Bergeron   juin 2003
*
*     Effectue la combinaison lineaire suivante :
*
*                gg=(scale*val + offset)*mult + add
*
*         scale   : "scale_factor" utilises dans l'algorythme de compression
*         offset  : "add_offset" utilises dans l'algorythme de compression
*
*         mult    : facteur multiplicatif changement d'unites
*         add     : facteur additif changement d'unites
*
*         fill_ok : variable logique pour le remplacement de fill_cdf
*         fill_cdf: Valeur de remplacement dans le fichier netcdf
*
*         invjj   : inverse l'ordre de l'indice "j" dans le vecteur de sortie.
*
* REVISIONS
*
* B.Dugas fev '14 :
* - Tenir compte du cas ou FILL_CDF est un NaN
* B.Dugas juin '13 :
* - Introduire une valeur de tolerance pour les valeurs manquantes
* - Enlever toutes references au mode missing_value et donc ...
* - Renommer a combline4 suite aux modifications de la sequence d'appel
* B. Dugas: mai 2012
* - Renommer a combline3 et ajouter les arguments MISS_ALL et
*   FILL_ALL qui indiquent respectivement au retour que le
*   champs courant a ete completer manque ou rempli
* B. Dugas: mai 2007
* - Ajouter l'argument laslon pour ajouter(+1) ou enlever (-1) une longitude
* - A cause du changement precedent on change le nom a combline2
* - Utiliser des "tolerances" pour les valeurs "miss_cdf" et "fill_cdf"
* - La fonction def_ijk est mise "inline"
* G. Bergeron: aout 2005
*       - Reorganisation des arguments d'appels
*       - Reintroduction des variables logiques en argument d'appel.
* A. Frigon: avril 2005
*     Elimine logical miss_ccc_def,fill_ccc_def car deja presents dans 
*     cdf2ccc.h
* C. Desrochers: mars 2005
*     Ajout de la possibilite de remplacer les valeurs de remplissage
* A. Frigon:  octobre 2003
*     Ajout de la possibilite de remplacer les valeurs manquantes
*
******

      real(8) fill_hold,fill_cdf_toler
      integer i,j,jmin,jmax,jstep,ijk,def_ijk,nif, fill_count

      def_ijk(i,j,k,ni,nj)=nj*ni*(k-1)+ni*(j-1)+i

      logical , external :: idnan
*-----------------------------------------------------------------------

      fill_count=0 ; fill_all=.false.
      if (fill_ok .and. .not.fill_cdf_nan)
     .    fill_cdf_toler=abs( fill_cdf ) * 0.001

      if (invjj) then
         jmin=nj
         jmax=1
         jstep=-1
      else
         jmin=1
         jmax=nj
         jstep=1
      end if

                        nif=ni
      if (laslon.eq.-1) nif=ni-1

      do j=jmin,jmax,jstep
         do i=1,nif

            indice=indice+1
            ijk=def_ijk(i,j,k,ni,nj)

            if (fill_ok .and. fill_cdf_nan) then

               if (idnan( value(ijk),.true. )) then
                  gg(indice)=fill_ccc
                  fill_ccc_oui=.true.
                  fill_count=fill_count+1
               else
                  gg(indice)=(scale*value(ijk)+offset)*mult + add
               endif

            else

               gg(indice)=(scale*value(ijk)+offset)*mult + add 

               if (fill_ok) then
                  fill_hold = abs( value(ijk) - fill_cdf )
                  if (fill_hold < fill_cdf_toler) then
                     gg(indice)=fill_ccc  
                     fill_ccc_oui=.true.
                     fill_count=fill_count+1
                  endif
               endif

            endif

         enddo

         if (laslon.eq.+1) then
            indice=indice+1
            gg(indice)=gg(indice-ni)
         endif

      enddo                     

      if (fill_count == nif*nj) fill_all=.true.

      return
*-----------------------------------------------------------------------

      end
