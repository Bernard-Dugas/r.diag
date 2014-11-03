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
      subroutine rdspec2 (NCID,FUNIT)

      implicit none

      include 'cdf2ccc.h'
      include 'dimmem.h'
      include 'infomem.h'
      include 'varmem.h'   
      include 'ibufmem.h'
      include 'specmem.h'
      include 'workmem.h'  

      integer ncid,funit

******
*
*AUTEUR Guy Bergeron       juin 2003
*
*     Traitement des variables definies dans l'espace spectral.
*
*REVISIONS
*
*  Bernard Dugas octobre 2014 :
*  - Declarations/initialisations locales des variables ktr et lmt
*  Bernard Dugas fevrier 2014 :
*  - Faire appel a COMBLINE5 (ajouter FILL_CDF_NAN=.false.)
*  Bernard Dugas juin 2013 :
*  - Faire appel a COMBLINE4 (enlever les miss_*)
*  Bernard Dugas aout 2012 :
*  - Remplacer COMBLINE2 par COMBLINE3 (deux nouveaux arguments)
*  B.Dugas mai 2012 : 
*  - Ne plus initialiser zdid ici, mais tout de meme,
*    lui donner la valeur de zid s'il le faut
*  B.Dugas automne 2007 : 
*  - Ne plus utiliser HPALLOC/HPDEALLC
*  - Remplacer la routine FXLRLMT par la fonction CLRLMT
*  - Tenir compte des IP1 dans IBUF(4) pour les fichiers CMC/RPN
*  G. Bergeron aout 2005 :  Reorganisation des arguments d'appels de combline
*  Guy Bergeron juin  2005 : parametre d'appel pour invj
*  Guy Bergeron avril 2004 : Declaration de coordonne en REAL*8
*
******

******netCDF :

      integer id,len
      real*8   miss_cdf,fill_cdf
      integer  start(maxdim),count(maxdim)
      logical  miss_all,fill_all

******CCCma :

      character cccname*4
      integer   i,n,nn,k,kk,first,indice,itime
      integer   ccctime,la,lm,lr,lrlmt, ktr,lmt

      real    bad

      integer ilevel(maxlev)

      integer  clrlmt
      external clrlmt

******
      data miss_cdf,fill_cdf /0.0,0.0/
      logical var_ok
*-----------------------------------------------------------------------
      if (invj) then
         write(6,5999)
         call                                        xit('rdspec2', -1)
      endif

      write(6,6000)

      nn=0

      do id=1,nvars

         if (list(id)%var_ok ) then

            list(id)%var_ok=.false.
            nn=nn+1
         
*           Definir le type derive var(id)

            call affecte_var(nn,list(id)%name,list(id)%type,
     .                      list(id)%ndim,list(id)%dimid,list(id)%nattr)
            
*           Initialisation :

            do i=1,ndims
               start(i)=1             
               count(i)=1
            enddo
            if (zdid < 1) zdid = zid

            count(numdid)=dim(numdid)%len            ! nbres coefficient spectraux
            count(zdid)=dim(zdid)%len                ! nbres niveaux


*           Lire les attributs :

            call get_attribut(ncid,id,var(nn)%nattr,var(nn)%name)
            call def_cccma(var(nn)%name,var(nn)%mult,var(nn)%add,
     .                                                      cccname,bad)   

*           Variables definissants la troncature (ktr) et le nombre d'onde (lmt):

            lmt=0 ; ktr=-1

            do n=1,var(nn)%nattr
               if(attr(n)%name.eq.'trunc_type') then
                  if (attr(n)%cvalue.eq.'Triangular') ktr=2
                  if (attr(n)%cvalue.eq.'Rhumboidal') ktr=0
               endif
               if(attr(n)%name.eq.'trunc_count') lmt=attr(n)%i2value(1)
            enddo

            lrlmt = clrlmt( lmt+1,lmt+1,ktr,.true. )    ! definit lrlmt
            call dimgt(ival,la,lr,lm,ktr,lrlmt)         ! definir la

            call def_level(ilevel,'encode')             ! definir les etiquettes ibuf(4)

*           Lire les valeurs de la variable :

            write(6,6020)
            do itime=1,dim(timedid)%len                  ! boucle temporelle

               start(timedid)=itime

               len=dim(numdid)%len*dim(zdid)%len
               call get_varda(ncid,id,var(nn)%type,start,count,len,
     .                                                   variable(1,nn))

*              Traitement des donnees :

               do kk=dim(zdid)%len,1,-1              ! Commence en haut

                  k=dim(zdid)%len-kk+1
                  indice=0
                  first=dim(level2did)%len*(itime-1)+2*kk-1

                  call combline5( variable(1,nn),dval,indice,kk,
     .                                 dim(numdid)%len,1,dim(zdid)%len,
     .                                   scale_fact(kk),add_offset(kk),
     .                                var(nn)%mult,var(nn)%add,.false.,
     .                             fill_cdf,.false.,0,fill_all,.false. )
     .                                               

                  dval(1)=mean(first)
                  dval(2)=mean(first+1)
                  

*                 Ecriture dans le fichier de sortie :
                  
                  call decodate(ncid,dcoordonne(itime,tid),ccctime)
                  call setlab(ibuf,'SPEC',ccctime,cccname,ilevel(k),
     .                                                 la,1,lrlmt,npack)

                  call putfld2(funit,dval,ibuf,maxpk)
                  if (itime.eq.dim(timedid)%len)
     .                                      call prtlab( ibuf )

               enddo            !level
            enddo               !time

         endif                  !var_ok

      enddo                     !variable
*-----------------------------------------------------------------------
 5999 format(/," L'INVERSION DE L'INDICE J POUR LE CAS SPECTRAL N'EST",/
     .         " PAS POSSIBLE !",/)
 6000 format(/,/,' VARIABLES :')
 6020 format(/)
*-----------------------------------------------------------------------
      end
