
/*
     $Log: divers.c,v $
     Revision 3.11  2013/02/07 21:38:21  bernard
     Ajouter la routine ETIME (getrusage).

     Revision 3.10  2010/02/10 18:22:22  dugas
     Ajouter CSORTL (pour entiers 'long long').

     Revision 3.9  2008/04/25 20:43:41  dugas
     Enlever le 'include <string.h>' et simplifier les 'define's.

     Revision 3.8  2004/11/08 20:42:31  dugas
     Ajouter la routine swap_endianness

     Revision 3.7  2003/10/24 21:05:48  dugas
     Implementer du code compatible RS6000

     Revision 3.6  2002/08/26 18:15:09  dugas
      Modification de nice19 suite a une suggestion
      de MFV de verifier le contenu de VERY_NICE_DIAG.

     Revision 3.5  2002/08/20 18:51:16  dugas
     Ajouter les fonctions appellables de fortran findlowcoreindex
        et findhighcoreindex, de meme que leurs versions en c, telles
        que fournies par Y.Chartier.

     Revision 3.4  1999/04/08 19:17:55  armnrbd
     Ajouter un macro pour identifier le code LINUX.

     Revision 3.3  1999/01/19 20:45:26  armnrbd
     Renommer SCOPY par DSCOPY.

 * Revision 3.2  1995/10/10  15:56:12  armnrbd
 * Remplacer le macro C4680 par le macro IRIX5.
 *
 * Revision 3.1  1994/11/17  14:14:18  armnrbd
 * Messages informatifs quand au passage de la version 2.x a 3.1...
 * 1) Les espaces en debut des noms de variables de sont plus pertinents.
 * 2) Les grilles complexes de type CMPL sont maintenant supportees.
 * 3) Les fichiers SQI sont reconnus, lus et ecrit directements.
 * 4) Plusieurs nouvelles cles sont disponibles au demarrage.
 *
 * Revision 3.0  94/11/17  13:56:22  13:56:22  armnrbd (Bernard Dugas)
 * *** empty log message ***
 * 
 * Revision 2.2  94/01/18  16:05:52  armnrbd
 * Ajouter la fonction NICE19.
 * 
 * Revision 2.1  94/01/12  11:33:17  armnrbd
 * Ajouter un "if defined (C4680)" pour cpp.
 * 
 * Revision 2.0  93/10/13  13:32:24  armnrbd
 * Premiere version compatible HP-UX.
 * 
 * Revision 1.9  93/08/23  11:49:27  11:49:27  armnrbd (Bernard Dugas)
 * Utiliser "rpnmacros.h" plutot que "macros.h".
 * 
 * Revision 1.8  93/08/23  10:20:42  armnrbd
 * Correction a GETENVC.
 * 
 * Revision 1.7  93/02/18  10:47:09  armnrbd
 * Permettre a GETENVC d'enlever les blancs en bout de parametres.
 * 
 * Revision 1.6  92/10/13  15:37:37  armnrbd
 * Enlever la definition du macro wordfloat.
 * 
 * Revision 1.5  92/06/19  13:43:33  armnrbd
 * Modifier la documentation a CSORT.
 * 
 * Revision 1.4  92/06/18  15:44:06  armnrbd
 * Ajouter des nouvelles versions de csort.
 * Ajouter le support du NEC (old-style paramaters).
 * 
 * Revision 1.3  92/04/02  11:01:34  armnrbd
 * Ajouter les routines CSORTR (reel) et CSORTD (double).
 * 
 * Revision 1.2  92/03/04  16:05:31  armnrbd
 * Correction a getenvc.
 * 
 * Revision 1.1  92/03/04  13:57:43  armnrbd
 * Ajouter la routine getenvc.
 * 
 * Revision 1.0  92/02/21  12:06:49  armnrbd
 * Initial revision
 * 
*/


#include <stdlib.h>
#include <rpnmacros.h>
#include <string.h>

/* 

 Les trois fichiers suivants etaient inclus dans
 la version YRC des fonctions findlowcoreindex
 et findhighcoreindex

 include <string.h>
 include <stdio.h>

*/

#if defined (HP)

#include <sys/times.h>
#include <unistd.h>

/*

   Fonction retournant la somme des temps usager et systeme d'un
   travail sur une machine HP.   ( BD, RPN - 08 octobre 1993 )

*/

wordfloat
f77name(second) ( )

{

   struct tms buffer;
   clock_t elapsed;
   wordint ticks;
   wordfloat hold;

   ticks = sysconf(_SC_CLK_TCK) ;
   elapsed = times(&buffer) ;
   hold = buffer.tms_utime + buffer.tms_stime ;

   return hold / ticks ;
   
}

#endif

/*

   Routine qui force un equivalent de "nice -19".

   Modifiee le 26 aout 2002 suite a une suggestion
   de MFV de verifier le contenu de VERY_NICE_DIAG.

*/

wordint
f77name(nice19) ( )

{
   char *nicefact =  (char *) getenv("VERY_NICE_DIAG");

   nicefact = nicefact ? nicefact : "19";
   return nice( atoi( nicefact ) ? atoi( nicefact ) : 19 );

}

/*

   Emulation de CVMGT du CRAY.  Les deux premiers 
   arguments sont traites en reel simple precision.

*/

wordfloat 
f77name(cvmgt) ( arg1, arg2, logic )
wordfloat *arg1, *arg2 ;
wordint   *logic ;

{
   return *logic ? *arg1 : *arg2 ;
}


/*

   La routine leadblk enleve les espaces pouvant se
   trouver au debut de la chaine de caracteres variab.
   Auteur: b. dugas - 31 aout 1994.

*/

void f77name(leadblk)( name,len )
wordint len;
char name[1];

{

   wordint i,j;
   unsigned wordint size;
   char *temp;

   size = len+1 ;
   temp = (char *) malloc( size ) ;

   for ( i=0 ; 
         i < len && name[i] == ' ' ; 
         i++ ) ;

   if ( i < len ) {

       for ( j=i     ; j < len ; j++ ) *(temp+j-i) = name[j] ;
       for ( j=len-i ; j < len ; j++ ) *(temp+j)   = ' ' ;
       for ( j=0     ; j < len ; j++ ) *(name+j)   = temp[j] ;

       }

   free( temp );

}

/*

   La routine getenvc fouille l'environnement pour une chaine
   de type "name=xxx", retournant "xxx" ou une chaine vide 
   dans value. ( BD, RPN - 04 mars 1992 )

*/

void
f77name(getenvc) ( name, value, len1, len2 )
wordint  len1, len2;
char name[1], value[1];

{

   wordint i;
   unsigned wordint size;
   char *temp, *hold;

/* Transfert de name dans une chaine C */

   size = len1+len2+1 ;
   temp = (char *) malloc( size ) ;

   for ( i=0 ; 
         i < len1 && name[i] != ' ' ; 
         i++ ) 
         *(temp+i) = name[i] ;

   *(temp+i) = '\0' ;

/* Appel a la fonction C getenv */

   hold = (char *) getenv( temp ) ;

/*
   Si la chaine n'a pas ete trouvee, getenv retourne un
   pointeur null. value sera alors entierement vide.
   Sinon, on copie le contenu de hold dans value.
 
*/

   for ( i=0 ; i < len2 ; i++ ) value[i] = ' ' ;

   if ( hold != 0 && size != 1 ) {
        size = strlen( hold ) ;
        for ( i=0 ; i < size ; i++ ) value[i] = *(hold+i) ;
        }

   free( temp );

}

/*

   La routine qqexit retourne un code EXIT au process parent.

*/

void 
f77name(qqexit) (val)
wordint *val ;

{
     exit(*val);
}

/*

   Emulation de la famille de routines CRAY ISRCHF{XX}
   ou XX peut etre LT,LE,GT ou GE.

   Appel FORTRAN ...

      CALL ISRCHF{XX}( n,sx,inc,sy )

   Parametres ...

      n   - Nombre d'elements a verifier. 
            Valeur de retour 0, si n <= 0.
      sx  - Tableau source de nombres reels (float).
      inc - Increments entre les elements successifs.
      sy  - Valeur comparative reelle.

*/

wordint
f77name(isrchflt) ( n, sx, inc, sy )
wordint    *n, *inc ;
wordfloat  *sx, *sy ;

{

   wordint count ;

   if ( *n <= 0 || *inc <= 0 ) return 0 ;

   for ( count = -(*n) ; count++ && (*sx >= *sy) ; )  sx += *inc ;

   return ( count + *n -1 ) * *inc + 1 ;

}

wordint
f77name(isrchfle) ( n, sx, inc, sy )
wordint    *n, *inc ;
wordfloat  *sx, *sy ;

{

   wordint count ;

   if ( *n <= 0 || *inc <= 0 ) return 0 ;

   for ( count = -(*n) ; count++ && (*sx > *sy) ; )  sx += *inc ;

   return ( count + *n -1 ) * *inc + 1 ;

}

wordint
f77name(isrchfgt) ( n, sx, inc, sy )
wordint    *n, *inc ;
wordfloat  *sx, *sy ;

{

   wordint count ;

   if ( *n <= 0 || *inc <= 0 ) return 0 ;

   for ( count = -(*n) ; count++ && (*sx <= *sy) ; )  sx += *inc ;

   return ( count + *n -1 ) * *inc + 1 ;

}

wordint
f77name(isrchfge) ( n, sx, inc, sy )
wordint    *n, *inc ;
wordfloat  *sx, *sy ;

{

   wordint count ;

   if ( *n <= 0 || *inc <= 0 ) return 0 ;

   for ( count = -(*n) ; count++ && (*sx < *sy) ; )  sx += *inc ;

   return ( count + *n -1 ) * *inc + 1 ;

}

/*

   Emulation de la routine SCOPY du CRAY.

   Appel FORTRAN ...

      CALL DSCOPY( n,sx,incx,sy,incy )

   Parametres ...

      n    - Nombre d'elements a copier. 
      sx   - Tableau source de nombres reels (float).
      incx - Increments entre les elements successifs de sx.
      sy   - Tableau destination de nombres reels (float).
      incy - Increments entre les elements successifs de sy.

*/

void 
f77name(dscopy) ( n, sx, incx, sy, incy  )
wordint   *n,  *incx, *incy ;
wordfloat      *sx,   *sy ;

{
   wordint count ;

   if ( *n > 0 ) 
   {

      count = *n ;

      while ( count-- )
      {
         *sy = *sx ;
         sy += *incy ;
         sx += *incx ;
      }

   }

}

#define SHRINKFACTOR 1.3

/* 

   Fonction Combsort11 extraite de BYTE, volume 16 numero 4
   (avril 1991), pp. 315-320,  par Richard Box et Stephen Lacey. 
   Modifiee par B. Dugas, RPN - 24 avril 1991.

      CALL CSORTC( list,index,size )  ...  version caracteres
      CALL CSORTD( list,index,size )  ...  version nombres reels "double"
      CALL CSORTE( list,index,size )  ...  version nombres entiers
      CALL CSORTL( list,index,size )  ...  version nombres entiers "long long"
      CALL CSORTR( list,index,size )  ...  version nombres reels "single"

   Parametres ...

      list  - Valeurs a trier (respectivement, de type caractere, 
              entier ou reel (simple ou double), selon que l'on fasse
              appel a csortc, csorte, csortr ou csortd).
      index - Valeurs entieres de 1 a "size". On ne deplace pas
              les elements de "list", mais plutot leurs indices.
      size  - Nombre entier d'elements dans "list" et "index".


   Notes ...

      Les routines CSORTE, CSORTD et CSORTR sont a toutes fins identiques.
      La seule difference reside dans la declaration des types contenus
      dans list.

*/

void 
f77name(csortr) ( list, index, size )
wordfloat *list ;
wordint *index, *size ;

{

   int switches = 0 ;
   wordint gap = *size, i = 0, j = 0, top = 0 ;

/* Faire le tri en comparant les elements de "list" */

   do
      {

      gap = (wordint)((wordfloat)gap/SHRINKFACTOR);
      switch (gap)
	 {
	 case 0:  /* The smallest gap is 1 - bubble sort */
	    gap = 1;
	    break;
	 case 9:  /* This is what makes this Combsort11 */
	 case 10:
 	    gap = 11;
	    break;
	 default:
	    break;
	 }

      switches = 0; /* Dirty Pass flag */
      top = *size-gap;
	    
      for ( i=0 ; i < top ; ++i )
	 {
	 j = i+gap;
	 if ( *(list + *(index + i) - 1) > *(list + *(index + j) - 1) )
	    { /* Swap indices */
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    *(index+j) = *(index+i) ^ *(index+j) ;
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    ++switches;
	    } /* End of Swap */
	 }    /* End of Pass */

      } while ( switches || (gap > 1) );
/* 
   Like the bubble sort and shell, we check for
   a clean pass, i.e. switches = 0 and gap = 1.
*/

}


void 
f77name(csortd) ( list, index, size )
double *list ;
wordint *index, *size ;

{

   int switches = 0 ;
   wordint gap = *size, i = 0, j = 0, top = 0 ;

/* Faire le tri en comparant les elements de "list" */

   do
      {

      gap = (wordint)((wordfloat)gap/SHRINKFACTOR);
      switch (gap)
	 {
	 case 0:  /* The smallest gap is 1 - bubble sort */
	    gap = 1;
	    break;
	 case 9:  /* This is what makes this Combsort11 */
	 case 10:
 	    gap = 11;
	    break;
	 default:
	    break;
	 }

      switches = 0; /* Dirty Pass flag */
      top = *size-gap;
	    
      for ( i=0 ; i < top ; ++i )
	 {
	 j = i+gap;
	 if ( *(list + *(index + i) - 1) > *(list + *(index + j) - 1) )
	    { /* Swap indices */
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    *(index+j) = *(index+i) ^ *(index+j) ;
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    ++switches;
	    } /* End of Swap */
	 }    /* End of Pass */

      } while ( switches || (gap > 1) );
/* 
   Like the bubble sort and shell, we check for
   a clean pass, i.e. switches = 0 and gap = 1.
*/

}


void 
f77name(csorte) ( list, index, size )
wordint *list, *index, *size ;

{

   int switches = 0 ;
   wordint gap = *size, i = 0, j = 0, top = 0 ;

/* Faire le tri en comparant les elements de "list" */

   do
      {

      gap = (wordint)((wordfloat)gap/SHRINKFACTOR);
      switch (gap)
	 {
	 case 0:  /* The smallest gap is 1 - bubble sort */
	    gap = 1;
	    break;
	 case 9:  /* This is what makes this Combsort11 */
	 case 10:
 	    gap = 11;
	    break;
	 default:
	    break;
	 }

      switches = 0; /* Dirty Pass flag */
      top = *size-gap;
	    
      for ( i=0 ; i < top ; ++i )
	 {
	 j = i+gap;
	 if ( *(list + *(index + i) - 1) > *(list + *(index + j) - 1) )
	    { /* Swap indices */
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    *(index+j) = *(index+i) ^ *(index+j) ;
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    ++switches;
	    } /* End of Swap */
	 }    /* End of Pass */

      } while ( switches || (gap > 1) );
/* 
   Like the bubble sort and shell, we check for
   a clean pass, i.e. switches = 0 and gap = 1.
*/

}


void 
f77name(csortl) ( list, index, size )
long long *list ;
wordint *index, *size ;

{

   int switches = 0 ;
   wordint gap = *size, i = 0, j = 0, top = 0 ;

/* Faire le tri en comparant les elements de "list" */

   do
      {

      gap = (wordint)((wordfloat)gap/SHRINKFACTOR);
      switch (gap)
	 {
	 case 0:  /* The smallest gap is 1 - bubble sort */
	    gap = 1;
	    break;
	 case 9:  /* This is what makes this Combsort11 */
	 case 10:
 	    gap = 11;
	    break;
	 default:
	    break;
	 }

      switches = 0; /* Dirty Pass flag */
      top = *size-gap;
	    
      for ( i=0 ; i < top ; ++i )
	 {
	 j = i+gap;
	 if ( *(list + *(index + i) - 1) > *(list + *(index + j) - 1) )
	    { /* Swap indices */
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    *(index+j) = *(index+i) ^ *(index+j) ;
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    ++switches;
	    } /* End of Swap */
	 }    /* End of Pass */

      } while ( switches || (gap > 1) );
/* 
   Like the bubble sort and shell, we check for
   a clean pass, i.e. switches = 0 and gap = 1.
*/

}

void 
f77name(csortc) ( list, index, size, len )
char *list ;
wordint *index, *size, len ;

{

   int switches = 0 ;
   wordint gap = *size, i = 0, j = 0, top = 0 ;

/* Faire le tri en comparant les elements de "list" */

   do
      {

      gap = (wordint)((wordfloat)gap/SHRINKFACTOR);
      switch (gap)
	 {
	 case 0:  /* The smallest gap is 1 - bubble sort */
	    gap = 1;
	    break;
	 case 9:  /* This is what makes this Combsort11 */
	 case 10:
 	    gap = 11;
	    break;
	 default:
	    break;
	 }

      switches = 0; /* Dirty Pass flag */
      top = *size-gap;
	    
      for ( i=0 ; i < top ; ++i )
	 {
	 j = i+gap;
	 if ( strncmp ( list + (*(index + i) - 1)*len,
                        list + (*(index + j) - 1)*len, len ) > 0 )
	    { /* Swap indices */
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    *(index+j) = *(index+i) ^ *(index+j) ;
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    ++switches;
	    } /* End of Swap */
	 }    /* End of Pass */

      } while ( switches || (gap > 1) );

/* Like the bubble sort and shell, we check for a clean pass */

} 

/*
   The following routine simply calls
   the rpnmacro swap_buffer_endianness
*/

void 
f77name(swap_endianness) ( buffer,size )
wordint *buffer, *size ;

{

  swap_buffer_endianness( buffer,*size );

}

/*
   The following two functions were borrowed from Yves Chartier
   in July 2002, and are also used in his version of delamineur
*/

#include <math.h>

int f77name(findlowcoreindex)(float ax[], int *ni);
int f77name(findhighcoreindex)(float ax[], int *ni);

int c_findLowCoreIndex(float ax[], int ni);
int c_findHighCoreIndex(float ax[], int ni);

int f77name(findlowcoreindex)(float ax[], int *ni)
{
  int lc;
  lc =  c_findLowCoreIndex(ax, *ni);
  return lc+1;
}

int f77name(findhighcoreindex)(float ax[], int *ni)
{
  int hc;
  hc = c_findHighCoreIndex(ax, *ni);
  return hc+1;
}


int c_findLowCoreIndex(float ax[], int ni)
{
  int is;
  float dx, dxref;
  
  is = ni/2;
  dxref = ax[is] - ax[is-1];
  dx = ax[is-1] - ax[is-2];
  is--;
  while (((fabs(dxref-dx)) < (0.001 * dxref)) && (is > 0))
    {
    is--;
    dx = ax[is] - ax[is-1];
    }
  if (is == 1) is = 0;
  return is;
}

int c_findHighCoreIndex(float ax[], int ni)
{
  int is;
  float dx, dxref;

  is = ni/2;
  dxref = ax[is+1] - ax[is];
  dx = ax[is+2] - ax[is+1];
  is++;
  while (((fabs(dxref-dx)) < (0.001 * dxref)) && (is < (ni-1)))
    {
    is++;
    dx = ax[is+1] - ax[is];
    }
  if (is == (ni-2)) is=ni-1;
  return is;
}

/*
   The following is an alternative to the etime fortran
   extension provided by Michel Valin in February 2013
*/

#include <sys/time.h>
#include <sys/resource.h>

#pragma weak etime_=etime
void etime(float *timearray, float *time)
{
  struct rusage usage;
  int status=getrusage(RUSAGE_SELF,&usage);
timearray[0]=usage.ru_utime.tv_sec+usage.ru_utime.tv_usec*.000001 ;
timearray[1]=usage.ru_stime.tv_sec+usage.ru_stime.tv_usec*.000001 ;
  *time=timearray[0]+timearray[1];
} 
