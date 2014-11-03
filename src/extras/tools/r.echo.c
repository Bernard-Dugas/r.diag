
/*

     Ce programme renvoie sur stdout la premiere chaine de caractere lue
     en argument ssi deux conditions sont verifiees:

     1) La variable d'environment DIAGNOSTIC_EXIT n'est pas definie ou bien,
        si elle l'est, le fichier auquel elle pointe est vide ou inexistant.
        De plus si la variable DIAGNOSTIC_EXIT n'est pas definie, le fichier
        DIAGNOSTIC_EXIT ne doit lui-meme pas exister ou bien, etre vide.

     2) La variable d'environment ECHO_COMMAND_LINE est definie non-vide.

     (par James Caveen, RPN/SI, 2 juillet 1993)

     Mis a jour le 22 decembre 2003 par Yves Chartier, pour AIX.

 */

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
     FILE *fichier;
     unsigned char *ptr_echo_command_line, *ptr_diagnostik_exit;

/*
 *   sauter le nom du programme
 */
       argv++;
/*
 *   verifier l'environnement
 */
       if ((ptr_diagnostik_exit = (unsigned char *)getenv("DIAGNOSTIC_EXIT")) == (unsigned char *) NULL) 
          ptr_diagnostik_exit="DIAGNOSTIC_EXIT";

       if(( fichier = fopen(ptr_diagnostik_exit,"r")) != (FILE *) NULL) 
       {
          if(fgetc(fichier) != EOF)
                return 0;
       }
          
       
       if((ptr_echo_command_line = (unsigned char *)getenv("ECHO_COMMAND_LINE")) != (unsigned char *) NULL)
       {     
       
           while (*argv)
           {
              printf("%s\n",*argv);
              argv++;
           }
        }
       return 0;
}
