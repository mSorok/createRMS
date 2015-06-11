#include <stdarg.h> /* WMB */
#include <c.h>
#include "hea.h"
 
EXTERN	t_pointer	malloc();

DEFINE t_err  hea_max(taille)
/************************************************************
        fixe pour le systeme la taille maximale du hea ]
        cf set_sbreak_size
        on ne peut l'appeler qu'une fois
         */
t_taille         taille;
{
      LOCAL    t_bool  prem= VRAI;

      if (prem != VRAI ) {
            xerr_c("hea max already fixed");
            return(ERROR);
            }
/*    else  set_sbrk_size(taille) ;  */ 
      prem = FAUX;
      return(OK);
      }
 
p_hea  hea_creer(t_taille taille_bloc,t_taille nb_std_cel, ...)
/************************************************************
         Si taille_bloc <= TAS_TAIL_BLOC_DEFAUT  alors on prend TAS_TAIL_BLOC_DEFAUT
         ATTENTION  la taille maximale d'une cellule est la taille du bloc/2
         La taille minimale d'une cellule standard est la taille de t_pointeur 
         Un type std est repere par sont rang dans la liste (de o a nb_std_cel)
         */
/* t_taille         taille_bloc; WMB */ 
/* t_taille         nb_std_cel; WMB */
/* va_dcl  WMB */ 
{
      p_hea          hea;
      t_taille       l_hea,l_bloc;
      p_bloc         bloc;
      t_taille       i,taille_std_i;
      t_pointeur     pt;
      va_list        liste_std;

      l_bloc = (taille_bloc < TAS_TAIL_BLOC_DEFAUT ? TAS_TAIL_BLOC_DEFAUT : taille_bloc+sizeof(t_bloc));
      l_hea = nb_std_cel*sizeof(t_std_cel) + sizeof(t_hea);
      pt = (t_pointeur)malloc(l_hea+l_bloc+2*sizeof(t_bloc));
                   /* pour le premier bloc */
      if (pt == NIL) {
             xerr_c("no more space");
             return(NIL);
             }
      hea = (p_hea) pt;
      hea->taille_bloc = l_bloc;
      hea->nb_std = nb_std_cel;
                  /* on initalise le tableau des std avec la liste des 
                     tailles en argument */
      va_start(liste_std,nb_std_cel);
      for(i=0 ; i< nb_std_cel ; i++ ) {
            taille_std_i = va_arg(liste_std,t_taille);
            if ( taille_std_i < sizeof(t_cellule)) {
	        printf("%d < %d\n",taille_std_i,sizeof(t_cellule));
                xerr_c ("heap too small");
                hea->std[i].taille = sizeof(t_cellule);
                }
            else hea->std[i].taille = taille_std_i;
            hea->std[i].pile = NIL;
            }
                   /* premier bloc */
      bloc = (p_bloc)(pt + l_hea);
      bloc->suc = NIL;
      bloc->deblib = (t_pointeur)(bloc+2);
      bloc->finlib = (t_pointeur)bloc ;
      bloc->finlib += l_bloc ;
      hea->pile = bloc;
      return(hea);
      }



DEFINE t_pointeur   hea_alloc(hea,taille)
/****************************************
         alloue une cellule de taille donnee non liberable
         La taille est modifier pour les pb de Word-alignement
         (taille est pair)
         */
p_hea            hea;
t_taille         taille;
{
      t_pointeur  cel;
      p_bloc      bloc=hea->pile;

      TAS_WORD_ALIGN(taille); taille += taille % 8;
      cel = bloc->deblib;
      bloc->deblib += taille;
      if (bloc->deblib <= bloc->finlib) {}
      else {
           unsigned long l_bloc=hea->taille_bloc;

           if (taille > l_bloc/2) {
	        printf("%d > %d\n",taille,l_bloc/2);
                xerr_c ("heap too big");
                bloc->deblib -= taille;
                return(NIL);
                }
           bloc = (p_bloc)malloc(l_bloc);
           if (bloc== NIL) {
                xerr_c("no more space");
                return(NIL);
                }
           bloc->suc = hea->pile;
           bloc->deblib = (t_pointeur)(bloc+2);
           bloc->finlib = (t_pointeur)bloc ;
           bloc->finlib += l_bloc ;
           hea->pile = bloc;
           cel = bloc->deblib;
           bloc->deblib += taille;
           }
       return(cel);
       }

DEFINE p_cellule   hea_std_alloc(hea,std)
/****************************************
         alloue une cellule standart de numero std (rang de la cellule
         dans la liste des tailles de hea_creer de 0 a n)
          */
p_hea            hea;
t_taille         std;
{
      p_cellule  cel;
      
      if (std >= hea->nb_std ) {
           xerr_c("no standard defined");
           return(NIL);
           }
      if (hea->std[std].pile != NIL) {
           cel = hea->std[std].pile;
           hea->std[std].pile = cel->suc;
           }
      else cel = (p_cellule)hea_alloc(hea,hea->std[std].taille);
      return( cel);
      }

DEFINE t_err   hea_std_free(hea,std,cel)
/****************************************
         libere une cellule standart numero std 
         on verifie que la cellule apartient au hea mais pas son type
         */
p_hea            hea;
t_taille         std;
p_cellule        cel;
{
      p_bloc      b;
      
      if (std >= hea->nb_std ) {
           xerr_c("no standard defined");
           return(NIL);
           }
      for (b =hea->pile ; b != NIL ; b = b->suc) {
             if (cel <(p_cellule) b || (p_cellule)(b->finlib) < cel ) continue;
             else {   /* la cellule apartient bien a un bloc du hea */
                  cel->suc  = hea->std[std].pile;
                  hea->std[std].pile = cel;
                  return(OK);
                  }
             }
                      /* la cellule n`apartient pas au hea !!!! */
      xerr_c("heap out of range");
      return(ERROR);
      }

DEFINE p_hea   hea_raz(hea)
/****************************************
         libere la place occupee par le hea  
         sauf le premeir bloc .toutes les cellules sont
         liberees Le hea est vide
         */
p_hea            hea;
{
      p_bloc       bloc=hea->pile;
      t_taille     i;

      LOOP{
            p_bloc    b=bloc->suc;

            if (b == NIL) ENDLOOP;
                   /* le dernier bloc de la pile est le premier cree
                      il a ete alloue avec le hea */
            free(bloc);
            bloc = b;
            }
      bloc->suc = NIL;
      bloc->deblib = (t_pointeur)(bloc+2);
      bloc->finlib = (t_pointeur)bloc ;
      bloc->finlib += hea->taille_bloc ;
      hea->pile = bloc;
      for(i=0 ; i< hea->nb_std ; i++) 
            hea->std[i].pile = NIL;
      return(hea);
      }



DEFINE t_err   hea_detruire(hea)
/****************************************
         libere la place occupee par le hea 
         */
p_hea            hea;
{
      p_bloc       bloc=hea->pile;

      LOOP{
            p_bloc    b=bloc->suc;

            if (b == NIL) ENDLOOP;
                   /* le dernier bloc de la pile est le premier cree
                      il a ete alloue avec le hea */
            free(bloc);
            bloc = b;
            }
      hea->pile = NIL;
      free(hea);
      return(OK);
      }
 

