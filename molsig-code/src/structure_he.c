/***********************************************************
  MODULE DE GESTION DES PILES ET DES TAS
       Jean-loup Faulon March 1991
       Modified :
                  May 91 Albuquerque NM 
	          Nov 91 PSU PA
***********************************************************/
#include <c.h>

EXTERN	t_pointer	hea_creer(),hea_alloc();

typedef        struct  PILE   {
               t_pointeur      info;
               struct PILE    *succ;
               } pile_t, *pile_p;

DEFINE t_err  heat_print_pile(pile)
/**********************************************************
 Print pile
***********************************************************/
pile_p         pile;
{
       t_integer	size = 0;
       pile_p  		p = pile;

       if (pile  == ((pile_p)ERREUR)) return(0);
       if (pile  == NIL) return(0);
       while ((p != (pile_p)ERROR) && (p)) { printf(" [@%X=%d]",p,p->info); ; p = p->succ; }
       printf("\n");
       return(size);
}

DEFINE t_pointeur      heat_creer_hea(taille)
/**********************************************************
***********************************************************/
t_compteur     taille;
{
       t_pointeur     TAS; 
LOCAL  t_compteur     taille_defaut = 4088;

       if (taille > 0)   taille_defaut = taille;
       TAS = (t_pointeur)hea_creer(taille_defaut,0);
       return(TAS);
       }


DEFINE t_pointeur      heat_creer_hea_pile(taille)
/**********************************************************
***********************************************************/
t_compteur     taille;
{
       t_pointeur     TAS; 
LOCAL  t_compteur     taille_defaut;

       if (taille > 0)   taille_defaut = taille;
       TAS = (t_pointeur)hea_creer(taille_defaut,0);
       return(TAS);
       }

DEFINE t_err   heat_detruire_hea_pile(TAS)
/**********************************************************
***********************************************************/
t_pointeur     TAS; 
{
       if (TAS == NIL) return(OK);
       return((t_err)hea_detruire(TAS));
       }


DEFINE t_err   heat_detruire_hea(TAS)
/**********************************************************
***********************************************************/
t_pointeur     TAS; 
{
       if (TAS == NIL) return(OK);
       return((t_err)hea_detruire(TAS));
       }


DEFINE pile_p  heat_creer_pile()
/**********************************************************
***********************************************************/
{
       return((pile_p)ERREUR);
       }

DEFINE t_bool  heat_vide_pile(pile)
/**********************************************************
***********************************************************/
pile_p         pile;
{
       if (pile == NIL) return(VRAI);
       if (pile == ((pile_p)ERREUR)) return(VRAI);
       return(FAUX);
}

DEFINE pile_p  heat_empile(element,pile,TAS)
/**********************************************************
 Insertion d'un element dans une pile
***********************************************************/
t_pointeur     element;
pile_p         pile;
t_pointeur     TAS;
{
       pile_p  p;

       if (TAS  == ((t_pointeur)ERREUR)) return((pile_p)ERREUR);
       
       p = (pile_p)hea_alloc(TAS,sizeof(pile_t));
       p->info = element;
       p->succ = pile;

       return(p);
       }

DEFINE t_bool  heat_inq_pile(element,pile)
/**********************************************************
 TRUE if element is in pile
***********************************************************/
t_pointeur     element;
pile_p         pile;
{
       pile_p  p = pile;

       while ((p != NIL) && (p!= (pile_p)ERROR)) {
         if (p->info == element) return(TRUE);
	 p = p->succ;
	 }
       return(FALSE);
       }

DEFINE pile_p  heat_insert_pile(element,pile,TAS)
/**********************************************************
 Insert element iff element is not already in pile
***********************************************************/
t_pointeur     element;
pile_p         pile;
t_pointeur     TAS;
{
       pile_p  p = pile;

       if (TAS  == ((t_pointeur)ERREUR)) return((pile_p)ERREUR);

       /* is element in pile ? */
       while ((p != NIL) && (p!= (pile_p)ERROR)) {
         if (p->info == element) return(pile);
	 p = p->succ;
	 }
       p = (pile_p)hea_alloc(TAS,sizeof(pile_t));
       p->info = element;
       p->succ = pile;

       return(p);
       }

DEFINE pile_p  heat_remove_pile(element,pile)
/**********************************************************
 Remove element from pile
***********************************************************/
t_pointeur     element;
pile_p         pile;
{
       pile_p  p = pile, pp = NIL;

       /* is element in pile ? */
       while ((p != NIL) && (p!= (pile_p)ERROR)) {
         if (p->info == element) {
	    if (pp == NIL) return(p->succ);
	    pp->succ = p->succ; return(pile);
	    }    
	 pp = p; p = p->succ;
	 }
       return(pile);
       }

DEFINE pile_p  heat_depile(pile)
/**********************************************************
 Depiler
***********************************************************/
pile_p         pile;
{
       if (pile == ((pile_p)ERREUR)) return((pile_p)ERREUR);
       if (pile == NIL) return((pile_p)ERROR);
       return(pile->succ);
}

DEFINE t_err  heat_size_pile(pile)
/**********************************************************
 Nombre d'element d'une pile
 PSU Nov 91.
***********************************************************/
pile_p         pile;
{
       t_integer	size = 0;
       pile_p  		p = pile;

       if (pile  == ((pile_p)ERREUR)) return(0);
       if (pile  == NIL) return(0);
       while ((p != (pile_p)ERROR) && (p)) { size++; p = p->succ; }
       return(size);
}

DEFINE t_pointer  heat_info_pile(pile)
/**********************************************************
 Info
***********************************************************/
pile_p         pile;
{
       if (pile == ((pile_p)ERREUR)) return((t_pointer)ERREUR);
       if (pile == NIL) return((t_pointer)ERROR);
       return(pile->info);
}

DEFINE pile_p  heat_elem_pile(pile,n)
/**********************************************************
 Return element n in the stack
***********************************************************/
pile_p         pile;
t_integer      n;
{
       t_integer	i;

       if (pile == ((pile_p)ERREUR)) return((pile_p)ERREUR);
       if (pile == NIL) return((pile_p)ERROR);
       if (n > heat_size_pile(pile)) return((pile_p)ERREUR);
       for (i = 0; i < n; i++) pile = pile->succ;
       return(pile);
}


DEFINE pile_p  heat_reverse_pile(pile)
/*********************************************************
 Reverse the order of a stack.
***********************************************************/
pile_p         pile;
{
       t_integer	size = 0,i,n;
       pile_p  		p1,p2;
       t_pointer        info;

       if (pile == ((pile_p)ERREUR)) return(pile);
       if (pile == NIL) return(pile);
       n = heat_size_pile(pile)-1;
       for (i = 0; i <= (n-1)/2; i++) {
	   p1 = heat_elem_pile(pile,i);
	   p2 = heat_elem_pile(pile,n-i);
	   info = p1->info; p1->info = p2->info; p2->info = info;
	   }
       return(pile);
}

DEFINE pile_p  heat_substract_pile(pile1,pile2,TAS)
/**********************************************************
 Return (pile2-pile1)
***********************************************************/
pile_p         pile1,pile2;
t_pointeur     TAS;
{
       pile_p  pile = heat_creer_pile();

       if (TAS  == ((t_pointeur)ERREUR)) return((pile_p)ERREUR);

       /* Insert all element of pile2 in pile that ARE NOT in pile1 */
       while (heat_vide_pile(pile2) == FALSE) {
	  if (heat_inq_pile(pile2->info,pile1) == FALSE)
	     heat_empile(pile2->info,pile,TAS);
	  pile2 = heat_depile(pile2);
	  }
       return(pile);
       }

DEFINE pile_p  heat_intersect_pile(pile1,pile2,TAS)
/**********************************************************
 Return (pile2 INTER pile1)
***********************************************************/
pile_p         pile1,pile2;
t_pointeur     TAS;
{
       pile_p  pile = heat_creer_pile();

       if (TAS  == ((t_pointeur)ERREUR)) return((pile_p)ERREUR);

       /* Insert all element of pile2 in pile that ARE in pile1 */
       while (heat_vide_pile(pile2) == FALSE) {
	  if (heat_inq_pile(pile2->info,pile1) == TRUE)
	     heat_empile(pile2->info,pile,TAS);
	  pile2 = heat_depile(pile2);
	  }
       return(pile);
       }

DEFINE pile_p  heat_merge_pile(pile1,pile2)
/**********************************************************
 Merge pile1 in pile2. return(pile2->pile1)
***********************************************************/
pile_p         pile1,pile2;
{
       pile_p  p = pile2;
       if ((pile1 == NIL) || (pile1 == (pile_p)ERROR)) return(pile2);
       if ((pile2 == NIL) || (pile2 == (pile_p)ERROR)) return(pile1);
       while ((p->succ != NIL) && (p->succ != (pile_p)ERROR)) p = p->succ;
       p->succ = pile1;
       return(pile2);
       }

