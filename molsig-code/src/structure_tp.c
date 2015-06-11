/***********************************************************
  File that deals with topological problems on the
  data structures.
  Jean-loup Faulon May 91 Albuquerque
  Modified PSU Aug. - Oct. 91, PSU Jan 92, April 93,
***********************************************************/
#include <general.h>
#include <eps.h>

typedef struct	ARG_CONNECT 	{
	t_bool		marq;
	bond_p		bond;
	f_traitement	f;
	t_pointer	farg;
	t_integer	distance,hight,number;
	} arg_connect_t, *arg_connect_p;

LOCAL	t_err	connectivity_bond(bond,arg)
/**********************************************************
  Determine and marks all bonds connected to the initial 
  atom (different from arg->bond).
***********************************************************/
bond_p		bond;
arg_connect_p	arg;
{
EXTERN	t_err	connectivity_atom();
	if (bond == arg->bond) return(OK);
	if (arg->marq) if (bond->marq == arg->marq) return(OK);
	bond->marq = arg->marq;
	connectivity_atom(bond->atom1,arg->marq,arg->bond,
                          arg->f,arg->farg,arg->distance,arg->hight,&arg->number);
        connectivity_atom(bond->atom2,arg->marq,arg->bond,
                          arg->f,arg->farg,arg->distance,arg->hight,&arg->number);
	return(OK);
	}

LOCAL	t_err	connectivity_bond2(bond,arg)
/**********************************************************
  Determine and mark TRUE all bonds connected to the initial 
  atom (different from arg->bond).
***********************************************************/
bond_p		bond;
arg_connect_p	arg;
{
EXTERN	t_err	connectivity_atom2();
	if (bond == arg->bond) return(OK);
	if (bond->marq) return(OK);
	bond->marq = TRUE;
	connectivity_atom2(bond->atom1,arg->bond,
                          arg->f,arg->farg,arg->distance,arg->hight,&arg->number);
        connectivity_atom2(bond->atom2,arg->bond,
                          arg->f,arg->farg,arg->distance,arg->hight,&arg->number);
	return(OK);
	}

LOCAL	t_err	connectivity_bond_degree(bond,arg)
/**********************************************************
  Determine and mark TRUE all bonds connected to the initial 
  atom (different from arg->bond).
***********************************************************/
bond_p		bond;
arg_connect_p	arg;
{
EXTERN	t_err	connectivity_atom_degree();
	if (bond == arg->bond) return(OK);
	if (bond->marq) return(OK);
	bond->marq = TRUE;
	connectivity_atom_degree(bond->atom1,
                          arg->f,arg->farg,arg->distance,&arg->number);
	connectivity_atom_degree(bond->atom2,
                          arg->f,arg->farg,arg->distance,&arg->number);
	return(OK);
	}

LOCAL	t_err	connectivity_bond_unmarq(bond)
/**********************************************************
  Unmarq bonds
***********************************************************/
bond_p		bond;
{
EXTERN	t_err	dast_connectivity_atom_unmarq();

	if (bond->marq == FALSE) return(OK);
	bond->marq = FALSE;
	dast_connectivity_atom_unmarq(bond->atom1);
        dast_connectivity_atom_unmarq(bond->atom2);
	return(OK);
	}

DEFINE	t_err	connectivity_atom(atom,marq,bond,
                                  f,farg,distance,hight,number)
/**********************************************************
  Determine and mark all atoms connected to atom.
  Without using the way by bond. The distance between this
  atom an the initial one must be less than distance.
  The function f is applied to atom.
  BUGG: distance not coded properly will work with 
        distance = 1 or 2 and MAXATOM
***********************************************************/
atom_p		atom;
t_bool		marq;
bond_p		bond;
f_traitement	f;
t_pointer	farg;
t_integer	distance,hight,*number;
{
	arg_connect_t	arg;
	if (atom == NIL) return(OK);
	if (marq) if (atom->marq == marq) return(OK);
        if (hight > distance) return(OK);
        if (f) farg = (t_pointer)f(atom,farg,hight); /* PSU Jan 92 */
        (*number)++; atom->marq = marq; 
        arg.marq = marq; arg.bond = bond; arg.f = f; arg.farg = farg; 
        arg.distance = distance; arg.hight = hight+1; arg.number = *number;
	dast_all_bond(NIL,NIL,atom,connectivity_bond,FALSE,&arg);
	*number = arg.number;
	return(OK);
	}

DEFINE	t_err	connectivity_atom2(atom,bond,f,farg,distance,hight,number)
/**********************************************************
  Determine and marq TRUE all atoms connected to atom.
  Without using the way by bond. The distance between this
  atom an the initial one must be less than distance.
  The function f is applied to atom.
  BUGG: distance not coded properly will work with 
        distance = 1 or 2 and MAXATOM
***********************************************************/
atom_p		atom;
bond_p		bond;
f_traitement	f;
t_pointer	farg;
t_integer	distance,hight,*number;
{
	arg_connect_t	arg;
	t_bool		marq = FALSE;

	if (atom == NIL) return(OK);
	if (atom->marq == TRUE) return(OK);
        if (hight > distance) return(OK);
        if (f) f(atom,farg,hight); /* Jan 97 */
        (*number)++; atom->marq = TRUE;
        arg.marq = marq; arg.bond = bond; arg.f = f; arg.farg = farg; 
        arg.distance = distance; arg.hight = hight+1; arg.number = *number;
	dast_all_bond(NIL,NIL,atom,connectivity_bond2,FALSE,&arg);
	*number = arg.number;
	return(OK);
	}

DEFINE	t_err	connectivity_atom_degree(atom,f,farg,distance,number)
/**********************************************************
  Determine and marq TRUE all atoms connected to atom.
  having a degree <= degree. 
  The function f is applied to atom.
***********************************************************/
atom_p		atom;
f_traitement	f;
t_pointer	farg;
t_integer	distance,*number;
{
	arg_connect_t	arg;
	t_bool		marq = FALSE;

	if (atom == NIL) return(OK);
	if (atom->marq == TRUE) return(OK);
	if (atom->degre > distance) return(OK);
        if (f) f(atom,farg); /* Jan 97 */
        (*number)++; atom->marq = TRUE;
        arg.marq = marq; arg.bond = NIL; arg.f = f; arg.farg = farg; 
        arg.distance = distance; arg.hight = 0; arg.number = *number;
	dast_all_bond(NIL,NIL,atom,connectivity_bond_degree,FALSE,&arg);
	*number = arg.number;
	return(OK);
	}

DEFINE	t_err	dast_connectivity_atom_unmarq(atom)
/**********************************************************
  Unmarq all marqued atoms attached to the initial atom
***********************************************************/
atom_p		atom;
{

	if (atom == NIL) return(OK);
	if (atom->marq == FALSE) return(OK);
	atom->marq = FALSE;
	dast_all_bond(NIL,NIL,atom,connectivity_bond_unmarq,FALSE);
	return(OK);
	}


LOCAL	t_err	connectivity_atom_marq(atom,f,farg)
/**********************************************************
  Apply function to all atoms attached to atom 
  that are marked
***********************************************************/
atom_p		atom;
f_traitement	f;
t_pointer	farg;
{
	set_bond_p	sb;
	bond_p		b;
	atom_p		a;

	if (atom == NIL) return(OK);
	if (atom->marq == FALSE) return(OK);
        if (f) farg = (t_pointer)f(atom,farg,0); 
	atom->marq = FALSE;
	sb = (set_bond_p)atom->SB;
	LOOP { /* all bonds in atom */
          if (sb == NIL) ENDLOOP;
          if (sb == (set_bond_p)ERROR) ENDLOOP;
          b = sb->bond;
	  a = ((b->atom1 == atom) ? b->atom2 : b->atom1);
          connectivity_atom_marq(a,f,farg);
          sb = sb->succ; 
          }
	return(OK);
	}

LOCAL	t_err	connectivity_atom_distance(atom,marq,bond,
                                  f,farg,distance,hight,number)
/**********************************************************
  Determine and mark all atoms connected to atom.
  The marq is the distance from the original atom.
  The function f is applied to atom.
***********************************************************/
atom_p		atom;
t_bool		marq;
bond_p		bond;
f_traitement	f;
t_pointer	farg;
t_integer	distance,hight,*number;
{
	set_bond_p	sb;
	bond_p		b;
	atom_p		a;

	if (atom == NIL) return(OK);
	if (marq) if (atom->marq)
	   if (hight > atom->marq) return(OK);
        if (hight > distance) return(OK);
        if (f) farg = (t_pointer)f(atom,farg,hight); /* PSU Jan 92 */
        if (atom->marq == 0) (*number)++; atom->marq = hight; 
	sb = (set_bond_p)atom->SB;
	LOOP { /* all bonds in atom */
          if (sb == NIL) ENDLOOP;
          if (sb == (set_bond_p)ERROR) ENDLOOP;
          b = sb->bond; b->marq = hight; 
	  a = ((b->atom1 == atom) ? b->atom2 : b->atom1);
          connectivity_atom_distance(a,marq,bond,\
                                     f,farg,distance,hight+1,number);
          sb = sb->succ; 
          }
	return(OK);
	}

DEFINE	t_err	dast_connectivity_atom(atom,marq,bond,f,farg,distance)
/**********************************************************
  Determine and marks all atoms connected to atom 
  without using the way by bond. The distance must be 
  less than "distance". The function f is applied
  to each atom.
  BUGG: distance not coded properly will work with 
        distance = 1 or 2 and MAXATOM
***********************************************************/
atom_p		atom;
t_bool		marq;
bond_p		bond;
f_traitement	f;
t_pointer	farg;
t_integer	distance;
{
	t_integer	number = 0;

	if (atom == NIL) return(OK);
        connectivity_atom(atom,marq,bond,f,farg,distance,0,&number);
	return(number);
	}

DEFINE	t_err	dast_connectivity_atom_marq(atom,f,farg)
/**********************************************************
  Apply function to all atoms attached to atom 
  that are marked
***********************************************************/
atom_p		atom;
f_traitement	f;
t_pointer	farg;
{
	if (atom == NIL) return(OK);
	connectivity_atom_marq(atom,f,farg);
	return(OK);
	}


DEFINE	t_err	dast_connectivity_atom2(atom,bond,f,farg,distance)
/**********************************************************
  Determine but does not marks all atoms connected to atom 
  without using the way by bond. The distance must be 
  less than "distance". The function f is applied
  to each atom.
  BUGG: distance not coded properly will work with 
        distance = 1 or 2 and MAXATOM
***********************************************************/
atom_p		atom;
bond_p		bond;
f_traitement	f;
t_pointer	farg;
t_integer	distance;
{
	t_integer	number = 0;

	if (atom == NIL) return(OK);
	connectivity_atom2(atom,bond,f,farg,distance,0,&number);
	dast_connectivity_atom_unmarq(atom);
	return(number);
	}

DEFINE	t_err	dast_connectivity_atom_degree(atom,f,farg,degree)
/**********************************************************
  Determine but does not marks all atoms connected to atom 
  having a degree <= degree. The function f is applied
  to each atom.
***********************************************************/
atom_p		atom;
f_traitement	f;
t_pointer	farg;
t_integer	degree;
{
	t_integer	number = 0;

	if (atom == NIL) return(OK);
	connectivity_atom_degree(atom,f,farg,degree,&number);
	dast_connectivity_atom_unmarq(atom);
	return(number);
	}

DEFINE	t_err	dast_connectivity_atom_distance(atom,marq,bond,f,farg,distance)
/**********************************************************
  Determine and marks all atoms connected to atom 
  without using the way by bond. The distance must be 
  less than "distance". The function f is applied
  to each atom.
***********************************************************/
atom_p		atom;
t_bool		marq;
bond_p		bond;
f_traitement	f;
t_pointer	farg;
t_integer	distance;
{
	t_integer	number = 0;

	if (atom == NIL) return(OK);
	connectivity_atom_distance(atom,marq,bond,f,farg,distance+1,1,&number);
	return(number);
	}

DEFINE	t_bool	dast_connected_atom(a1,a2,bond)
/**********************************************************
 True if a1 and a2 are connected without bond
***********************************************************/
atom_p		a1,a2;
bond_p		bond;
{
	molecule_p	molecule;
	t_bool		answer = FALSE;

	if ((a1 == NIL) || (a2 == NIL)) return(FALSE);
	molecule = (molecule_p)((group_p)a1->group)->molecule;
	dast_connectivity_atom(a1,TRUE,bond,NIL,NIL,MAXINTEGER);
	if (a1->marq == a2->marq) answer = TRUE;
	dast_set_marq(molecule,NIL,NIL,FALSE);
	return(answer);
	}

#define	MAXCYCLESIZE	14
LOCAL	t_err	HMIN = MAXCYCLESIZE+1;
LOCAL	t_err	path_atom(a1,a2,bond,h)
/**********************************************************
 Return the path length between a1 and a2 are connected 
 without using bond
***********************************************************/
atom_p		a1,a2;
bond_p		bond;
t_integer	h;
{
	molecule_p	molecule;
	t_integer	l,lmin = MAXCYCLESIZE+1;
	atom_p		a = a1;
	set_bond_p	s;

	if ((a->marq) && (a->marq < h+1)) return(MAXCYCLESIZE+1);
	if (h > MAXCYCLESIZE) return(MAXCYCLESIZE+1);
	if (a == a2) return(0);
	a->marq = h+1;
	s = (set_bond_p)a->SB;           
	while (s) {
	  bond_p	b = s->bond; 
	  if (b == NIL) break;
	  if (b == bond) {s = s->succ; continue;}
          a = ((b->atom1 == a1) ? b->atom2 : b->atom1);
	  l = path_atom(a,a2,bond,h+1) + 1;
	  if (l < lmin) lmin = l;
	  s = s->succ;
          }
	return(lmin);
	}

DEFINE	t_err	dast_path_atom(a1,a2,bond)
/**********************************************************
 Return the path length between a1 and a2 are connected 
 without using bond
***********************************************************/
atom_p		a1,a2;
bond_p		bond;
{
	molecule_p	molecule;
	t_integer	l = 0;

	if ((a1 == NIL) || (a2 == NIL)) return(ERROR);
	molecule = (molecule_p)((group_p)a1->group)->molecule;
	l = path_atom(a1,a2,bond,0);
	dast_set_marq(molecule,NIL,NIL,FALSE);
	return((l<MAXCYCLESIZE)?l:-1);
	}

DEFINE   t_bool  dast_cycle_atom(atom)
/**********************************************************************
  TRUE if atom is in a cycle.
  Atom is in a cycle if any of its neighbour can be reached
  without using the bond linking the two atoms.
**********************************************************************/
atom_p          atom;
{
	set_bond_p	sb;
	bond_p		b;
	atom_p		a;

	if (atom == NIL) return(FALSE);
	if (!(dast_number_bond_atom(atom) > 1)) return(FALSE); 

	sb = (set_bond_p)atom->SB;

	LOOP { /* all bonds in atom */
          if (sb == NIL) ENDLOOP;
          if (sb == (set_bond_p)ERROR) ENDLOOP;
          b = sb->bond;
	  a = ((b->atom1 == atom) ? b->atom2 : b->atom1);
	  if (dast_number_bond_atom(a) > 1) 
	     if (dast_connected_atom(atom,a,b)) return(TRUE);
          sb = sb->succ; 
          }
	return(FALSE);
        }

DEFINE   t_bool  dast_cycle_atom_pair(molecule,ax,ay,b)
/**********************************************************************
  TRUE if pair ax-ay add a cycle.
**********************************************************************/
molecule_p      molecule;
atom_p          ax,ay;
bond_p          b;
{
	t_bool	RT;
        dast_set_marq(molecule,NIL,NIL,FALSE);
        dast_connectivity_atom(ax,'x',b,NIL,NIL,MAXATOM);
        if (ay->marq == ax->marq) RT = TRUE; else RT = FALSE;
        dast_set_marq(molecule,NIL,NIL,FALSE);
        return(RT);
        }


LOCAL   t_bool  inq_cycle_atom(atom,bond,L,l,CYCLE)
/**********************************************************************
  Recursive routine finding any cycle attached to atom of length L
**********************************************************************/
atom_p          atom;
bond_p		bond;
t_integer	L,l;
atom_p		CYCLE[];
{
	set_bond_p	sb;
	bond_p		b;
	atom_p		a;

	if (dast_number_bond_atom(atom) < 2) return(FALSE);
	if (atom->marq) return(FALSE);

	if (l == L) {
	   if (atom == CYCLE[0]) return(TRUE); return(FALSE);
	   }
	atom->marq = TRUE;
	sb = (set_bond_p)atom->SB;

	LOOP { /* all bonds of atom */
          if (sb == NIL) ENDLOOP;
          if (sb == (set_bond_p)ERROR) ENDLOOP;
          b = sb->bond;
	  if (b == bond) {sb = sb->succ; continue;}
	  a = ((b->atom1 == atom) ? b->atom2 : b->atom1);
	  if (inq_cycle_atom(a,b,L,l+1,CYCLE)) {
	     CYCLE[l] = atom; return(TRUE);
             }
          sb = sb->succ; 
          }
	return(FALSE);
        }

DEFINE   t_bool  dast_inq_cycle_atom(atom,L,CYCLE)
/**********************************************************************
  TRUE if atom is in a cycle of lenght L, FALSE otherwise.
  If TRUE any cycle of lenght L is returned in CYCLE[0]..CYCLE[L-1].
**********************************************************************/
atom_p          atom;
t_integer	L;
atom_p		CYCLE[];
{
	set_bond_p	sb;
	bond_p		b;
	atom_p		a;
	t_integer	i;
	t_bool		err = FALSE;
	molecule_p	molecule;

	if (atom == NIL) return(FALSE);
	if (!(dast_number_bond_atom(atom) > 1)) return(FALSE); 
	if (L < 3) return(FALSE);
	molecule = (molecule_p)((group_p)atom->group)->molecule;
        dast_set_marq(molecule,NIL,NIL,FALSE);
	for (i = 0; i < L ; i++) CYCLE[i] = NIL;
	sb = (set_bond_p)atom->SB;
	CYCLE[0] = atom;

	LOOP { /* all bonds of atom */
          if (sb == NIL) ENDLOOP;
          if (sb == (set_bond_p)ERROR) ENDLOOP;
          b = sb->bond;
	  a = ((b->atom1 == atom) ? b->atom2 : b->atom1);
	  if ((err = inq_cycle_atom(a,b,L,1,CYCLE))) break;
          sb = sb->succ; 
          }
        dast_set_marq(molecule,NIL,NIL,FALSE);
	return(err);
        }

        
