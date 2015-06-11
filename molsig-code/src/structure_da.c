/***********************************************************
  File for walk in data structure for signature program
  Jean-loup Faulon May 91 Albuquerque
  Modified PSU Aug. - Oct. 91, PSU Jan 92, April 93,
***********************************************************/
#include <general.h>
#include <eps.h>

LOCAL	t_err	reset_all_bond_visit(molecule,group,atom)
/**********************************************************
  A bond is visited only once.
  Created by Jean-Loup Faulon, March 1996.
***********************************************************/
molecule_p	molecule;
group_p		group;
atom_p		atom;
{
	molecule_p	m;
	set_group_p	sg;
	group_p		g;
	set_atom_p	sa;
	atom_p		a;
	set_bond_p	sb;
	bond_p		b;

	if ((molecule == NIL) && (group == NIL) && (atom == NIL)) return(ERROR);
	if (molecule != NIL) {
	   m = molecule;
	   sg = m->SG;
	   g = sg->group; 
           sa = g->SA; 
	   }
	else if (group != NIL) {
	   g = group;
           m = (molecule_p)g->molecule;
	   sg = NIL; 
           sa = g->SA;
	   }
	else {
	   a = atom;
           m = (molecule_p)((group_p)(a->group))->molecule;
	   g = (group_p)a->group;
	   sg = NIL; sa = NIL;
	   }

	LOOP { /* all groups in the molecule */
          if (sg) { g = sg->group; sa = g->SA; }
	  LOOP { /* all atoms in the group */
            if (sa) a = sa->atom;
            sb = (set_bond_p)a->SB;
	    LOOP { /* all bonds in atom */
              if (sb == NIL) ENDLOOP;
              if (sb == (set_bond_p)ERROR) ENDLOOP;
              b = sb->bond; if (b) b->visit = FALSE;
              sb = sb->succ; 
              }
            if (sa) sa = sa->succ; else ENDLOOP;
            if (sa == NIL) ENDLOOP;
	    }
	  if (sg) sg = sg->succ; else ENDLOOP;
          if (sg == NIL) ENDLOOP;
	  }
	return(OK);
	}

LOCAL	t_err	all_dihedral(bond,_function,marq,arg)
/**********************************************************
  Walk in all dihedrals of a bond. Run
  function with argument arg.
***********************************************************/
bond_p		bond;
f_traitement	_function;
t_bool		marq;
t_pointer	arg;
{
	atom_p		a1,a2,a3,a4;
	set_bond_p	sb1,sb2;
	bond_p		b1,b2;

	if (bond == NIL) return(ERROR);
	a1 = bond->atom1; a2 = bond->atom2;
	if ((a1 == NIL) || (a2 == NIL)) return(OK);

        sb1 = (set_bond_p)a1->SB;
	LOOP { /* all bonds connected to a1 */
          if (sb1 == NIL) ENDLOOP;
          if (sb1 == (set_bond_p)ERROR) ENDLOOP;
          b1 = sb1->bond; if (b1 == bond) {sb1 = sb1->succ; continue;}
	  a3 = ((b1->atom1 == a1)? b1->atom2: b1->atom1);
	  sb2 = (set_bond_p)a2->SB;
	  LOOP {
            if (sb2 == NIL) ENDLOOP;
            if (sb2 == (set_bond_p)ERROR) ENDLOOP;
            b2 = sb2->bond; if (b2 == bond) {sb2 = sb2->succ; continue;}
	    a4 = ((b2->atom1 == a2)? b2->atom2: b2->atom1);
	    if ((a3) && (a4))
  	       if (_function(a3,a1,a2,a4,arg) == ERROR) return(ERROR);
	    sb2 = sb2->succ;
	    }
          sb1 = sb1->succ; 
          }
	  return(OK);
	}

DEFINE	t_err	dast_all_dihedral(molecule,group,atom,bond,
                                  _function,marq,arg)
/**********************************************************
  Walk in all dihedrals of a molecule a group an atom
  or a bond. Run function with argument arg.
***********************************************************/
molecule_p	molecule;
group_p		group;
atom_p		atom;
bond_p		bond;
f_traitement	_function;
t_bool		marq;
t_pointer	arg;
{
	molecule_p	m;
	set_group_p	sg;
	group_p		g;
	set_atom_p	sa;
	atom_p		a;
	set_bond_p	sb;
	bond_p		b;

	if (bond) {
	   if ((marq == FALSE) || (b->marq == marq))
              return(all_dihedral(bond,_function,marq,arg));
	   return(OK);
	   }

	if ((molecule == NIL) && (group == NIL) && (atom == NIL)) return(ERROR);
	reset_all_bond_visit(molecule,group,atom);
	if (molecule != NIL) {
	   m = molecule;
	   sg = m->SG;
	   g = sg->group; 
           sa = g->SA; 
	   }
	else if (group != NIL) {
	   g = group;
           m = (molecule_p)g->molecule;
	   sg = NIL; 
           sa = g->SA;
	   }
	else {
	   a = atom;
           m = (molecule_p)((group_p)(a->group))->molecule;
	   g = (group_p)a->group;
	   sg = NIL; sa = NIL;
	   }

	LOOP { /* all groups in the molecule */
          if (sg) { g = sg->group; sa = g->SA; }
	  LOOP { /* all atoms in the group */
            if (sa) a = sa->atom;
            sb = (set_bond_p)a->SB;
	    LOOP { /* all bonds in atom */
              if (sb == NIL) ENDLOOP;
              if (sb == (set_bond_p)ERROR) ENDLOOP;
              b = sb->bond;
	      if (b->visit == FALSE)
              if ((marq == FALSE) || (b->marq == marq))
  	         if (all_dihedral(b,_function,marq,arg) == ERROR) {
	            reset_all_bond_visit(molecule,group,atom);
                    return(ERROR);
	            }
	      b->visit = TRUE;
              sb = sb->succ; 
              }
            if (sa) sa = sa->succ; else ENDLOOP;
            if (sa == NIL) ENDLOOP;
	    }
	  if (sg) sg = sg->succ; else ENDLOOP;
          if (sg == NIL) ENDLOOP;
	  }
	reset_all_bond_visit(molecule,group,atom);
	return(OK);
	}

DEFINE	t_err	dast_all_bond(molecule,group,atom,_function,
                              marq,arg)
/**********************************************************
  Walk in all bonds of a molecule a group or an atom. Run 
  function with argument arg.
  A bond is vivited only one (J. L. Faulon March 1996).
***********************************************************/
molecule_p	molecule;
group_p		group;
atom_p		atom;
f_traitement	_function;
t_bool		marq;
t_pointer	arg;
{
	molecule_p	m;
	set_group_p	sg;
	group_p		g;
	set_atom_p	sa;
	atom_p		a;
	set_bond_p	sb;
	bond_p		b;

	if ((molecule == NIL) && (group == NIL) && (atom == NIL)) return(ERROR);
	reset_all_bond_visit(molecule,group,atom);
	if (molecule != NIL) {
	   m = molecule;
	   sg = m->SG;
	   g = sg->group; 
           sa = g->SA; 
	   }
	else if (group != NIL) {
	   g = group;
           m = (molecule_p)g->molecule;
	   sg = NIL; 
           sa = g->SA;
	   }
	else {
	   a = atom;
           m = (molecule_p)((group_p)(a->group))->molecule;
	   g = (group_p)a->group;
	   sg = NIL; sa = NIL;
	   }

	LOOP { /* all groups in the molecule */
          if (sg) { g = sg->group; sa = g->SA; }
	  LOOP { /* all atoms in the group */
            if (sa) a = sa->atom;
            sb = (set_bond_p)a->SB;
	    LOOP { /* all bonds in atom */
              if (sb == NIL) ENDLOOP;
              if (sb == (set_bond_p)ERROR) ENDLOOP;
              b = sb->bond;
	      if (b->visit == FALSE)
              if ((marq == FALSE) || (b->marq == marq))
  	         if (_function(b,arg) == ERROR) {
	            reset_all_bond_visit(molecule,group,atom);
                    return(ERROR);
	            }
	      b->visit = TRUE;
              sb = sb->succ; 
              }
            if (sa) sa = sa->succ; else ENDLOOP;
            if (sa == NIL) ENDLOOP;
	    }
	  if (sg) sg = sg->succ; else ENDLOOP;
          if (sg == NIL) ENDLOOP;
	  }
	reset_all_bond_visit(molecule,group,atom);
	return(OK);
	}

LOCAL	t_err	all_angle(atom,_function,marq,arg)
/**********************************************************
  Walk in all angles of an atom. Run
  function with argument arg.
***********************************************************/
atom_p		atom;
f_traitement	_function;
t_bool		marq;
t_pointer	arg;
{
	atom_p		a1,a2;
	set_bond_p	sb1,sb2;
	bond_p		b1,b2;

	if (atom == NIL) return(ERROR);
        sb1 = (set_bond_p)atom->SB;
	LOOP { /* all triplets a1 - atom - a2 */
          if (sb1 == NIL) ENDLOOP;
          if (sb1 == (set_bond_p)ERROR) ENDLOOP;
          b1 = sb1->bond; a1 = ((b1->atom1 == atom)? b1->atom2: b1->atom1);
	  sb2 = sb1->succ;
	  LOOP {
            if (sb2 == NIL) ENDLOOP;
            if (sb2 == (set_bond_p)ERROR) ENDLOOP;
            b2 = sb2->bond; a2 = ((b2->atom1 == atom)? b2->atom2: b2->atom1);
	    if ((a1) && (a2))
  	       if (_function(a1,atom,a2,arg) == ERROR) return(ERROR);
	    sb2 = sb2->succ;
	    }
          sb1 = sb1->succ; 
          }
	  return(OK);
	}

DEFINE	t_err	dast_all_angle(molecule,group,atom,_function,marq,arg)
/**********************************************************
  Walk in all angle in a molecule a group or an atom. Run
  function with argument arg.
***********************************************************/
molecule_p	molecule;
group_p		group;
atom_p		atom;
f_traitement	_function;
t_bool		marq;
t_pointer	arg;
{
	molecule_p	m;
	set_group_p	sg;
	group_p		g;
	set_atom_p	sa;
	atom_p		a;

	if (atom) {
	   if ((marq == FALSE) || (a->marq == marq))
              return(all_angle(atom,_function,marq,arg));
	   return(OK);
	   }
	if ((molecule == NIL) && (group == NIL)) return(ERROR);
	if (molecule != NIL) {
	   m = molecule;
	   sg = m->SG; 
	   g = sg->group; 
	   sa = g->SA; 
	   }
	else  {
	   g = group;
           m = (molecule_p)g->molecule;
	   sg = NIL; 
           sa = g->SA;
	   }

	LOOP { /* all groups in the molecule */
          if (sg) { g = sg->group; sa = g->SA; }
	  LOOP { /* all atoms in the group */
            if (sa == NIL) ENDLOOP;
            a = sa->atom;
            if ((marq == FALSE) || (a->marq == marq))
	       if (all_angle(a,_function,marq,arg) == ERROR) return(ERROR);
            sa = sa->succ;
	    }
	  if (sg) sg = sg->succ; else ENDLOOP;
          if (sg == NIL) ENDLOOP;
	  }
	return(OK);
	}

DEFINE	t_err	dast_all_atom(molecule,group,_function,marq,arg)
/**********************************************************
  Walk in all atoms in a molecule or a group. Run
  function with argument arg.
***********************************************************/
molecule_p	molecule;
group_p		group;
f_traitement	_function;
t_bool		marq;
t_pointer	arg;
{
	molecule_p	m;
	set_group_p	sg;
	group_p		g;
	set_atom_p	sa;
	atom_p		a;

	if ((molecule == NIL) && (group == NIL) ) return(ERROR);
	if (molecule != NIL) {
	   m = molecule;
	   sg = m->SG; 
	   g = sg->group; 
	   sa = g->SA; 
	   }
	else  {
	   g = group;
           m = (molecule_p)g->molecule;
	   sg = NIL; 
           sa = g->SA;
	   }

	LOOP { /* all groups in the molecule */
          if (sg) { g = sg->group; sa = g->SA; }
	  LOOP { /* all atoms in the group */
            if (sa == NIL) ENDLOOP;
            a = sa->atom;
            if ((marq == FALSE) || (a->marq == marq))
	       if (_function(a,arg) == ERROR) return(ERROR);
            sa = sa->succ;
	    }
	  if (sg) sg = sg->succ; else ENDLOOP;
          if (sg == NIL) ENDLOOP;
	  }
	return(OK);
	}

DEFINE	t_err	dast_all_group(molecule,_function,marq,arg)
/**********************************************************
  Walk in all groups in a molecule. Run
  function with argument arg.
***********************************************************/
molecule_p	molecule;
f_traitement	_function;
t_bool		marq;
t_pointer	arg;
{
	molecule_p	m;
	set_group_p	sg;
	group_p		g;

        if (molecule == NIL) return(ERROR);
	m = molecule;

	sg = m->SG; 
	LOOP { /* all groups in the molecule */
          if (sg == NIL) ENDLOOP;
          g = sg->group;
          if ((marq == FALSE) || (g->marq == marq))
	     if (_function(g,arg) == ERROR) return(ERROR);
	  sg = sg->succ; 
	  }
	return(OK);
	}

DEFINE	t_err	dast_all_molecule(SM,_function,marq,arg)
/**********************************************************
  Walk in all molecule in a set. Run
  function with argument arg.
***********************************************************/
set_molecule_p	SM;
f_traitement	_function;
t_bool		marq;
t_pointer	arg;
{
	molecule_p	m;

        if (SM == NIL) return(ERROR);

	LOOP { /* all molecule in the set */
          if (SM == NIL) ENDLOOP;
          m = SM->molecule;
          if ((marq == FALSE) || (m->marq == marq))
	     if (_function(m,arg) == ERROR) return(ERROR);
	  SM = SM->succ; 
	  }
	return(OK);
	}

typedef	struct	ARG {
		group_p		group;
		t_bool		marq;
		f_traitement	f;
		t_pointer	farg;
		} arg_t, *arg_p;
LOCAL	t_err	all_bonding_site(atom,arg)
/**********************************************************
  Apply _function to all the bonding sites
***********************************************************/
atom_p		atom;
arg_p		arg;
{
	group_p		g;

        if (atom->degre - dast_number_bond_atom(atom) < 1) return(OK);
        g = (group_p)atom->group;
	if (arg->group) {
	   if (strcmp(arg->group->comment,((group_p)atom->group)->comment)) return(OK);
           }
	atom->marq = arg->marq;
	if (arg->f) arg->f(atom,arg->farg);
	return(OK);
	}

DEFINE	t_err	dast_all_bonding_site(molecule,group,atom,marq,
					_function,farg)
/**********************************************************
  Apply _function to all the bonding sites
***********************************************************/
molecule_p	molecule;
group_p		group;
atom_p		atom;
f_traitement	_function;
t_bool		marq;
t_pointer	farg;
{
	arg_t	arg;

	arg.marq = marq; arg.f = _function; arg.farg = farg;
	if (group) molecule = (molecule_p)group->molecule; arg.group = group;
	if (atom) all_bonding_site(atom,&arg);
	else dast_all_atom(molecule,NIL,all_bonding_site,FALSE,&arg);
	return(OK);
	}

LOCAL	t_err	set_bond_marq(b,m)
/**********************************************************
 no comment.
***********************************************************/
bond_p	b;
t_bool	m;
{	if (b) b->marq = m; return(OK); 
	}

LOCAL	t_err	set_atom_marq(a,m)
/**********************************************************
 no comment.
***********************************************************/
atom_p	a;
t_bool	m;
{	if (a) a->marq = m; return(OK); 
	}

LOCAL	t_err	set_group_marq(g,m)
/**********************************************************
 no comment.
***********************************************************/
group_p	g;
t_bool	m;
{	if (g) g->marq = m; return(OK); 
	}

DEFINE	t_err	dast_set_marq(molecule,group,atom,marq)
/**********************************************************
  Set the marq for all elements of a molecule a group 
  or an atom.
***********************************************************/
molecule_p	molecule;
group_p		group;
atom_p		atom;
t_bool		marq;
{
	if (molecule) molecule->marq = marq;
        if (molecule) dast_all_group(molecule,set_group_marq,FALSE,marq);
	if (group)    group->marq = marq;
        if ((molecule) || (group))
           dast_all_atom(molecule,group,set_atom_marq,FALSE,marq);
        if (atom)     atom->marq = marq;
        dast_all_bond(molecule,group,atom,set_bond_marq,FALSE,marq);
	return(OK);
	}

typedef struct ARG_EXIST {
	bond_p	bond;
	atom_p	atom;
	} arg_exist_t, *arg_exist_p;
LOCAL	t_err	exist_bond(b,arg)
/**********************************************************
 ERROR if bond b is connected to a
***********************************************************/
bond_p		b;
arg_exist_p	arg;
{
	if ((b == NIL) || (arg == NIL)) return(OK);

	arg->bond = b;
	if (b->atom1 == arg->atom) return(ERROR);
	if (b->atom2 == arg->atom) return(ERROR);
	return(OK);
	}

DEFINE	bond_p	dast_inq_bond(atom1,atom2)
/**********************************************************
  TEST and return if a bond exist beetween atom1 and atom2.
***********************************************************/
atom_p	atom1,atom2;
{
	arg_exist_t	arg;
	bond_p		bond1 = NIL, bond2 = NIL;

	arg.atom = atom2; arg.bond = NIL; 
        if (dast_all_bond(NIL,NIL,atom1,exist_bond,FALSE,&arg) == ERROR) 
        bond1 = arg.bond;
	arg.atom = atom1; arg.bond = NIL;  
        if (dast_all_bond(NIL,NIL,atom2,exist_bond,FALSE,&arg) == ERROR) 
        bond2 = arg.bond;
	return((bond1 == bond2) ? bond1 : NIL);
	}

DEFINE	bond_p	dast_inq_bond_directed(atom1,atom2)
/**********************************************************
  TEST and return if a bond exist from atom1 to atom2.
***********************************************************/
atom_p	atom1,atom2;
{
	arg_exist_t	arg;
	bond_p		bond1 = NIL, bond2 = NIL;

	arg.atom = atom2; arg.bond = NIL; 
        if (dast_all_bond(NIL,NIL,atom1,exist_bond,FALSE,&arg) == ERROR) 
        bond1 = arg.bond;
	return(bond1);
	}

DEFINE	t_err	dast_inq_bonding_site(atom1,atom2,bds1,bds2)
/**********************************************************
  TEST  if a bond site exist in atom1 and atom2.
  Return the two bonding site.
***********************************************************/
atom_p	atom1,atom2;
bond_p	*bds1,*bds2;

{
	arg_exist_t		arg;

	if (bds1) *bds1 = NIL; if (bds2) *bds2 = NIL;
	if (atom1) {
	arg.atom = NIL; arg.bond = NIL; 
        if (dast_all_bond(NIL,NIL,atom1,exist_bond,FALSE,&arg) == ERROR) 
        *bds1 = arg.bond;
	}
	if (atom2) {
	arg.atom = NIL; arg.bond = NIL;  
        if (dast_all_bond(NIL,NIL,atom2,exist_bond,FALSE,&arg) == ERROR) 
        *bds2 = arg.bond;
	}
        if ((atom1) && (*bds1 == NIL)) return(ERROR);
        if ((atom2) && (*bds2 == NIL)) return(ERROR);
	return(OK);
	}

        
