/***********************************************************
  File of data structure allocation for signature program
  Jean-loup Faulon May 91 Albuquerque
  MODIFIED JULY 92 Albq : INTRODUCTION OF COPY CONNECTIVITY
  MODIFIED OCT. 92 Albq : INTRODUCTION OF BONDING SITES
                          IN CREATE AND DELETE BOND
  MODIFIED SEP. 96 Albp:  THE *_large ROUTINES WERE ADDED
			  IN ORDER TO DUPLICATE LARGE 
			  SYSTEMS (> 10,000 ATOMS) 
***********************************************************/
#include <general.h>

EXTERN 	t_pointer 	dast_inq_bond();
EXTERN 	t_pointer 	hea_alloc();
EXTERN	t_pointer	heat_empile();
EXTERN	t_pointer	datw_inq_twin();

LOCAL	t_integer	m_ID = 1,g_ID,a_ID = 1;
LOCAL	atom_p		A = NIL, ATOM = NIL;


DEFINE	bond_p	daal_create_bond(TAS,a1,a2,b1,b2,geom1,geom2,
                                 order,comment,marq)
/**********************************************************
  Creation of a bond in the heap memory TAS.
  The bond is add in the two sets of bonds of a1 and a2
  If one of the two atoms a1 and a2 is NIL the bond is
  a bonding site (half bond) and the degre of the 
  corresponding atom is not changed.
***********************************************************/
t_pointer	TAS; 
atom_p		a1,a2;
bond_p		b1,b2; /* the two bonding sites */
t_coord		geom1,geom2;
t_bool		order;
t_bufstring	comment;
t_bool		marq;
{
EXTERN	t_err	daal_delete_bond();
LOCAL	t_coord	ZERO = {0,0,0};
       	bond_p		bond = NIL;
	group_p		g1,g2;

       	if (TAS  == ((t_pointer)ERROR)) return((bond_p)ERROR);
        if ((a1 == NIL) && (a2 == NIL)) return((bond_p)ERROR);

	/* the bond is created only if it doesn't exist already */
        if ((bond = (bond_p)dast_inq_bond(a1,a2))) return(bond);
        bond = (bond_p)hea_alloc(TAS,sizeof(bond_t));

        bond->atom1 = a1; bond->atom2 = a2;
        if (geve_equal_vector(&geom1,&ZERO) == OK)
	if ((a1) && (a2))  
             geve_compute_vector(&a1->geom,&a2->geom,&geom1);
        else geom1 = ZERO;

        if (geve_equal_vector(&geom2,&ZERO) == OK)
	if ((a1) && (a2)) 
	     geve_compute_vector(&a2->geom,&a1->geom,&geom2);   
        else geom2 = ZERO;
 	bond->geom1 = geom1; bond->geom2 = geom2;

        bond->order = order;
        if (a1) { 
           g1 = (group_p)a1->group; 
           a1->SB = (t_pointer)heat_empile(bond,a1->SB,TAS); 
           }
        if (a2) {
	   g2 = (group_p)a2->group; 
           a2->SB = (t_pointer)heat_empile(bond,a2->SB,TAS); 
           }
	strcpy(bond->comment,comment); 
	bond->marq = marq; bond->visit = FALSE;

	/* delete the bonding sites */
	daal_delete_bond(b1); daal_delete_bond(b2);
	return(bond);
        }

DEFINE	bond_p	daal_create_bond_directed(TAS,a1,a2,b1,b2,geom1,geom2,
                                 order,comment,marq)
/**********************************************************
  Creation of a directed bond in the heap memory TAS.
  The bond is added in the setof bonds of a1 
***********************************************************/
t_pointer	TAS; 
atom_p		a1,a2;
bond_p		b1,b2; /* the two bonding sites */
t_coord		geom1,geom2;
t_bool		order;
t_bufstring	comment;
t_bool		marq;
{
EXTERN	t_err	daal_delete_bond();
LOCAL	t_coord	ZERO = {0,0,0};
       	bond_p		bond = NIL;
	group_p		g1,g2;

       	if (TAS  == ((t_pointer)ERROR)) return((bond_p)ERROR);
        if ((a1 == NIL) || (a2 == NIL)) return((bond_p)ERROR);

	/* the bond is created only if it doesn't exist already */
        if ((bond = (bond_p)dast_inq_bond_directed(a1,a2))) return(bond);
        bond = (bond_p)hea_alloc(TAS,sizeof(bond_t));

        bond->atom1 = a1; bond->atom2 = a2;
        if (geve_equal_vector(&geom1,&ZERO) == OK)
	if ((a1) && (a2))  
             geve_compute_vector(&a1->geom,&a2->geom,&geom1);
        else geom1 = ZERO;

        if (geve_equal_vector(&geom2,&ZERO) == OK)
	if ((a1) && (a2)) 
	     geve_compute_vector(&a2->geom,&a1->geom,&geom2);   
        else geom2 = ZERO;
 	bond->geom1 = geom1; bond->geom2 = geom2;

        bond->order = order;
        if (a1) { 
           g1 = (group_p)a1->group; 
           a1->SB = (t_pointer)heat_empile(bond,a1->SB,TAS); 
           }
	strcpy(bond->comment,comment); 
	bond->marq = marq; bond->visit = FALSE;

	/* delete the bonding sites */
	daal_delete_bond(b1); daal_delete_bond(b2);
	return(bond);
        }

DEFINE	bond_p	daal_create_bond_large(TAS,a1,a2,b1,b2,geom1,geom2,
                                 order,comment,marq)
/**********************************************************
  Creation of a bond in the heap memory TAS.
  The bond is add in the two sets of bonds of a1 and a2
  If one of the two atoms a1 and a2 is NIL the bond is
  a bonding site (half bond) and the degre of the 
  corresponding atom is not changed.
***********************************************************/
t_pointer	TAS; 
atom_p		a1,a2;
bond_p		b1,b2; /* the two bonding sites */
t_coord		geom1,geom2;
t_bool		order;
t_bufstring	comment;
t_bool		marq;
{
EXTERN	t_err	daal_delete_bond();
LOCAL	t_coord	ZERO = {0,0,0};
       	bond_p		bond = NIL;
	group_p		g1,g2;

       	if (TAS  == ((t_pointer)ERROR)) return((bond_p)ERROR);
        if ((a1 == NIL) && (a2 == NIL)) return((bond_p)ERROR);

	/* the bond is created only if it doesn't exist already */
        if ((bond = (bond_p)dast_inq_bond(a1,a2))) return(bond);
        if ((bond = (bond_p)malloc(sizeof(bond_t))) == NIL) {
	   xerr_c("malloc no more space");
	   }
        bond->atom1 = a1; bond->atom2 = a2;
        if (geve_equal_vector(&geom1,&ZERO) == OK)
	if ((a1) && (a2))  
             geve_compute_vector(&a1->geom,&a2->geom,&geom1);
        else geom1 = ZERO;

        if (geve_equal_vector(&geom2,&ZERO) == OK)
	if ((a1) && (a2)) 
	     geve_compute_vector(&a2->geom,&a1->geom,&geom2);   
        else geom2 = ZERO;
 	bond->geom1 = geom1; bond->geom2 = geom2;

        bond->order = order;
        if (a1) { 
           g1 = (group_p)a1->group; 
           a1->SB = (t_pointer)heat_empile(bond,a1->SB,TAS); 
           }
        if (a2) {
	   g2 = (group_p)a2->group; 
           a2->SB = (t_pointer)heat_empile(bond,a2->SB,TAS); 
           }
	strcpy(bond->comment,comment); 
	bond->marq = marq; bond->visit = FALSE;

	/* delete the bonding sites */
	daal_delete_bond(b1); daal_delete_bond(b2);
	return(bond);
        }

DEFINE	t_err	daal_delete_bond(bond)
/**********************************************************
  Delete a bond, and create the bonding sites.
***********************************************************/
bond_p	bond;
{
LOCAL	t_coord	ZERO = {0,0,0};
	set_bond_p	sbs,sbp;
	t_pointer	HEAP = NIL;

	if (bond == NIL) return(ERROR);

	/* delete bond in the SB's lists */
        /* in b->atom1 */
        if (bond->atom1) {
           sbp = (set_bond_p)bond->atom1->SB; sbs = sbp;
           while (sbs) {
                 if (sbs->bond == bond) break;
                 sbp = sbs; sbs = sbs->succ;
                 }
	   if (sbs != NIL) {
              sbp->succ = sbs->succ;
              if (sbp == sbs) bond->atom1->SB = (t_pointer)sbp->succ;
              }
	   bond->atom1 = NIL;
           }
        /* in b->atom2 */
        if (bond->atom2) {
           sbp = (set_bond_p)bond->atom2->SB; sbs = sbp;
           while (sbs) {
                 if (sbs->bond == bond) break;
                 sbp = sbs; sbs = sbs->succ;
                 }
	   if (sbs != NIL) {
              sbp->succ = sbs->succ;
              if (sbp == sbs) bond->atom2->SB = (t_pointer)sbp->succ;
              }
	   bond->atom2 = NIL;
           }

        /* create the bonding sites */
	if (bond->atom1 == NIL) return(OK);
	if (bond->atom2 == NIL) return(OK);
	HEAP = (t_pointer)((molecule_p)((group_p)bond->atom1->group)->molecule)->HEAP;
	daal_create_bond(HEAP,bond->atom1,NIL,NIL,NIL,bond->geom1,ZERO,1,"",FALSE);
	daal_create_bond(HEAP,NIL,bond->atom2,NIL,NIL,ZERO,bond->geom2,1,"",FALSE);
 	return(OK);
	}

LOCAL	bond_p	copy_bond(bond,m)
/**********************************************************
  Copy of a bond. 
***********************************************************/
bond_p		bond;
molecule_p	m;
{
        t_pointer	TAS; 
	atom_p		a1 = NIL,a2 = NIL;
	bond_p		b;

        if (bond->atom1) a1 = (atom_p)datw_inq_twin(bond->atom1->ID,NIL,NIL);
	if (bond->atom2) a2 = (atom_p)datw_inq_twin(bond->atom2->ID,NIL,NIL);
        TAS = (t_pointer)(m->HEAP);
	b = daal_create_bond(TAS,a1,a2,NIL,NIL,bond->geom1,bond->geom2,bond->order,bond->comment,bond->marq);
	return((b == (bond_p)ERROR) ? (bond_p)NIL : b);
        }

LOCAL	bond_p	copy_bond_large(bond,m)
/**********************************************************
  Copy of a bond. 
  Same as above + set b->marq to FALSE.
***********************************************************/
bond_p		bond;
molecule_p	m;
{
        t_pointer	TAS; 
	atom_p		a1 = NIL,a2 = NIL;
	bond_p		b;

        if (bond->atom1) a1 = (atom_p)datw_inq_twin(bond->atom1->ID,NIL,NIL);
	if (bond->atom2) a2 = (atom_p)datw_inq_twin(bond->atom2->ID,NIL,NIL);
        TAS = (t_pointer)(m->HEAP);
	b = daal_create_bond_large(TAS,a1,a2,NIL,NIL,bond->geom1,bond->geom2,bond->order,bond->comment,FALSE);
	return((b == (bond_p)ERROR) ? (bond_p)NIL : b);
        }

DEFINE	atom_p	daal_create_atom(TAS,name,ID,group, SB, up, down,
                                 geom,potential_type,charge,degre,
                                 comment,marq)
/**********************************************************
  Creation of an atom in the heap memory TAS. The atom
  is add to the set of atom of its group.
***********************************************************/
t_pointer	TAS; 
t_name		name;
t_integer	ID;
group_p		group;
set_bond_p	SB;
atom_p		up,down;
t_coord		geom;  /* geometry */
t_name		potential_type;  /* chemistry */
t_real		charge;
t_integer	degre;
t_bufstring	comment;
t_bool		marq;     
{
       	atom_p	atom;

       	if (TAS  == ((t_pointer)ERROR)) return((atom_p)ERROR);
        if (group == NIL) return((atom_p)ERROR);
       	atom = (atom_p)hea_alloc(TAS,sizeof(atom_t));
       	strcpy(atom->name,name); strcpy(atom->comment,comment); 
       	strcpy(atom->potential_type,potential_type);
       	atom->ID = ID; atom->geom = geom; atom->group = (t_pointer)group;
       	atom->charge = charge; atom->degre = degre; atom->SB = (t_pointer)SB;
	atom->up = up; atom->down = down;
	group->SA = (set_atom_p)heat_empile(atom,group->SA,TAS);
        atom->marq = marq;
       	return(atom);
        }

DEFINE	atom_p	daal_create_atom_large(TAS,name,ID,group, SB, up, down,
                                 geom,potential_type,charge,degre,
                                 comment,marq)
/**********************************************************
  Creation of an atom in the heap memory TAS. The atom
  is add to the set of atom of its group.
***********************************************************/
t_pointer	TAS; 
t_name		name;
t_integer	ID;
group_p		group;
set_bond_p	SB;
atom_p		up,down;
t_coord		geom;  /* geometry */
t_name		potential_type;  /* chemistry */
t_real		charge;
t_integer	degre;
t_bufstring	comment;
t_bool		marq;     
{
       	atom_p	atom;

       	if (TAS  == ((t_pointer)ERROR)) return((atom_p)ERROR);
        if (group == NIL) return((atom_p)ERROR);
        if ((atom = (atom_p)malloc(sizeof(atom_t))) == NIL) {
	   xerr_c("malloc no more space");
	   }

       	strcpy(atom->name,name); strcpy(atom->comment,comment); 
       	strcpy(atom->potential_type,potential_type);
       	atom->ID = ID; atom->geom = geom; atom->group = (t_pointer)group;
       	atom->charge = charge; atom->degre = degre; atom->SB = (t_pointer)SB;
	atom->up = up; atom->down = down;
	group->SA = (set_atom_p)heat_empile(atom,group->SA,TAS);
        atom->marq = marq;
       	return(atom);
        }

LOCAL	t_err	delete_atom_bond(bond,atom)
/**********************************************************
  Delete atom in the bond bond.
***********************************************************/
bond_p	bond;
atom_p	atom;
{
	if (bond->atom1 == atom)  bond->atom1 = NIL;
	if (bond->atom2 == atom)  bond->atom2 = NIL;
	return(OK);
	}

DEFINE	t_err	daal_delete_atom(atom)
/**********************************************************
  Delete an atom.
***********************************************************/
atom_p	atom;
{
EXTERN	t_err		daal_delete_group();
	group_p		g;
	set_atom_p	sas,sap;
	if (atom == NIL) return(ERROR);

	/* delete in the SA list */
        g = (group_p)atom->group;
	sap = g->SA; sas = sap;
        while (sas) {
              if (sas->atom == atom) break;
              sap = sas; sas = sas->succ;
              }
	if (sas == NIL) return(ERROR);
        sap->succ = sas->succ;
        if (sas == sap) g->SA = sap->succ;
	
	/* if it was the only atom of the group */
	if (g->SA == NIL) daal_delete_group(g);

	/* in the bonds */
        return(dast_all_bond(NIL,NIL,atom,delete_atom_bond,FALSE,atom));
	}

LOCAL	t_bool	connectivity = FALSE;
LOCAL	atom_p	copy_atom(atom,g)
/**********************************************************
  Creation of an atom in the heap memory TAS. The atom
  is add to the set of atom of its group.
***********************************************************/
atom_p	atom;
group_p	g;
{
        t_pointer	TAS = (t_pointer)(((molecule_p)(g->molecule))->HEAP); 
	atom_p		a;

	if (connectivity) if (atom->marq == FALSE) return(NIL);
	a_ID = atom->ID;
	a = daal_create_atom(TAS,atom->name,a_ID,g,NIL,atom->up,atom->down,
                             atom->geom,atom->potential_type,atom->charge,atom->degre,
                             atom->comment,atom->marq);
	if ((ATOM) && (atom == ATOM)) A = a; /* for dast_copy_connectivity */
	datw_set_twin(a_ID,a);
       	return(a);
        }

LOCAL	atom_p	copy_atom_large(atom,g)
/**********************************************************
  Creation of an atom in the heap memory TAS. The atom
  is add to the set of atom of its group.
  Same as above + Set m->size and a->ID.
  Created 09/20/96 (large)
***********************************************************/
atom_p	atom;
group_p	g;
{
        t_pointer	TAS = (t_pointer)(((molecule_p)(g->molecule))->HEAP); 
	atom_p		a;
	molecule_p	m = (molecule_p)g->molecule;

	if (connectivity) if (atom->marq == FALSE) return(NIL);
	++m->size; a_ID = m->size; 
	a = daal_create_atom_large(TAS,atom->name,a_ID,g,NIL,atom->up,atom->down,
                             atom->geom,atom->potential_type,atom->charge,atom->degre,
                             atom->comment,FALSE);
	if ((ATOM) && (atom == ATOM)) A = a; /* for dast_copy_connectivity */
	datw_set_twin(atom->ID,a);
       	return(a);
        }

DEFINE	group_p	daal_create_group_large(TAS,name,ID,molecule,
                                  degre,SA,SB,comment,marq,work)
/**********************************************************
  Creation of a group by malloc!
  The group is add to the set of group of its molecule.
***********************************************************/
t_pointer	TAS; 
t_name		name;
t_integer	ID;
molecule_p	molecule;
t_integer	degre;
set_atom_p	SA;
set_bond_p	SB;
t_bufstring	comment;
t_bool		marq;
t_pointer	work;
{
       	group_p	group;

       	if (TAS  == ((t_pointer)ERROR)) return((group_p)ERROR);
	if (molecule == NIL) return((group_p)ERROR);
	if ( (group = (group_p)malloc(sizeof(group_t))) == NIL) {
	   xerr_c("malloc no more space");
	   }
       	strcpy(group->name,name); strcpy(group->comment,comment); 
       	group->ID = ID; group->molecule = (t_pointer)molecule; 
	group->SA = SA; group->SB = SB; group->degre = degre;
        group->marq = marq; group->work = work;
	molecule->SG = (set_group_p)heat_empile(group,molecule->SG,TAS);
	return(group);
        }

DEFINE	group_p	daal_create_group(TAS,name,ID,molecule,
                                  degre,SA,SB,comment,marq,work)
/**********************************************************
  Creation of a group in the heap memory TAS.
  The group is add to the set of group of its molecule.
***********************************************************/
t_pointer	TAS; 
t_name		name;
t_integer	ID;
molecule_p	molecule;
t_integer	degre;
set_atom_p	SA;
set_bond_p	SB;
t_bufstring	comment;
t_bool		marq;
t_pointer	work;
{
       	group_p	group;

       	if (TAS  == ((t_pointer)ERROR)) return((group_p)ERROR);
	if (molecule == NIL) return((group_p)ERROR);
       	group = (group_p)hea_alloc(TAS,sizeof(group_t));
       	strcpy(group->name,name); strcpy(group->comment,comment); 
       	group->ID = ID; group->molecule = (t_pointer)molecule; 
	group->SA = SA; group->SB = SB; group->degre = degre;
        group->marq = marq; group->work = work;
	molecule->SG = (set_group_p)heat_empile(group,molecule->SG,TAS);
	return(group);
        }

DEFINE	t_err	daal_delete_group(group)
/**********************************************************
  Delete a group.
***********************************************************/
group_p	group;
{
	molecule_p	m;
	set_group_p	sgs,sgp;

	if (group == NIL) return(ERROR);

	/* delete in the SG list */
        m = (molecule_p)group->molecule;
	sgp = m->SG; sgs = sgp;
        while (sgs) {
              if (sgs->group == group) break;
              sgp = sgs; sgs = sgs->succ;
              }
	if (sgs == NIL) return(ERROR);
        sgp->succ = sgs->succ;
        if (sgs == sgp) m->SG = sgp->succ;

	/* the atoms */
        return(dast_all_atom(NIL,group,daal_delete_atom,FALSE,NIL));
	}


LOCAL	group_p	copy_group(group,m)
/**********************************************************
  Creation of a group in the heap memory TAS.
  The group is add to the set of group of its molecule.
***********************************************************/
group_p		group;
molecule_p	m;
{
       	group_p	g;

	connectivity = FALSE;
        g = daal_create_group(m->HEAP,group->name,group->ID,m,0,NIL,NIL,
                              group->comment,group->marq,group->work);
	dast_all_atom(NIL,group,copy_atom,FALSE,g);
	return(g);
        }

DEFINE	group_p	daal_copy_group(group,m)
/**********************************************************
  The group is copied and add to m.
  The ID of atoms and groups are recompute.
***********************************************************/
group_p		group;
molecule_p	m;
{
EXTERN	t_err		dast_set_ID_atom();
EXTERN	t_err		dast_set_ID_group();
       	group_p		g;
	t_integer	ID;

	datw_init_twin(m->size);
	g = copy_group(group,m);
	dast_all_bond(NIL,group,NIL,copy_bond,FALSE,m);
	ID = 0; dast_all_atom(m,NIL,dast_set_ID_atom,FALSE,&ID);
        m->size = ID;
        ID = 1; dast_all_group(m,dast_set_ID_group,FALSE,&ID);   
	return(g);
        }

LOCAL	molecule_p	copy_connectivity(a,m,h)
/**********************************************************
  Creation of a group in the heap memory TAS.
  The group is add to the set of group of its molecule.
***********************************************************/
atom_p		a;
molecule_p	m;
t_integer	h;
{
       	group_p		g,group;

	group = (group_p)a->group;
	if (group->marq) return(m); group->marq = TRUE;
        g = daal_create_group(m->HEAP,group->name,group->ID,m,0,NIL,NIL,
                              group->comment,FALSE,group->work);
	dast_all_atom(NIL,group,copy_atom,FALSE,g);
	return(m);
        }

LOCAL	molecule_p	copy_connectivity_large(a,m,h)
/**********************************************************
  Creation of a group in the heap memory TAS.
  The group is add to the set of group of its molecule.
  Created 09/20/96 (large)
***********************************************************/
atom_p		a;
molecule_p	m;
t_integer	h;
{
       	group_p		g,group;

	group = (group_p)a->group;
	if (group->marq) return(m); group->marq = TRUE;
        g = daal_create_group_large(m->HEAP,group->name,group->ID,m,0,NIL,NIL,
                              group->comment,FALSE,group->work);
	dast_all_atom(NIL,group,copy_atom_large,FALSE,g);
	return(m);
        }


DEFINE	atom_p	daal_copy_connectivity(atom,m)
/**********************************************************
  Copy all the atoms connected to atom and add these atoms 
  to m. 
  Return the copy of atom.
***********************************************************/
molecule_p	m;
atom_p		atom;
{
EXTERN	t_err		dast_set_ID_atom();
EXTERN	t_err		dast_set_ID_group();
EXTERN	t_err		dast_set_name_atom();
       	molecule_p	molecule;
	t_integer	ID;

	ATOM = atom; connectivity = TRUE;
	molecule = (molecule_p)((group_p)atom->group)->molecule;
	datw_init_twin(molecule->size);
	dast_set_marq(molecule,NIL,NIL,FALSE);
	dast_connectivity_atom(atom,TRUE,NIL,NIL,m,MAXINTEGER);
	dast_all_atom(molecule,NIL,copy_connectivity,TRUE,m);
	dast_all_bond(molecule,NIL,NIL,copy_bond,TRUE,m);
	dast_set_marq(molecule,NIL,NIL,FALSE);
	connectivity = FALSE;
	ID = 0; dast_all_atom(m,NIL,dast_set_ID_atom,FALSE,&ID);
        m->size = ID;
        ID = 1; dast_all_group(m,dast_set_ID_group,FALSE,&ID);
	dast_set_marq(m,NIL,NIL,FALSE);
       	return(A);
        }


DEFINE	atom_p	daal_copy_connectivity_large(atom,m)
/**********************************************************
  Copy all the atoms connected to atom and add these atoms 
  to m. 
  Return the copy of atom.
  Created 09/20/96 (large)
***********************************************************/
molecule_p	m;
atom_p		atom;
{
EXTERN	t_err		dast_set_ID_atom();
EXTERN	t_err		dast_set_ID_group();
EXTERN	t_err		dast_set_name_atom();
       	molecule_p	molecule;
	t_integer	ID;


	ATOM = atom; connectivity = TRUE;
	molecule = (molecule_p)((group_p)atom->group)->molecule;
	datw_init_twin(molecule->size);
	dast_connectivity_atom(atom,TRUE,NIL,NIL,m,MAXINTEGER);
	dast_all_atom(molecule,NIL,copy_connectivity_large,TRUE,m);
	dast_all_bond(molecule,NIL,NIL,copy_bond_large,TRUE,m);
	connectivity = FALSE;
	dast_set_marq(molecule,NIL,NIL,FALSE);
       	return(A);
        }

DEFINE	molecule_p	daal_create_molecule(TAS,name,ID,
                                             size,SG,
                                             signature,lmg,lib,
                                             comment,marq)
/**********************************************************
  Creation of a molecule in the heap memory TAS  
***********************************************************/
t_pointer	TAS; 
t_name		name;
t_integer	ID;
t_integer	size;
set_group_p	SG;
t_pointer       signature,lmg,lib;
t_bufstring	comment;
t_bool		marq;
{
       	molecule_p	molecule;
	
       	if (TAS  == ((t_pointer)ERROR)) return((molecule_p)ERROR);
       	molecule = (molecule_p)hea_alloc(TAS,sizeof(molecule_t));
        molecule->HEAP = TAS; molecule->marq = marq;
       	strcpy(molecule->name,name); strcpy(molecule->comment,comment);
        molecule->ID = ID; molecule->size = size; molecule->SG = SG;
        molecule->signature = signature; molecule->lmg = lmg; molecule->lib = lib;
		a_ID = g_ID = 0;
       	return(molecule);
        }


DEFINE	molecule_p	daal_copy_molecule(TAS,molecule)
/**********************************************************
  Creation of a molecule m from molecule  
***********************************************************/
t_pointer	TAS;
molecule_p	molecule;
{
       	molecule_p	m;

        m = daal_create_molecule(TAS,molecule->name,++m_ID,
                                 molecule->size,NIL,
                                 molecule->signature,molecule->lmg,molecule->lib,
                                 molecule->comment,molecule->marq);
	datw_init_twin(m->size);
	dast_all_group(molecule,copy_group,FALSE,m);
	dast_all_bond(molecule,NIL,NIL,copy_bond,FALSE,m);
       	return(m);
        }



