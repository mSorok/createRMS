/***********************************************************
  File for counting in data structure of the 
  signature program
  Jean-loup Faulon May 91 Albuquerque
  Modified PSU Aug. - Oct. 91, PSU Jan 92, April 93, 
                      April 97
***********************************************************/
#include <general.h>

LOCAL	t_err	number_atom(atom,n)
/**********************************************************
 no comment
***********************************************************/
atom_p		atom;
t_integer	*n;
{
	if (atom == NIL) return(ERROR);
	*n = (*n) + 1;
	return(OK);
	}

DEFINE	t_err	dast_number_atom(molecule,group,atom)
/**********************************************************
  Compute the number of atoms in a molecule or a group.
  If atom != NIL only the atom connected to atom are counted.
***********************************************************/
molecule_p	molecule;
group_p		group;
atom_p		atom;
{
	molecule_p	m;
	t_integer	n = 0;
	t_bool		marq;

	if ((molecule == NIL) && (group == NIL)) return(ERROR);
	m = ((molecule) ? molecule : (molecule_p)group->molecule);
	if (atom) dast_connectivity_atom2(atom,NIL,number_atom,&n,MAXATOM);
	else      dast_all_atom(molecule,group,number_atom,FALSE,&n);
	return(n);
	}

EXTERN	t_real	dast_weight_element_atom();
LOCAL	t_err	weight_atom(atom,mw)
/**********************************************************
 no comment
***********************************************************/
atom_p		atom;
t_real		*mw;
{
	if (atom == NIL) return(ERROR);
	*mw = *mw + dast_weight_element_atom(atom);
	return(OK);
	}

DEFINE	t_real	dast_molecular_weight(molecule,group,atom)
/**********************************************************
  Compute the molecular weight of  a molecule a group or
  an atom. if atom != NIL then compute the mass of the
  fragment containing atom.
***********************************************************/
molecule_p	molecule;
group_p		group;
atom_p		atom;
{
	t_real		mw = 0;
	molecule_p	m;
	t_integer	n;
	t_bool		marq = FALSE;

	if ((molecule == NIL) && (group == NIL) && (atom == NIL)) return(ERROR);
	if ((molecule == NIL) && (group == NIL)) return(dast_weight_element_atom(atom));
	m = ((molecule) ? molecule : (molecule_p)group->molecule);
	if (atom) dast_connectivity_atom2(atom,NIL,weight_atom,&mw,MAXATOM);
	else      dast_all_atom(molecule,group,weight_atom,FALSE,&mw);
	return(mw);
	}

LOCAL	t_err	mass_center_atom(atom,center)
/**********************************************************
 no comment
***********************************************************/
atom_p		atom;
p_coord		center;
{
	if (atom == NIL) return(ERROR);
	center->x += atom->geom.x;
	center->y += atom->geom.y;
	center->z += atom->geom.z;
	return(OK);
	}

DEFINE	t_err	dast_molecular_mass_center(molecule,group,atom,center)
/**********************************************************
  Compute the center of mass of a molecule a group or
  an atom. if atom != NIL then compute the center of the
  fragment containing atom.
***********************************************************/
molecule_p	molecule;
group_p		group;
atom_p		atom;
p_coord		center;
{
	molecule_p	m;
	t_integer	n;
	t_bool		marq = FALSE;

	if (center == NIL) return(ERROR); center->x = center->y = center->z = 0;
	if ((molecule == NIL) && (group == NIL) && (atom == NIL)) return(ERROR);
	if ((molecule == NIL) && (group == NIL)) return(mass_center_atom(atom,center));
	m = ((molecule) ? molecule : (molecule_p)group->molecule);

	n = molecule->size;
	if (atom) n = dast_connectivity_atom2(atom,NIL,mass_center_atom,center,MAXATOM);
	else      dast_all_atom(molecule,group,mass_center_atom,FALSE,center);
	if (n < 1) return(ERROR);
	center->x = center->x/(t_real)n;
	center->y = center->y/(t_real)n;
	center->z = center->z/(t_real)n;
	return(OK);
	}

DEFINE	t_err	dast_number_bond_atom(atom)
/**********************************************************
 Jean-Loup Faulon 
 BUG (Apr 93) : In case of intergoup bond n = n + b->order.
 On can create more than one bond between two atoms, 
 and one can change the atom type c -> c= for example.
 This was modified because of the data structure and hence
 the molecular files read by Signature:
 With Biosym, MSI, Pcmodel there are no 1.5 bond types therefore
 the valences of cp is 3
                 c= is 2
                 o= is 1
 However, when creating a bond ("intergroup bond") one can change
 the type of an atom c -> c= or cp, therefore it is necessary
 to include an order for the bonds that are created.
 In short for initial bonds order = 1 always, for the bonds 
 created order can be > 1.
***********************************************************/
atom_p		atom;
{
	set_bond_p	s;
	bond_p		b;
	t_real		n = 0;

	if (atom == NIL) return(ERROR);
        s = (set_bond_p)atom->SB;	
	while (s) {
	  b = s->bond;
          if (b == (bond_p)ERROR) break;
	  if (b) if ((b->atom1) && (b->atom2)) {
	     if ((atom->potential_type[1] == '_') &&
                 (strcmp(b->comment,"intergroup_bond") == OK))
	        n = n +  b->order;
	     else n = n + 1;
             }
	  s = s->succ;
          }
	return((t_err)n);
	}

LOCAL	t_err	number_link_atom(atom,N)
/**********************************************************
 Jean-Loup Faulon 
 BUG (Apr 93) : Number of bond without order
***********************************************************/
atom_p		atom;
t_integer	*N;
{
	set_bond_p	s;
	bond_p		b;
	t_integer	n = 0;

	if ((atom == NIL) || (N == NIL)) return(ERROR);

        s = (set_bond_p)atom->SB;
	while (s) {
	  b = s->bond;
          if (b == (bond_p)ERROR) break;
	  if (b) if ((b->atom1) && (b->atom2)) n++;
	  s = s->succ;
          }
	*N += n;
	return(OK);
	}

DEFINE	t_err	dast_number_link_atom(atom)
/**********************************************************
 Jean-Loup Faulon 
 BUG (Apr 93) : Number of bond without order
***********************************************************/
atom_p		atom;
{
	set_bond_p	s;
	bond_p		b;
	t_integer	n = 0;

	if (atom == NIL) return(ERROR);
	number_link_atom(atom,&n);
	return(n);
	}

DEFINE	t_err	dast_number_link(molecule,group,atom)
/**********************************************************
  Compute the number of links (bond without order)
  in a molecule a group or an atom.
  If molecule or group and atom != NIL only the atom 
  connected to atom are counted.
***********************************************************/
molecule_p	molecule;
group_p		group;
atom_p		atom;
{
	molecule_p	m;
	t_integer	n,number = 0;
	t_bool		marq = FALSE;

	if ((molecule == NIL) && (group == NIL) && (atom == NIL))
           return(ERROR);
	if ((molecule == NIL) && (group == NIL)) 
           return(dast_number_link_atom(atom)); 
	m = ((molecule) ? molecule : (molecule_p)group->molecule);
	if (atom) dast_connectivity_atom2(atom,NIL,number_link_atom,&number,MAXATOM);
        else      dast_all_atom(molecule,group,number_link_atom,FALSE,&number);
	return(number/2);
	}

LOCAL	t_err	print_bonding_site(a,number)
/**********************************************************
 Jean-Loup Faulon
 OK April 93
***********************************************************/
atom_p		a;
t_integer	*number;
{
	t_integer	n;

	n = a->degre - dast_number_bond_atom(a);
if (n > 0) printf("deg(%s) = %d ",a->potential_type,n);
	*number += n;
	return(OK);
	}

LOCAL	t_err	number_bonding_site(a,number)
/**********************************************************
 Jean-Loup Faulon
 OK April 93
***********************************************************/
atom_p		a;
t_integer	*number;
{
	t_integer	n;

	n = a->degre - dast_number_bond_atom(a);
if (n) {
}
	*number += n;
	return(OK);
	}

LOCAL	t_err	number_H_bonding_site(a,number)
/**********************************************************
 Jean-Loup Faulon
 OK April 93
***********************************************************/
atom_p		a;
t_integer	*number;
{
EXTERN	t_err		dast_element_atom();
	t_integer	n;
	t_name		e;

	if (dast_element_atom(a,"H") != OK) return(OK);
	n = a->degre - dast_number_bond_atom(a);
	*number += n;
	return(OK);
	}

DEFINE	t_err	dast_number_bonding_site(molecule,group,atom)
/**********************************************************
  Compute the number of bonding site in a molecule a group
  or an atom.
  If molecule or group and atom != NIL only the atom 
  connected to atom are counted.
***********************************************************/
molecule_p	molecule;
group_p		group;
atom_p		atom;
{
	molecule_p	m;
	t_integer	n,number = 0;
	t_bool		marq = FALSE;

	if ((molecule == NIL) && (group == NIL) && (atom == NIL))
           return(ERROR);

	if ((molecule == NIL) && (group == NIL)) {
           number_bonding_site(atom,&number); 
           return(number);
           }
	m = ((molecule) ? molecule : (molecule_p)group->molecule);
/*
	if (atom) marq = TRUE; else marq = FALSE;
        if (atom) n = dast_connectivity_atom(atom,marq,NIL,NIL,NIL,MAXATOM);
        dast_all_atom(molecule,group,number_bonding_site,marq,&number);
	dast_set_marq(m,NIL,NIL,FALSE);
*/
	if (atom) dast_connectivity_atom2(atom,NIL,number_bonding_site,&number,MAXATOM);
        else      dast_all_atom(molecule,group,number_bonding_site,FALSE,&number);

	return(number);
	}

DEFINE	t_err	dast_print_bonding_site(molecule,group,atom)
/**********************************************************
  Print bonding site in a molecule a group
  or an atom.
  If molecule or group and atom != NIL only the atom 
  connected to atom are counted.
***********************************************************/
molecule_p	molecule;
group_p		group;
atom_p		atom;
{
	molecule_p	m;
	t_integer	n,number = 0;
	t_bool		marq = FALSE;

	if ((molecule == NIL) && (group == NIL) && (atom == NIL))
           return(ERROR);

	if ((molecule == NIL) && (group == NIL)) {
           number_bonding_site(atom,&number); 
           return(number);
           }
printf("LIST BS : ");
	m = ((molecule) ? molecule : (molecule_p)group->molecule);
	if (atom) marq = TRUE; else marq = FALSE;
        if (atom) n = dast_connectivity_atom(atom,marq,NIL,NIL,NIL,MAXATOM);
        dast_all_atom(molecule,group,print_bonding_site,marq,&number);
	dast_set_marq(m,NIL,NIL,FALSE);
printf("number = %d\n",number);
	return(number);
	}

DEFINE	t_err	dast_number_H_bonding_site(molecule,group,atom)
/**********************************************************
  Compute the number of H that are 
  bonding site in a molecule a group or an atom.
  If molecule or group and atom != NIL only the atom 
  connected to atom are counted.
***********************************************************/
molecule_p	molecule;
group_p		group;
atom_p		atom;
{
	molecule_p	m;
	t_integer	n,number = 0;
	t_bool		marq = FALSE;

	if ((molecule == NIL) && (group == NIL) && (atom == NIL))
           return(ERROR);

	if ((molecule == NIL) && (group == NIL)) {
           number_bonding_site(atom,&number); 
           return(number);
           }

	m = ((molecule) ? molecule : (molecule_p)group->molecule);
        if (atom) dast_connectivity_atom2(atom,NIL,number_H_bonding_site,&number,MAXATOM);
        else      dast_all_atom(molecule,group,number_H_bonding_site,FALSE,&number);
	return(number);
	}

LOCAL	t_err	marq_hydrogen(atom)
/**********************************************************
  no comment.
***********************************************************/
atom_p	atom;
{
EXTERN	t_err	dast_element_atom();

	if (dast_element_atom(atom,"H") == OK) atom->marq = MAXATOM;
	return(OK);
	}

LOCAL	t_err	non_connect_atom(atom,a)
/**********************************************************
  no comment.
***********************************************************/
atom_p	atom,*a;
{
	*a = NIL;
	if (atom->marq == FALSE) { *a = atom; return(ERROR); }
	return(OK);
	}
	
DEFINE	t_err	dast_connectivity(molecule,group,hydrogen)
/**********************************************************
  Determine the number of connected 
 components in a molecule
  or a group.
  If hydrogen == TRUE, hydrogen are not counted.
***********************************************************/
molecule_p	molecule;
group_p		group;
t_bool		hydrogen;
{
	molecule_p	m;
	atom_p		a;
	t_integer	n,number = 0;

	m = ((molecule) ? molecule : (molecule_p)group->molecule);
	a = m->SG->group->SA->atom;
	dast_set_marq(m,NIL,NIL,FALSE);
	if (hydrogen) dast_all_atom(m,group,marq_hydrogen,FALSE,NIL);
	LOOP {
	  dast_all_atom(m,group,non_connect_atom,FALSE,&a);
          if (a == NIL) break;
	  n = dast_connectivity_atom(a,(t_bool)(number + 1),NIL,NIL,NIL,MAXATOM); number++;
	  }
	dast_set_marq(m,NIL,NIL,FALSE);	  
	return(number);
	}
	
LOCAL	t_err	component_atom(atom,SIZEDIS)
/**********************************************************
 no comment.
 ***********************************************************/
atom_p		atom;
t_integer	SIZEDIS[];
{
	SIZEDIS[atom->ID] = atom->marq;
	return(OK);
}


DEFINE	t_err	dast_connectivity_distribution(molecule,group,hydrogen,SIZEDIS)
/**********************************************************
 Determine the number of connected 
 components in a molecule and return the distribution of connected 
 component sizes.
 If hydrogen == TRUE, hydrogen are not counted.
 JLF 02/2012
 ***********************************************************/
molecule_p	molecule;
group_p		group;
t_bool		hydrogen;
t_integer	SIZEDIS[];
{
	molecule_p	m;
	atom_p		a;
	t_integer	n,number = 0,i;
	
	if (SIZEDIS == NIL) return(ERROR);
	m = ((molecule) ? molecule : (molecule_p)group->molecule);
	for (i = 0; i < m->size; i++) SIZEDIS[i] = 0;
	a = m->SG->group->SA->atom;
	dast_set_marq(m,NIL,NIL,FALSE);
	if (hydrogen) dast_all_atom(m,group,marq_hydrogen,FALSE,NIL);
	LOOP {
		dast_all_atom(m,group,non_connect_atom,FALSE,&a);
		if (a == NIL) break;
		n = dast_connectivity_atom(a,(t_bool)(number + 1),NIL,NIL,NIL,MAXATOM); number++;
// dast_print_molecule(molecule);
	}
	dast_all_atom(m,group,component_atom,FALSE,SIZEDIS);
	dast_set_marq(m,NIL,NIL,FALSE);	  
	return(number);
}

DEFINE	t_err	dast_connectivity_function(molecule,group,
		                f1,farg1,f2,farg2,hydrogen)
/**********************************************************
  Apply the function f1 to all connex component of molecule or 
  group. The function return the number of components.
  If hydrogen == TRUE, hydrogen are not counted.
  and the function f is not applied to hydrogen.
***********************************************************/
molecule_p	molecule;
group_p		group;
f_traitement	f1,f2;
t_pointer	farg1,farg2;
t_bool		hydrogen;
{
	molecule_p	m;
	atom_p		a;
	t_integer	n,number = 0;

	m = ((molecule) ? molecule : (molecule_p)group->molecule);
	a = m->SG->group->SA->atom;
	dast_set_marq(m,NIL,NIL,FALSE);
	if (hydrogen) dast_all_atom(m,group,marq_hydrogen,FALSE,NIL);
	LOOP {
	  dast_all_atom(m,group,non_connect_atom,FALSE,&a);
	  if (a == NIL) break;
	  n = dast_connectivity_atom(a,(t_bool)(number + 1),NIL,f1,farg1,MAXATOM); number++;
	  if (f2) f2(farg1,farg2);
	  }
	dast_set_marq(m,NIL,NIL,FALSE);	  
	return(number);
	}

        
