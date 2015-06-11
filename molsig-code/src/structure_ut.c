/***********************************************************
  Utility File on data structure for signature program
  Jean-loup Faulon May 91 Albuquerque
  Modified PSU Aug. - Oct. 91, PSU Jan 92, April 93.
  Start to modify the file for working with super-atom
  (field up and down).
***********************************************************/
#include <general.h>
#include <eps.h>

LOCAL   t_coord         ZERO_VECTOR = {0,0,0};

typedef struct ARG_ROOT	{
	atom_p		atom;
	t_coord		geom;
	} arg_root_t, *arg_root_p;
	 
LOCAL	t_err	root_bond(b,arg)
/**********************************************************
  no comment.
***********************************************************/
bond_p		b;
arg_root_p	arg;
{
	if ((arg->atom = b->atom1) != NIL) 
           if (arg->atom->degre > 1) { 
              arg->geom = b->geom1; return(ERROR);
	      }
	if ((arg->atom = b->atom2) != NIL) 
           if (arg->atom->degre > 1) { 
              arg->geom = b->geom2; return(ERROR);
	      }
	arg->atom = NIL; return(OK);
	}

DEFINE	atom_p	dast_root_atom(atom,geom)
/**********************************************************
  Determine the non-leaf atom connected to atom.
  A non leaf atom is an atom with a degre > 1.
***********************************************************/
atom_p		atom;
t_coord		*geom;
{
	arg_root_t	arg;

	if (atom == NIL) return(NIL);
	geom->x = geom->y = geom->z = 0;
	if (atom->degre > 1) return(atom);
	dast_all_bond(NIL,NIL,atom,root_bond,FALSE,&arg);
        *geom = arg.geom;
	return(arg.atom);
	}

DEFINE	t_err	dast_comment_atom(atom,comment)
/**********************************************************
  OK if atom->comment == comment.
***********************************************************/
atom_p		atom;
t_bufstring	comment;
{
	if (atom == NIL) return(ERROR);
	if (strcmp(atom->comment,comment) == OK) return(OK);
	return(ERROR);
	}

DEFINE	t_err	dast_geom_atom(atom,geom)
/**********************************************************
  Return the coordinates of atom. 
***********************************************************/
atom_p		atom;
p_coord		geom;
{
	if (atom == NIL) return(ERROR);
	if (geom == NIL) return(ERROR);
	*geom = atom->geom;
	return(OK);
	}

DEFINE	t_err	dast_hybridization(atom)
/**********************************************************
  Return the hybridization of atom. 
***********************************************************/
atom_p		atom;
{
	if (atom == NIL) return(ERROR);
	if (atom->potential_type[1] == 't') return(1);
	if (atom->potential_type[1] == '=') return(2);
	if (atom->potential_type[1] == 'p') return(2);
	return(3);
	}

LOCAL	t_err	bond_order(b,order)
/**********************************************************
***********************************************************/
bond_p		b;
t_integer	*order;
{
	if (b == NIL) return(ERROR);
	if (b->order > *order) *order = b->order;
	return(OK);
	}

DEFINE	t_err	dast_degree_atom(atom)
/**********************************************************
  Modify the degree of an atom. 
  The potential type of an atom may change if 
  a double or triple bond is created.
  Before modifying the atom must be saturated.
***********************************************************/
atom_p		atom;
{
	t_name		pt;
	t_integer	order;

	if (atom == NIL) return(ERROR);
	if (dast_number_bonding_site(NIL,NIL,atom) != 0) return(OK);
	atom->degre = dast_number_link_atom(atom);
	return(OK);
	}

DEFINE	t_err	dast_element_atom(atom,element)
/**********************************************************
  if element == "" return OK and the element of atom.
  Z is a dummy atom.
  ERROR if atom->potential_type != element.
	MODIFIED June 1993 by Jean-Loup Faulon
***********************************************************/
atom_p		atom;
t_name		element;
{
	t_name		e,pt;

	if (atom == NIL) return(ERROR);
	strcpy(pt,atom->potential_type);
		   if (strcmp(pt,"h_") == OK)      strcpy(e,"H");
           else if (strcmp(pt,"c_") == OK) strcpy(e,"C");
           else if (strcmp(pt,"cp") == OK) strcpy(e,"C");
           else if (strcmp(pt,"c=") == OK) strcpy(e,"C");
           else if (strcmp(pt,"ct") == OK) strcpy(e,"C");
           else if (strcmp(pt,"o_") == OK) strcpy(e,"O"); 
           else if (strcmp(pt,"o'") == OK) strcpy(e,"O"); 
           else if (strcmp(pt,"o=") == OK) strcpy(e,"O"); 
           else if (strcmp(pt,"n_") == OK) strcpy(e,"N");
           else if (strcmp(pt,"n2") == OK) strcpy(e,"N");
           else if (strcmp(pt,"np") == OK) strcpy(e,"N");
           else if (strcmp(pt,"nt") == OK) strcpy(e,"N");  
           else if (strcmp(pt,"s_") == OK) strcpy(e,"S");    
           else if (strcmp(pt,"si") == OK) strcpy(e,"Si");
           else if (strcmp(pt,"f_") == OK) strcpy(e,"F");
           else if (strcmp(pt,"cl") == OK) strcpy(e,"Cl");
           else if (strcmp(pt,"br") == OK) strcpy(e,"Br");
           else if (strcmp(pt,"i_") == OK) strcpy(e,"I");
           else if (strcmp(pt,"at") == OK) strcpy(e,"At");
           else if (strcmp(pt,"._") == OK) strcpy(e,".");
           else if (strcmp(pt,"*_") == OK) strcpy(e,"*");
           else if (pt[0] == 'Z')          strcpy(e,"Z");
           else strcpy(e,"C");
	if (element[0] == '\0') strcpy(element,e);
	if (strcmp(e,element) == OK) return(OK);
	return(ERROR);
	}

DEFINE	t_real	dast_weight_element_atom(atom)
/**********************************************************
  Return the weight of the element correponding to atom
  Jean-Loup Faulon New function introduced in April 1993 
  Added radical (.) and unknown (*).
  At = Cl
  Modified 05/2012
***********************************************************/
atom_p		atom;
{
	t_name		element;
	t_real		mass;
	t_integer	i;
	
	if (atom == NIL) return(0.0);
	strcpy(element,atom->potential_type);
	for (i = 0; element[i]; i++)
		if ((element[i] < 'A') || (element[i] > 'z')) break;
	element[i] = '\0';	
	if (strcmp(element,"H") == OK) mass=1.00794;
	else if (strcmp(element,"He") == OK) mass=4.002602;
	else if (strcmp(element,"Li") == OK) mass=6.941;
	else if (strcmp(element,"Be") == OK) mass=9.012182;
	else if (strcmp(element,"B") == OK) mass=10.811;
	else if (strcmp(element,"C") == OK) mass=12.0107;
	else if (strcmp(element,"N") == OK) mass=14.0067;
	else if (strcmp(element,"O") == OK) mass=15.9994;
	else if (strcmp(element,"F") == OK) mass=18.9984032;
	else if (strcmp(element,"Ne") == OK) mass=20.1797;
	else if (strcmp(element,"Na") == OK) mass=22.98976928;
	else if (strcmp(element,"Mg") == OK) mass=24.3050;
	else if (strcmp(element,"Al") == OK) mass=26.9815386;
	else if (strcmp(element,"Si") == OK) mass=28.0855;
	else if (strcmp(element,"P") == OK) mass=30.973762;
	else if (strcmp(element,"S") == OK) mass=32.065;
	else if (strcmp(element,"Cl") == OK) mass=35.453;
	else if (strcmp(element,"Ar") == OK) mass=39.948;
	else if (strcmp(element,"K") == OK) mass=39.0983;
	else if (strcmp(element,"Ca") == OK) mass=40.078;
	else if (strcmp(element,"Sc") == OK) mass=44.955912;
	else if (strcmp(element,"Ti") == OK) mass=47.867;
	else if (strcmp(element,"V") == OK) mass=50.9415;
	else if (strcmp(element,"Cr") == OK) mass=51.9961;
	else if (strcmp(element,"Mn") == OK) mass=54.938045;
	else if (strcmp(element,"Fe") == OK) mass=55.845;
	else if (strcmp(element,"Ni") == OK) mass=58.6934;
	else if (strcmp(element,"Co") == OK) mass=58.933195;
	else if (strcmp(element,"Cu") == OK) mass=63.546;
	else if (strcmp(element,"Zn") == OK) mass=65.38;
	else if (strcmp(element,"Ga") == OK) mass=69.723;
	else if (strcmp(element,"Ge") == OK) mass=72.64;
	else if (strcmp(element,"As") == OK) mass=74.92160;
	else if (strcmp(element,"Se") == OK) mass=78.96;
	else if (strcmp(element,"Br") == OK) mass=79.904;
	else if (strcmp(element,"Kr") == OK) mass=83.798;
	else if (strcmp(element,"Rb") == OK) mass=85.4678;
	else if (strcmp(element,"Sr") == OK) mass=87.62;
	else if (strcmp(element,"Y") == OK) mass=88.90585;
	else if (strcmp(element,"Zr") == OK) mass=91.224;
	else if (strcmp(element,"Nb") == OK) mass=92.90638;
	else if (strcmp(element,"Mo") == OK) mass=95.96;
	else if (strcmp(element,"Tc") == OK) mass=98;
	else if (strcmp(element,"Ru") == OK) mass=101.07;
	else if (strcmp(element,"Rh") == OK) mass=102.90550;
	else if (strcmp(element,"Pd") == OK) mass=106.42;
	else if (strcmp(element,"Ag") == OK) mass=107.8682;
	else if (strcmp(element,"Cd") == OK) mass=112.411;
	else if (strcmp(element,"In") == OK) mass=114.818;
	else if (strcmp(element,"Sn") == OK) mass=118.710;
	else if (strcmp(element,"Sb") == OK) mass=121.760;
	else if (strcmp(element,"Te") == OK) mass=127.60;
	else if (strcmp(element,"I") == OK) mass=126.90447;
	else if (strcmp(element,"Xe") == OK) mass=131.293;
	else if (strcmp(element,"Cs") == OK) mass=132.9054519;
	else if (strcmp(element,"Ba") == OK) mass=137.327;
	else if (strcmp(element,"La") == OK) mass=138.90547;
	else if (strcmp(element,"Ce") == OK) mass=140.116;
	else if (strcmp(element,"Pr") == OK) mass=140.90765;
	else if (strcmp(element,"Nd") == OK) mass=144.242;
	else if (strcmp(element,"Pm") == OK) mass=145;
	else if (strcmp(element,"Sm") == OK) mass=150.36;
	else if (strcmp(element,"Eu") == OK) mass=151.964;
	else if (strcmp(element,"Gd") == OK) mass=157.25;
	else if (strcmp(element,"Tb") == OK) mass=158.92535;
	else if (strcmp(element,"Dy") == OK) mass=162.500;
	else if (strcmp(element,"Ho") == OK) mass=164.93032;
	else if (strcmp(element,"Er") == OK) mass=167.259;
	else if (strcmp(element,"Tm") == OK) mass=168.93421;
	else if (strcmp(element,"Yb") == OK) mass=173.054;
	else if (strcmp(element,"Lu") == OK) mass=174.9668;
	else if (strcmp(element,"Hf") == OK) mass=178.49;
	else if (strcmp(element,"Ta") == OK) mass=180.94788;
	else if (strcmp(element,"W") == OK) mass=183.84;
	else if (strcmp(element,"Re") == OK) mass=186.207;
	else if (strcmp(element,"Os") == OK) mass=190.23;
	else if (strcmp(element,"Ir") == OK) mass=192.217;
	else if (strcmp(element,"Pt") == OK) mass=195.084;
	else if (strcmp(element,"Au") == OK) mass=196.966569;
	else if (strcmp(element,"Hg") == OK) mass=200.59;
	else if (strcmp(element,"Tl") == OK) mass=204.3833;
	else if (strcmp(element,"Pb") == OK) mass=207.2;
	else if (strcmp(element,"Bi") == OK) mass=208.98040;
	else if (strcmp(element,"Po") == OK) mass=210;
	else if (strcmp(element,"At") == OK) mass=210;
	else if (strcmp(element,"Rn") == OK) mass=222;
	else if (strcmp(element,"Fr") == OK) mass=223;
	else if (strcmp(element,"Ra") == OK) mass=226;
	else if (strcmp(element,"Ac") == OK) mass=227;
	else if (strcmp(element,"Th") == OK) mass=232.03806;
	else if (strcmp(element,"Pa") == OK) mass=231.03588;
	else if (strcmp(element,"U") == OK) mass=238.02891;
	else if (strcmp(element,"Np") == OK) mass=237;
	else if (strcmp(element,"Pu") == OK) mass=244;
	else if (strcmp(element,"Am") == OK) mass=243;
	else if (strcmp(element,"Cm") == OK) mass=247;
	else if (strcmp(element,"Bk") == OK) mass=247;
	else if (strcmp(element,"Cf") == OK) mass=251;
	else if (strcmp(element,"Es") == OK) mass=252;
	else if (strcmp(element,"Fm") == OK) mass=257;
	else if (strcmp(element,"Md") == OK) mass=258;
	else if (strcmp(element,"No") == OK) mass=259;
	else if (strcmp(element,"Lr") == OK) mass=262;
	else if (strcmp(element,"Rf") == OK) mass=261;
	else if (strcmp(element,"Db") == OK) mass=262;
	else if (strcmp(element,"Sg") == OK) mass=266;
	else if (strcmp(element,"Bh") == OK) mass=264;
	else if (strcmp(element,"Hs") == OK) mass=277;
	else if (strcmp(element,"Mt") == OK) mass=268;
	else if (strcmp(element,"Ds") == OK) mass=271;
	else if (strcmp(element,"Rg") == OK) mass=272;
	else if (strcmp(element,"Uub") == OK) mass=285;
	else if (strcmp(element,"Uut") == OK) mass=284;
	else if (strcmp(element,"Uuq") == OK) mass=289;
	else if (strcmp(element,"Uup") == OK) mass=288;
	else if (strcmp(element,"Uuh") == OK) mass=292;
	else if (strcmp(element,"Uuo") == OK) mass=294;
	else mass=300.0;
	return(mass);
}

DEFINE	t_real	dast_Z_element_atom(atom)
/**********************************************************
 Return the Z number of an atom JLF 05/2012
 ***********************************************************/
atom_p		atom;
{
	t_name		element;
	t_real		Z;
	t_integer	i;
	
	if (atom == NIL) return(0.0);
	strcpy(element,atom->potential_type);
	for (i = 0; element[i]; i++)
		if ((element[i] < 'A') || (element[i] > 'z')) break;
	element[i] = '\0';
	if (strcmp(element,"H") == OK) Z = 1;
	else if (strcmp(element,"He") == OK) Z = 2;
	else if (strcmp(element,"Li") == OK) Z = 3;
	else if (strcmp(element,"Be") == OK) Z = 4;
	else if (strcmp(element,"B") == OK) Z = 5;
	else if (strcmp(element,"C") == OK) Z = 6;
	else if (strcmp(element,"N") == OK) Z = 7;
	else if (strcmp(element,"O") == OK) Z = 8;
	else if (strcmp(element,"F") == OK) Z = 9;
	else if (strcmp(element,"Ne") == OK) Z = 10;
	else if (strcmp(element,"Na") == OK) Z = 11;
	else if (strcmp(element,"Mg") == OK) Z = 12;
	else if (strcmp(element,"Al") == OK) Z = 13;
	else if (strcmp(element,"Si") == OK) Z = 14;
	else if (strcmp(element,"P") == OK) Z = 15;
	else if (strcmp(element,"S") == OK) Z = 16;
	else if (strcmp(element,"Cl") == OK) Z = 17;
	else if (strcmp(element,"Ar") == OK) Z = 18;
	else if (strcmp(element,"K") == OK) Z = 19;
	else if (strcmp(element,"Ca") == OK) Z = 20;
	else if (strcmp(element,"Sc") == OK) Z = 21;
	else if (strcmp(element,"Ti") == OK) Z = 22;
	else if (strcmp(element,"V") == OK) Z = 23;
	else if (strcmp(element,"Cr") == OK) Z = 24;
	else if (strcmp(element,"Mn") == OK) Z = 25;
	else if (strcmp(element,"Fe") == OK) Z = 26;
	else if (strcmp(element,"Ni") == OK) Z = 27;
	else if (strcmp(element,"Co") == OK) Z = 28;
	else if (strcmp(element,"Cu") == OK) Z = 29;
	else if (strcmp(element,"Zn") == OK) Z = 30;
	else if (strcmp(element,"Ga") == OK) Z = 31;
	else if (strcmp(element,"Ge") == OK) Z = 32;
	else if (strcmp(element,"As") == OK) Z = 33;
	else if (strcmp(element,"Se") == OK) Z = 34;
	else if (strcmp(element,"Br") == OK) Z = 35;
	else if (strcmp(element,"Kr") == OK) Z = 36;
	else if (strcmp(element,"Rb") == OK) Z = 37;
	else if (strcmp(element,"Sr") == OK) Z = 38;
	else if (strcmp(element,"Y") == OK) Z = 39;
	else if (strcmp(element,"Zr") == OK) Z = 40;
	else if (strcmp(element,"Nb") == OK) Z = 41;
	else if (strcmp(element,"Mo") == OK) Z = 42;
	else if (strcmp(element,"Tc") == OK) Z = 43;
	else if (strcmp(element,"Ru") == OK) Z = 44;
	else if (strcmp(element,"Rh") == OK) Z = 45;
	else if (strcmp(element,"Pd") == OK) Z = 46;
	else if (strcmp(element,"Ag") == OK) Z = 47;
	else if (strcmp(element,"Cd") == OK) Z = 48;
	else if (strcmp(element,"In") == OK) Z = 49;
	else if (strcmp(element,"Sn") == OK) Z = 50;
	else if (strcmp(element,"Sb") == OK) Z = 51;
	else if (strcmp(element,"Te") == OK) Z = 52;
	else if (strcmp(element,"I") == OK) Z = 53;
	else if (strcmp(element,"Xe") == OK) Z = 54;
	else if (strcmp(element,"Cs") == OK) Z = 55;
	else if (strcmp(element,"Ba") == OK) Z = 56;
	else if (strcmp(element,"La") == OK) Z = 57;
	else if (strcmp(element,"Ce") == OK) Z = 58;
	else if (strcmp(element,"Pr") == OK) Z = 59;
	else if (strcmp(element,"Nd") == OK) Z = 60;
	else if (strcmp(element,"Pm") == OK) Z = 61;
	else if (strcmp(element,"Sm") == OK) Z = 62;
	else if (strcmp(element,"Eu") == OK) Z = 63;
	else if (strcmp(element,"Gd") == OK) Z = 64;
	else if (strcmp(element,"Tb") == OK) Z = 65;
	else if (strcmp(element,"Dy") == OK) Z = 66;
	else if (strcmp(element,"Ho") == OK) Z = 67;
	else if (strcmp(element,"Er") == OK) Z = 68;
	else if (strcmp(element,"Tm") == OK) Z = 69;
	else if (strcmp(element,"Yb") == OK) Z = 70;
	else if (strcmp(element,"Lu") == OK) Z = 71;
	else if (strcmp(element,"Hf") == OK) Z = 72;
	else if (strcmp(element,"Ta") == OK) Z = 73;
	else if (strcmp(element,"W") == OK) Z = 74;
	else if (strcmp(element,"Re") == OK) Z = 75;
	else if (strcmp(element,"Os") == OK) Z = 76;
	else if (strcmp(element,"Ir") == OK) Z = 77;
	else if (strcmp(element,"Pt") == OK) Z = 78;
	else if (strcmp(element,"Au") == OK) Z = 79;
	else if (strcmp(element,"Hg") == OK) Z = 80;
	else if (strcmp(element,"Tl") == OK) Z = 81;
	else if (strcmp(element,"Pb") == OK) Z = 82;
	else if (strcmp(element,"Bi") == OK) Z = 83;
	else if (strcmp(element,"Po") == OK) Z = 84;
	else if (strcmp(element,"At") == OK) Z = 85;
	else if (strcmp(element,"Rn") == OK) Z = 86;
	else if (strcmp(element,"Fr") == OK) Z = 87;
	else if (strcmp(element,"Ra") == OK) Z = 88;
	else if (strcmp(element,"Ac") == OK) Z = 89;
	else if (strcmp(element,"Th") == OK) Z = 90;
	else if (strcmp(element,"Pa") == OK) Z = 91;
	else if (strcmp(element,"U") == OK) Z = 92;
	else if (strcmp(element,"Np") == OK) Z = 93;
	else if (strcmp(element,"Pu") == OK) Z = 94;
	else if (strcmp(element,"Am") == OK) Z = 95;
	else if (strcmp(element,"Cm") == OK) Z = 96;
	else if (strcmp(element,"Bk") == OK) Z = 97;
	else if (strcmp(element,"Cf") == OK) Z = 98;
	else if (strcmp(element,"Es") == OK) Z = 99;
	else if (strcmp(element,"Fm") == OK) Z = 100;
	else if (strcmp(element,"Md") == OK) Z = 101;
	else if (strcmp(element,"No") == OK) Z = 102;
	else if (strcmp(element,"Lr") == OK) Z = 103;
	else if (strcmp(element,"Rf") == OK) Z = 104;
	else if (strcmp(element,"Db") == OK) Z = 105;
	else if (strcmp(element,"Sg") == OK) Z = 106;
	else if (strcmp(element,"Bh") == OK) Z = 107;
	else if (strcmp(element,"Hs") == OK) Z = 108;
	else if (strcmp(element,"Mt") == OK) Z = 109;
	else if (strcmp(element,"Ds") == OK) Z = 110;
	else if (strcmp(element,"Rg") == OK) Z = 111;
	else if (strcmp(element,"Uub") == OK) Z = 112;
	else if (strcmp(element,"Uut") == OK) Z = 113;
	else if (strcmp(element,"Uuq") == OK) Z = 114;
	else if (strcmp(element,"Uup") == OK) Z = 115;
	else if (strcmp(element,"Lv") == OK) Z = 116;
	else if (strcmp(element,"Uus") == OK) Z = 117;
	else if (strcmp(element,"Uuo") == OK) Z = 118;
	else Z=0.0;
	return(Z);
}


DEFINE	t_err	dast_delete_hydrogen(atom)
/******************************************************************************
 Delete hydrogen atom.
******************************************************************************/
atom_p		atom;
{
	if (dast_element_atom(atom,"H") != ERROR)
           return(daal_delete_atom(atom));
        return(OK);
	}


DEFINE	t_bool	dast_intergroup_bond(b)
/******************************************************************************
 TRUE if b is an intergroup bond.
******************************************************************************/
bond_p	b;
{
	if (b == NIL) return(FALSE);
	if (strcmp(b->comment,"intergroup_bond")) return(FALSE);
	return(TRUE);
	}

DEFINE	t_err	dast_inq_size_molecule(molecule)
/******************************************************************************
 return the size of molecule
******************************************************************************/
molecule_p	molecule;
{
 	if (molecule == NIL) return(ERROR);
	return(molecule->size);
	}

DEFINE	t_err	dast_set_size_molecule(molecule,size)
/******************************************************************************
 Fix the size of molecule
******************************************************************************/
molecule_p	molecule;
t_integer	size;
{
 	if (molecule == NIL) return(ERROR);
	molecule->size = size;
        return(OK);
	}

DEFINE	t_err	dast_set_ID_atom(atom,ID)
/******************************************************************************
 Fix the ID of an atom
******************************************************************************/
atom_p		atom;
t_integer	*ID;
{
        atom->ID = *ID; *ID += 1;
        return(OK);
	}


DEFINE	t_bool	dast_inq_degree_atom(a)
/******************************************************************************
 return the order of a bond
******************************************************************************/
atom_p		a;
{
	if (a == NIL) return(ERROR);
        return(a->degre);
	}

DEFINE	t_err	dast_set_degree_atom(a,deg)
/******************************************************************************
 Fix the order of a bond
******************************************************************************/
atom_p		a;
t_integer	deg;
{
	if (a == NIL) return(ERROR);
	a->degre = deg;
        return(OK);
	}

DEFINE	t_bool	dast_inq_order_bond(b)
/******************************************************************************
 return the order of a bond
******************************************************************************/
bond_p		b;
{
	if (b == NIL) return(ERROR);
        return(b->order);
	}

DEFINE	t_err	dast_set_order_bond(b,order)
/******************************************************************************
 Fix the order of a bond
******************************************************************************/
bond_p	b;
t_bool	order;
{
	if (b == NIL) return(ERROR);
	b->order = order;
        return(OK);
	}

DEFINE	t_err	dast_inq_potential_type_atom(atom,pt)
/******************************************************************************
 return the pt of an atom
******************************************************************************/
atom_p		atom;
t_name		pt;
{
	if ((atom == NIL) || (pt == NIL)) return(ERROR);
        strcpy(pt,atom->potential_type);
        return(OK);
	}

DEFINE	t_err	dast_set_potential_type_atom(atom,pt)
/******************************************************************************
 Fix the pt of an atom
******************************************************************************/
atom_p		atom;
t_name		pt;
{
	if (atom == NIL) return(ERROR);
        strcpy(atom->potential_type,pt);
        return(OK);
	}

DEFINE	t_err	dast_set_name_atom(atom,ID)
/******************************************************************************
 Fix the ID of an atom
******************************************************************************/
atom_p		atom;
t_integer	*ID;
{
	t_name	a_name,a_idname;

	a_name[0] = '\0', a_idname[0] = '\0';
        dast_element_atom(atom,a_name); 
        sprintf(a_idname,"%X",atom->ID+1);
	strcat(a_name,a_idname); strcpy(atom->name,a_name);
        return(OK);
	}

DEFINE	t_err	dast_set_ID_group(group,ID)
/******************************************************************************
 Fix the ID of a group
******************************************************************************/
group_p		group;
t_integer	*ID;
{
        group->ID = *ID; *ID += 1;
        return(OK);
	}

DEFINE	t_err	dast_integer_string(r,s)
/******************************************************************************
 Convert r in 4 char string ijkl
******************************************************************************/
t_integer	r;
char		s[];
{
	t_err	i = 0, j = 0, k = 0, l = 0;

	i = r/1000; 
	j = (r - 1000*i)/100;
	k = (r - 1000*i - 100*j)/10;
	l = (r - 1000*i - 100*j - 10*k)/1;
	if (i) {s[0] = i + 0x30; s[1] = j + 0x30; s[2] = k + 0x30; s[3] = l + 0x30; s[4] = '\0';}
	else if (j) {s[0] = j + 0x30; s[1] = k + 0x30; s[2] = l + 0x30; s[3] = '\0';}
	else if (k) {s[0] = k + 0x30; s[1] = l + 0x30; s[2] = '\0';}
	else {s[0] = l + 0x30; s[1] = '\0';} 
        return(OK);
	}

DEFINE	t_err	dast_string_integer(s)
/******************************************************************************
 Convert s in interger .
******************************************************************************/
char		s[];
{
	t_err	i = 0,j = 0,r = 0,n = 0,pow = 1;

	while (s[n] != '\0') n++; 
	for (i = 0; i < n; i++) {
            pow = 1;
            for (j = i; j < n-1; j++) pow = pow * 10;
	    j = s[i] - 0x30; j = j * pow ;
            r = r + j;
            }
	return(r);
	}

DEFINE t_real	dast_string_real(s)
/******************************************************************************
 Convert s in real number .
******************************************************************************/
char		s[];
{
	t_err	i = 0,j = 0,r = 0,n = 0,pow = 1,d = 1,signe;
	t_bool	flip = FALSE;

	if (s[0] == '-') signe = -1; else signe = 1;
        n = strlen(s); if (signe < 0) for (i = 0; i < n; i++) s[i] = s[i+1];
	n = strlen(s); flip = FALSE;	
        for (i = 0; i < strlen(s); i++) {
	    if (flip) d = 10 * d;
	    if (s[i] == '.') { n = n -1; flip = TRUE; }
	    if (flip) s[i] = s[i+1];
	    }

	/* the string is now transformed in integer */
	for (i = 0; i < n; i++) {
            pow = 1;
            for (j = i; j < n-1; j++) pow = pow * 10;
	    j = s[i] - 0x30; j = j * pow ;
            r = r + j;
            }
	return((t_real)(signe*(t_real)r/(t_real)d));
	}

DEFINE	t_err	dast_copy_file(file1,file2)
/******************************************************************************
 Copy file1 into file2 . (cp -r file1 file2)
*****************************************************************************/
t_bufstring	file1,file2;
{
	t_bufstring	DUMS;
	FILE		*FD1, *FD2;

	if ((FD1 = fopen(file1,"r")) == NULL) return(ERROR);   
	if ((FD2 = fopen(file2,"w")) == NULL) return(ERROR);   
        while (fscanf(FD1," %s",DUMS) != EOF) fprintf(FD2,"%s ",DUMS);
	fclose(FD1); fclose(FD2);
	return(OK);
	}


LOCAL	t_err	saturate_atom(a,pt)
/******************************************************************************
 Saturate the molecule with the atom of potential type pt.
******************************************************************************/
atom_p		a;
t_name		pt;
{
EXTERN	t_real		dara_random_real();
EXTERN	t_err		dast_print_atom();
LOCAL	t_coord	        CZERO = {0,0,0};
LOCAL	t_integer	NATM = 0;
LOCAL	group_p		group = NIL;
	t_name		name,number;
	atom_p		aa;
	bond_p		bs;
	molecule_p	molecule = NIL;
	t_coord		geom;
	t_integer	i;

	if ((a == NIL) || (pt == NIL)) return(ERROR);
	if ((group_p)a->group != group) NATM = 0;
	group = (group_p)a->group;
	molecule = (molecule_p)group->molecule;

	/* saturate the atom */
	LOOP {
	    if ((i=dast_number_bond_atom(a)) == a->degre) break;
if (i > 4) {
printf("%s %d comment = %s nb = %d\n",a->name,a->ID,a->comment,i);
dast_print_atom(a);
xerr_c("bugg in saturate atom");
}
	    dast_inq_bonding_site(a,NIL,&bs,NIL);
	    if (bs) {
	       geom = ((bs->atom1) ? bs->geom1 : bs->geom2);
               geom.x = a->geom.x + VDW * geom.x;
	       geom.y = a->geom.y + VDW * geom.y;
	       geom.z = a->geom.z + VDW * geom.z;
	       }
            else { /* adds on: all dummy atoms are located at (0,0,0) */
	       if (strcmp(pt,"*_") == OK) geom = CZERO;
	       else {
                 geom.x = a->geom.x + dara_random_real((t_real)VDW);
	         geom.y = a->geom.y + dara_random_real((t_real)VDW);
	         geom.z = a->geom.z + dara_random_real((t_real)VDW);
	         }
	       }
	    strcpy(name,pt); for (i = 0; name[i] != '\0'; i++) if (name[i] > 90) name[i] -= 32;
	    NATM++; sprintf(number,"%X",NATM); strcat(name,number);
	    aa = (atom_p)daal_create_atom(molecule->HEAP,name,a->ID + MAXATOM,a->group,NIL,
                                          NIL,NIL,geom,pt,(t_real)0,1,"",FALSE);
            daal_create_bond(molecule->HEAP,a,aa,bs,NIL,CZERO,CZERO,1,"",FALSE);
	    }
	    return(OK);
	}

DEFINE	t_err	dast_saturate_molecule(molecule,pt)
/******************************************************************************
 Saturate the molecule with the atom pt.
******************************************************************************/
molecule_p	molecule;
t_name		pt;
{
	if (pt == NIL) return(ERROR);
	dast_all_atom(molecule,NIL,saturate_atom,FALSE,pt);
	return(OK);
	}

LOCAL	t_err	delete_atom_name(atom)
/**********************************************************
  Delete atom name
***********************************************************/
atom_p	atom;
{
	if (atom == NIL) return(ERROR);
	atom->name[0] = '\0';
 	return(OK);
	}

LOCAL	t_err	delete_unsature(atom)
/**********************************************************
  Delete the atom if unsaturate
***********************************************************/
atom_p	atom;
{
	if (atom == NIL) return(ERROR);
	if (atom->name[0] == '\0') return(OK);
	if (dast_number_bonding_site(NIL,NIL,atom) > 0) {
           delete_atom_name(atom);
	   dast_connectivity_atom(atom,TRUE,NIL,delete_atom_name,NIL,MAXATOM);
	   }
	return(OK);
	}

DEFINE	t_err	dast_delete_unsature(molecule,group)
/**********************************************************
  Delete all atoms in molecule or group that are unsaturate
***********************************************************/
molecule_p	molecule;
group_p		group;
{

	if ((molecule == NIL) && (group == NIL)) return(OK);
        dast_all_atom(molecule,group,delete_unsature,FALSE,NIL);
	return(OK);
	}

typedef struct BOX {
	t_coord	   min,max;
	} box_t, *box_p;

LOCAL	t_err	boundary_atom(a,arg)
/**********************************************************
***********************************************************/
atom_p	a;
box_p	arg;
{
	if (a == NIL) return(ERROR);
	if (EPS_LT(a->geom.x,arg->min.x)) arg->min.x = a->geom.x;
	if (EPS_LT(a->geom.y,arg->min.y)) arg->min.y = a->geom.y;
	if (EPS_LT(a->geom.z,arg->min.z)) arg->min.z = a->geom.z;
	if (EPS_GT(a->geom.x,arg->max.x)) arg->max.x = a->geom.x;
	if (EPS_GT(a->geom.y,arg->max.y)) arg->max.y = a->geom.y;
	if (EPS_GT(a->geom.z,arg->max.z)) arg->max.z = a->geom.z;
	return(OK);
	}

DEFINE	t_err	dast_boundary_box(molecule,group,bond,min,max)
/**********************************************************
  Compute the boundary box of molecule or group or bond.
***********************************************************/
molecule_p	molecule;
group_p		group;
bond_p		bond;
p_coord		min, max;
{
	box_t	arg;

	if ((molecule == NIL) && (group == NIL) && (bond == NIL)) return(ERROR);
	if ((min == NIL) || (max == NIL)) return(ERROR);
	arg.min.x = INFINI; arg.min.y = INFINI; arg.min.z = INFINI; 
	arg.max.x = -INFINI; arg.max.y = -INFINI; arg.max.z = -INFINI; 
	if ((molecule) || (group)) dast_all_atom(molecule,group,boundary_atom,FALSE,&arg);
	else {
	   boundary_atom(bond->atom1,&arg);
	   boundary_atom(bond->atom2,&arg);
	   }
	*min = arg.min;
	*max = arg.max;
	return(OK);
	}

LOCAL	t_bool	PBC = FALSE, PBCX = 1, PBCY = 1, PBCZ = 1;
LOCAL	t_bool	COORD_IN_PBC = FALSE;
LOCAL	t_coord	LENGTH,ANGLE;
LOCAL	t_real	WALLX = ERROR,WALLY = ERROR,WALLZ = ERROR;
DEFINE	t_bool	dast_pbc() { return(PBC); }
DEFINE	t_err	dast_nopbc() { PBC = FALSE; return(PBC); }
DEFINE	t_err	dast_inq_periodicity(x,y,z)
/**********************************************************
  Return the periodicity in x,y and z directions
  LAMMPS convention: x = 1 mean non periodic in x direction
***********************************************************/
t_bool	*x,*y,*z;
{
	if ((x == NIL) || (y == NIL) || (z == NIL)) return(ERROR);
	*x = PBCX; *y = PBCY; *z = PBCZ;
	return(OK);
	}

DEFINE	t_err	dast_set_periodicity(x,y,z)
/**********************************************************
  Fix the periodicity in x,y and z directions
  LAMMPS convention: x = 1 mean non periodic in x direction
***********************************************************/
t_bool	x,y,z;
{
	PBCX = x; PBCY = y; PBCZ = z;
	return(OK);
	}

DEFINE	t_err	dast_inq_wall(x,y,z)
/**********************************************************
  Return the number of chain's ends attached to the walls
  x = 0, x = LENGTH.x,y = 0, y = LENGTH.y,
  and z = 0, z = LENGTH.z 
***********************************************************/
t_real	*x,*y,*z;
{
	if ((x == NIL) || (y == NIL) || (z == NIL)) return(ERROR);
	*x = WALLX; *y = WALLY; *z = WALLZ;
	return(OK);
	}

DEFINE	t_err	dast_set_wall(x,y,z)
/**********************************************************
  Fix the number of chain's ends attached to the walls
   x = 0, x = LENGTH.x,y = 0, y = LENGTH.y,
  and z = 0, z = LENGTH.z
**********************************************************/
t_real	x,y,z;
{
	WALLX = x; WALLY = y; WALLZ = z;
	return(OK);
	}

DEFINE	t_err	dast_inq_pbc(lenght,angle)
/**********************************************************
  Inquire the PBC box parameters (group is P1).
***********************************************************/
p_coord		lenght,angle;
{
	if (PBC == FALSE) return(ERROR);
	if ((lenght == NIL) || (angle == NIL)) return(ERROR);
	*lenght = LENGTH; *angle = ANGLE;
	return(OK);
	}

DEFINE	t_err	dast_set_pbc(lenght,angle)
/**********************************************************
  Set the PBC box (group is P1).
***********************************************************/
p_coord		lenght,angle;
{
	if ((lenght == NIL) || (angle == NIL)) return(ERROR);
	PBC = TRUE; LENGTH = *lenght; ANGLE = *angle;
	return(OK);
	}

DEFINE	t_bool	dast_inq_coord_in_pbc() { return(COORD_IN_PBC); }
DEFINE	t_err	dast_set_coord_in_pbc(val)
/**********************************************************
  Set the PBC box (group is P1).
***********************************************************/
t_bool	val;
{
	COORD_IN_PBC = val;
	return(OK);
	}

DEFINE	t_err	dast_coord_pbc(a,coord)
/**********************************************************
  Transfert all the coord in PBC box
***********************************************************/
atom_p		a;
p_coord		coord;
{
	t_coord		min,max,angle;

	if ((a == NIL) || (coord == NIL)) return(ERROR);
	*coord = a->geom;
	if (PBC == FALSE) return(OK);
	if (COORD_IN_PBC == FALSE) return(OK);
	min = max = ZERO_VECTOR;
	dast_inq_pbc(&max,&angle);
	while (coord->x < 0) coord->x += max.x;
	while (coord->y < 0) coord->y += max.y;
	while (coord->z < 0) coord->z += max.z;
	while (EPS_GT(coord->x,max.x)) coord->x -= max.x;
	while (EPS_GT(coord->y,max.y)) coord->y -= max.y;
	while (EPS_GT(coord->z,max.z)) coord->z -= max.z;
	if (EPS_EQ(coord->x,max.x)) coord->x = ZERO;
	if (EPS_EQ(coord->y,max.y)) coord->y = ZERO;
	if (EPS_EQ(coord->z,max.z)) coord->z = ZERO;
	if (EPS_EQ(coord->x,ZERO)) coord->x = ZERO;
	if (EPS_EQ(coord->y,ZERO)) coord->y = ZERO;
	if (EPS_EQ(coord->z,ZERO)) coord->z = ZERO;
	return(OK);
	}

DEFINE	t_err	print_group(group)
/**********************************************************
***********************************************************/
group_p	group;
{
	printf("\nGROUP NAME = %s, ID = %d, comment = %s\n",
                group->name,group->ID,group->comment);
	return(OK);
	}


DEFINE	t_err	dast_print_atom_name(atom)
/**********************************************************
***********************************************************/
atom_p	atom;
{
	if (atom == NIL) return(OK);
       	printf("ATOM = %20s, ID = %6d, degre = %6d pt = %4s\n",\
                atom->name,atom->ID,atom->degre,atom->potential_type);
        return(OK);
	}

DEFINE	t_err	dast_print_atom_debug(atom)
/**********************************************************
***********************************************************/
atom_p	atom;
{
	t_integer	n;
	set_bond_p	s;
	bond_p		b;


	if (atom == NIL) return(OK);

        n = dast_number_bond_atom(atom);


       	printf("ATOM = %s, ID = %d, degre = %d pt = %s\n",\
                atom->name,atom->ID,atom->degre,atom->potential_type);
        printf("BONDS : ");
        s = (set_bond_p)atom->SB;           
	while (s) {
	      b = s->bond; if (b == NIL) break;
              if (b->atom1) if (b->atom1 != atom) printf(" %d",b->atom1->ID);
              if (b->atom2) if (b->atom2 != atom) printf(" %d",b->atom2->ID);
	      s = s->succ;
              }
        printf("\n");

        return(OK);
	}

DEFINE	t_err	dast_print_atom(atom)
/**********************************************************
***********************************************************/
atom_p	atom;
{
	t_integer	n;
	set_bond_p	s;
	bond_p		b;


	if (atom == NIL) return(OK);
        n = dast_number_bond_atom(atom);
//	if (n == atom->degre) return(OK);
//       	printf("%s ",atom->potential_type);
 

       	printf(" ATOM = %s, ID = %d, marq = %d ",
                atom->name,atom->ID,atom->marq);
        printf(" BONDS : ");
        s = (set_bond_p)atom->SB;           
	while (s) {
	      b = s->bond; if (b == NIL) break;
              if (b->atom1) printf("%s-",b->atom1->name);
              else printf("NIL-");
              if (b->atom2) printf("%s ",b->atom2->name);
              else printf("NIL ");
	      s = s->succ;
              }
        printf("\n");

        return(OK);
	}

DEFINE	t_err	dast_print_bond(bond)
/**********************************************************
***********************************************************/
bond_p	bond;
{
	if (bond == NIL) return(OK);
	if (bond->comment == NIL) return(OK);
	if (strcmp(bond->comment,"") == OK) return(OK); 
	return(OK);
	}

DEFINE	t_err	dast_print_group(group)
/**********************************************************
 
***********************************************************/
group_p	group;
{
	print_group(group);
	dast_all_atom(NIL,group,dast_print_atom,FALSE,NIL); 
	return(OK);
	}

DEFINE	t_err	dast_print_molecule_debug(molecule)
/**********************************************************
 
***********************************************************/
molecule_p	molecule;
{
	dast_all_atom(molecule,NIL,dast_print_atom_debug,FALSE,NIL); 
	return(OK);
	}

DEFINE	t_err	dast_print_molecule(molecule)
/**********************************************************
 
***********************************************************/
molecule_p	molecule;
{
	dast_all_atom(molecule,NIL,dast_print_atom,FALSE,NIL); 
	return(OK);
	}

DEFINE	t_err	dast_print_connectivity(atom)
/**********************************************************
  Print all atoms connected to atom.
***********************************************************/
atom_p		atom;
{
	molecule_p	m;

	m = (molecule_p)((group_p)atom->group)->molecule;
	dast_set_marq(m,NIL,NIL,FALSE);
	dast_connectivity_atom(atom,TRUE,NIL,NIL,NIL,MAXATOM); 
        dast_all_atom(m,NIL,dast_print_atom,TRUE,NIL);
	dast_set_marq(m,NIL,NIL,FALSE);
	return(OK);
	}

DEFINE	t_err	dast_print_connectivity_debug(atom)
/**********************************************************
  Print all atoms connected to atom.
***********************************************************/
atom_p		atom;
{
	molecule_p	m;

	m = (molecule_p)((group_p)atom->group)->molecule;
	dast_set_marq(m,NIL,NIL,FALSE);
	dast_connectivity_atom(atom,TRUE,NIL,NIL,NIL,MAXATOM); 
        dast_all_atom(m,NIL,dast_print_atom_debug,TRUE,NIL);
	dast_set_marq(m,NIL,NIL,FALSE);
	return(OK);
	}



        
