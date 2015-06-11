/***********************************************************
  Computation and manipulation of the signature.
  Create Jean-loup Faulon Aug. 91 
  Modified Pennstate dec. 91. Jan 92.
  The signature is composed of real number (ratio for
  100 C). To compute integer number the number of
  carbon is necessary.
  MODIFIED JULY 92
  INTRODUCED THE NOTION OF MOLECULAR GROUPS COMPOSED
  EVENTUALY OF SEVERAL GROUPS. ONE CAN COMPUTE THE
  NUMBER OF BONDING SITES IN A MOLECULAR GROUP.
  Modified Sept. 93 (if no carbon max and min carbon
  is max and min number of atoms).
***********************************************************/
#include <general.h>
#include <eps.h>
#include "signature.h"
#define MAXSIGSTRING	10000

EXTERN	t_pointer	hea_alloc(),heat_creer_hea();
EXTERN	t_pointer	sial_create_atom_signature(),sial_create_signature();

typedef	struct	ARG {
		t_pointer	HEAP;
		group_p		group;
                atom_p		root;
		signature_p	signature;
                t_integer	number;
		} arg_t, *arg_p;

LOCAL	t_integer	COLOR,MINC,MAXC,MINXL,MAXXL;
LOCAL	t_integer	DIMENSION,DIMENSION_MIN = -1, DIMENSION_MAX = -1;
LOCAL	t_integer	MAXDIM	= 5;

LOCAL	t_err	dimension_size(n)
/**********************************************************
 Give the size of the string for a signature
 of dimension n.
***********************************************************/
t_integer	n;
{
	t_integer	s,p3,i;

	if (n == 0) return(1);
	s = 1; p3 = 1;
	for (i = 1; i < n; i++) {
            p3 = p3 * 3; s = s + p3;	
	}
        return(4*s+1);
	}

LOCAL	t_err	size_dimension(s)
/**********************************************************
 reverse of the previous function.
***********************************************************/
p_string	s;
{
	t_integer	l,n;

	if (s == NIL) return(ERROR);
	l = strlen(s);
	for (n = 1; n <= MAXDIM; n++) 
	  if (l < dimension_size(n)) return(n-1);
	return(MAXDIM);
	}

DEFINE	t_err  sisi_set_dimension(d)
/**********************************************************
 no comment.
***********************************************************/
t_integer	d;
{
	DIMENSION = d;
	return(OK);
	}
	
DEFINE	t_err  sisi_set_dimension_min(d)
/**********************************************************
 no comment.
***********************************************************/
t_integer	d;
{
	DIMENSION_MIN = d;
	return(OK);
	}
DEFINE	t_err  sisi_set_dimension_max(d)
/**********************************************************
 no comment.
***********************************************************/
t_integer	d;
{
	DIMENSION_MAX = d;
	return(OK);
	}


DEFINE	t_err  sisi_inq_dimension(s)
/**********************************************************
 no comment.
***********************************************************/
p_string	s;
{
	if (s == NIL) return(DIMENSION);
	return(size_dimension(s));
	}
	
DEFINE	t_err  sisi_inq_dimension_min()
/**********************************************************
 no comment.
***********************************************************/
{
	if (DIMENSION_MIN < 0) return(DIMENSION);
	return(DIMENSION_MIN);
	}
	
DEFINE	t_err  sisi_inq_dimension_max()
/**********************************************************
 no comment.
***********************************************************/
{
	if (DIMENSION_MAX < 0) return(DIMENSION);
	return(DIMENSION_MAX);
	}
	

DEFINE	t_err  sisi_set_color(d)
/**********************************************************
 no comment.
***********************************************************/
t_integer	d;
{
	COLOR = d;
	return(OK);
	}


DEFINE	t_err  sisi_inq_color()
/**********************************************************
 no comment.
***********************************************************/
{
	return(COLOR);
	}

DEFINE	t_err  sisi_set_carbon(min,max)
/**********************************************************
 no comment. 
***********************************************************/
t_integer	min,max;
{
	MINC = min; MAXC = max;
	if (MINC <= 0) MINC = 0;
	if (MAXC <= 0) MAXC = MAXINTEGER;
	return(OK);
	}

DEFINE	t_err  sisi_inq_min_carbon()
/**********************************************************
 no comment.
***********************************************************/
{	return(MINC);
	}

DEFINE	t_err  sisi_inq_max_carbon()
/**********************************************************
 no comment.
***********************************************************/
{	return(MAXC);
	}


typedef	struct SAST {
	t_pointer	O,E; /* Origin, Extremity */
	t_err		R; /* Relation */
	} SAS_t, *SAS_p; /* LOCAL Data struct for relation signature */

LOCAL	t_pointer	allocate_SAS(HEAP,size)
/**********************************************************
	LOCAL SAS data struct allocation
	CREATED Jean-Loup Faulon, May 1993            
***********************************************************/
t_pointer	HEAP;
t_integer	size;
{
	if (HEAP == NIL) return(NIL);
	return((t_pointer)hea_alloc(HEAP,size * sizeof(SAS_t)));
	}

LOCAL	t_err	inq_SAS(F,SAS,i)
/**********************************************************
	LOCAL SAS data struct manipulation
	CREATED Jean-Loup Faulon, May 1993            
***********************************************************/
t_bufstring     F; /* Field = Origin, Extermity, Relation */
SAS_t		SAS[];
t_err		i;
{
	if (i == ERROR) return(ERROR);
	if (SAS == NIL) return(ERROR);
	if (strcmp(F,"O") == OK) return((t_err)SAS[i].O);
	if (strcmp(F,"E") == OK) return((t_err)SAS[i].E);
	if (strcmp(F,"R") == OK) return((t_err)SAS[i].R);
	return(ERROR);
	}

LOCAL	t_err	set_SAS(SAS,i,O,E,R)
/**********************************************************
	LOCAL SAS data struct manipulation
	CREATED Jean-Loup Faulon, May 1993            
***********************************************************/
SAS_t		SAS[];
t_err		i,R;
t_pointer	O,E;
{
	if (i == ERROR) return(ERROR);
	if (SAS == NIL) return(ERROR);
	SAS[i].O = O; SAS[i].E = E; SAS[i].R = R;
	return(OK);
	}



DEFINE	t_integer sisi_relation_signature(s1,s2)
/**********************************************************
 Relation between s1 and s2 : NE GT LT EQ.
 The lower relation is returned.
 s1 is in the list s2 is in the graph.
 For each element of s1 (the list) an element of s2 
 is searched.
 
 There was a bug in this function
 imagine : as1 = 2cc - 2c 
           and as2 = cccpcph + ccchh - ccpcph - cchh
 as1 <= as2 but v1 > v2 the function answered NE
 MODIFICATION by Jean-loup Faulon Sept. 92 :
     Do not take acount of v1 and v2 for bonds.
 MODIFIED in May 1993. Introduction of SAS.           
***********************************************************/
signature_p	s1,s2;
{
LOCAL	t_pointer		hea_SAS = NIL;
	t_integer		i,i1 = 0,i2 = 0, i12 = ERROR;
	t_real			v1,v2;
	atom_signature_p	as1,as2;
	t_err			r,r12;
	t_pointer		SAS1,SAS2;
	t_integer		NAS1 = 0,NAS2 = 0;

	/* initialize heap_SAS */
	if (hea_SAS == NIL)
	if ((hea_SAS = (t_pointer)heat_creer_hea(MAXHEAP)) == NIL) return(ERROR);
	hea_raz(hea_SAS);

	if ((s1 == NIL) || (s2 == NIL)) return(NE);

	/* allocate SAS1 calculate size of s1 */
        i1 = NAS1 = 0;
	LOOP {
	  siut_inq_signature(s1,i1,&v1,&as1);
          if ((as1 == NIL) || (v1 < 0)) break;
	  NAS1 += (t_integer)v1;
          i1++;
	  }
	SAS1 = allocate_SAS(hea_SAS,NAS1);

	/* allocate SAS2 calculate size of s2 */
        i2 = NAS2 = 0;
	LOOP {
	  siut_inq_signature(s2,i2,&v2,&as2);
          if ((as2 == NIL) || (v2 < 0)) break;
	  NAS2 += (t_integer)v2;
          i2++;
	  }
	SAS2 = allocate_SAS(hea_SAS,NAS2);

        /* load SAS1, only positive value */
        i1 = NAS1 = 0;
	LOOP {
	  siut_inq_signature(s1,i1,&v1,&as1);
          if (as1 == NIL) break; 
	  for (i = 0; i < (t_integer)v1; i++) set_SAS(SAS1,NAS1++,as1,NIL,NE);
          i1++;
	  }
	set_SAS(SAS1,NAS1,NIL,NIL,NE);

        /* load SAS2, only positive value */
        i2 = NAS2 = 0;
	LOOP {
	  siut_inq_signature(s2,i2,&v2,&as2);
          if (as2 == NIL) break; 
	  for (i = 0; i < (t_integer)v2; i++) set_SAS(SAS2,NAS2++,as2,NIL,NE);
          i2++;
	  }
	set_SAS(SAS2,NAS2,NIL,NIL,NE);

        /* SAS1[i].R = MAX (relation(SAS1[i],SAS2[]) */
	for (i1 = 0; i1 < NAS1; i1++) {
          r12 = NE; i12 = ERROR;
	  for (i2 = 0; i2 < NAS2; i2++) {
            if (inq_SAS("E",SAS2,i2)) continue; /* only one association */
            r = siut_relation_atom_signature(inq_SAS("O",SAS1,i1),inq_SAS("O",SAS2,i2));
            if (r > r12) { 
               i12 = i2; r12 = r; set_SAS(SAS1,i1,inq_SAS("O",SAS1,i1),inq_SAS("O",SAS2,i2),r12);
	       }
            }

	  /* mark SAS2[i12] */
          if (i12 != ERROR) set_SAS(SAS2,i12,inq_SAS("O",SAS2,i12),inq_SAS("O",SAS1,i1),r12);
          }

	/* r = MIN (SAS1[].R) */
	r = EQ;
	for (i1 = 0; i1 < NAS1; i1++)
          if (inq_SAS("R",SAS1,i1) < r) r = inq_SAS("R",SAS1,i1);
	if ((r == EQ) && (NAS1 < NAS2)) r = LT;
	if ((r == EQ) && (NAS1 > NAS2)) r = GT;
	return(r);
	}


DEFINE	t_bool sisi_null_signature(signature,dim,as,value)
/**********************************************************
 TRUE if all values are equal to 0.
 Only the atom signature where as->dim <= dim are tested.
***********************************************************/
signature_p		signature;
t_integer		dim;
atom_signature_p	*as;
t_real			*value;
{
	t_integer		i = 0;
	t_real			v,DELTA = 0;

	if (signature == NIL) return(TRUE);

	i = 0; LOOP {
	  siut_inq_signature(signature,i,value,as);
	  if ((*as) == NIL) break;
	  if ((*as)->dim > dim) { i++; continue; }
          DELTA = ((*as)->max - (*as)->min)/(double)2;
          v = ((*value > 0) ? *value : -*value);
          if (EPS_GT(v,DELTA)) {
             return(FALSE);
             }
	  i++;
	  }
	return(TRUE);
	}

DEFINE	t_bool sisi_positive_signature(signature,dim,as,value)
/**********************************************************
 TRUE if all values are >= 0.
 Only the atom signature where as->dim <= dim are tested.
***********************************************************/
signature_p		signature;
t_integer		dim;
atom_signature_p	*as;
t_real			*value;
{
	t_integer		i = 0;
	t_real			v,DELTA = 0;

	if (signature == NIL) return(TRUE);

	i = 0; LOOP {
	  siut_inq_signature(signature,i,value,as);
	  if ((*as) == NIL) break;
	  if ((*as)->dim > dim) { i++; continue; }
          DELTA = ((*as)->max - (*as)->min)/(double)2;
          if (EPS_LT(*value,ZERO)) if (EPS_GT(-*value,DELTA)) {
             return(FALSE);
             }
	  i++;
	  }
	return(TRUE);
	}


DEFINE	signature_p sisi_sub_bond_signature(HEAP,s1,k,s2,relation)
/**********************************************************
 s1 = s1 + k * s2.
 if relation(s1->as,s2->as) >= relation 
 s1->value += k * s2->value.

 WARNING : THIS FUNCTION IS NOT PROPER. IT WAS WRITTEN TO
	   HANDLE THE RESOLUTION OF THE SIGNATURE EQUATION
	   WHEN DEALING WITH INTERFRAGMENTS BONDS DEFINED
	   ONLY WITH TWO ATOMS (DIMENSION 0)
	   IN THIS CASE WE NEED TO SUBSTRACT THE BOND
           ONLY FROM THE REMAINING NON-ZERO PARAMETER
           OF THE SIGNATURE.
	   THERE IS STILL A BUG IF TWO PARAMETER ARE
  	   NON ZERO AND THE BOND SHOULD ONLY INFLUENCE 1;
           BOTH OF THEM WILL BE MODIFIED.
	   FOR EXAMPLE 2 cpcp - 2 cp WILL CHANGE NON-ZERO
           SIGNATURE FOR :
	             cp(cpcpcp), cp(cpcpo), cp(cpcpc),
                     cp(cpcph), cp(cpcpn) and so on...

	   THE REAL SOLUTION IS TO MAKE A LIBRARY
           OF INTERFRAGMENT BONDS AND KEEP IN MIND THAT
           A BOND HAS TO INFLUENCE ONLY ONE PARAMETER
           OF THE SIGNATURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   THEREFORE THE FUNCTION BELOW IS COOCKING RECIPE AND
           MUST BE REMOVED.
***********************************************************/
t_pointer	HEAP;
signature_p	s1,s2;
t_real		k;
t_err		relation;
{
	t_integer		i1 = 0,i2 = 0;
	t_real			v1,v2;
	atom_signature_p	as1,as2;
	t_bool			add = FALSE;
	t_err			r;

	LOOP {
	  siut_inq_signature(s2,i2,&v2,&as2); if (as2 == NIL) break; 
	  if (HEAP == NIL) if (EPS_LT(v2,ZERO)) { i2++; continue; }
	  i1 = 0; add = FALSE;
	  LOOP {
	    siut_inq_signature(s1,i1,&v1,&as1); if (as1 == NIL) break;      
	    if (strcmp(as2->s,"cpcp*_*_*_") == OK) 
	       if (strcmp(as1->s,"cpcpcpcp*_") != OK) { i1++; continue; }
            r = siut_relation_atom_signature(as1,as2);
            if (r == relation) {
               if (HEAP) break; 
	       if ((r != GT) || (EPS_GT(v1,ZERO))) siut_set_signature(s1,i1,v1 + k * v2,as1);
               }
            i1++;
	    }
          if (HEAP) {
             if (as1 == NIL) as1 =(atom_signature_p)sial_create_atom_signature
                                  (HEAP,as2->dim,as2->dis,as2->s,(t_real)1,(t_real)1); 
	     if (v1 + k * v2) siut_set_signature(s1,i1,v1 + k * v2,as1);
             }
          i2++;
	  }
	return(s1);
	}

DEFINE	signature_p sisi_add_signature(HEAP,s1,k,s2,relation)
/**********************************************************
 s1 = s1 + k * s2.
 if relation(s1->as,s2->as) >= relation 
 s1->value += k * s2->value.
***********************************************************/
t_pointer	HEAP;
signature_p	s1,s2;
t_real		k;
t_err		relation;
{
	t_integer		i1 = 0,i2 = 0;
	t_real			v1,v2;
	atom_signature_p	as1,as2;
	t_bool			add = FALSE;
	t_err			r;

	LOOP {
	  siut_inq_signature(s2,i2,&v2,&as2); if (as2 == NIL) break; 
	  if (HEAP == NIL) if (EPS_LT(v2,ZERO)) { i2++; continue; }
	  i1 = 0; add = FALSE;
	  LOOP {
	    siut_inq_signature(s1,i1,&v1,&as1); if (as1 == NIL) break;
	    if   (strcmp(as1->s,as2->s) == OK) r = EQ;      
            else r = siut_relation_atom_signature(as1,as2);
	    /* r = siut_relation_atom_signature(as1,as2); */
            if (r == relation) {
               if (HEAP) break; else siut_set_signature(s1,i1,v1 + k * v2,as1);
               }
            i1++;
	    }
          if (HEAP) {
             if (as1 == NIL) as1 =(atom_signature_p)sial_create_atom_signature
                                  (HEAP,as2->dim,as2->dis,as2->s,(t_real)1,(t_real)1); 
	     if (v1 + k * v2) siut_set_signature(s1,i1,v1 + k * v2,as1);
             }
          i2++;
	  }
	return(s1);
	}

LOCAL	p_string 	string(hight,father,son,s)
/**********************************************************
 Recursive routine to compute the signature string for
 - an atom.
***********************************************************/
t_integer	hight;
atom_p		father,son;
p_string	s;
{
	t_integer	dim = sisi_inq_dimension(NIL);
	p_string	ss;
	t_integer	nb,NB,i;
	set_bond_p	SB;
	bond_p		b;

	if (son == NIL) return(s);
        if (hight > dim) return(s);
	if (father == son) return(s);
	strcpy(s,son->potential_type); s++; s++; 
        if (hight == dim) return(s); 
        SB = (set_bond_p)son->SB; nb = 0; ss = s;
	while (SB) {
	      b = SB->bond;
              if ((b->atom1 != son) && (b->atom1 != father)) 
                 if ((son) && (b->atom1)) ss = string(hight+1,son,b->atom1,s);
              if ((b->atom2 != son) && (b->atom2 != father))       
                 if ((son) && (b->atom2)) ss = string(hight+1,son,b->atom2,s);
              if (ss != s) nb++; SB = SB->succ; s = ss;
	      } /* nb is the number of bonds */

	/* look for radical until NB = son->degre */
	NB = ((hight == 0) ? son->degre : son->degre-1);
        for (i = nb; i < NB; i++) {
            *s = '.'; s++; *s = '_'; s++;
            } *s = '\0';
	nb = i;

        /* every atoms must have NB = VALENCE (4) - 1 new neighbour,
           the first atom has NB = VALENCE (4) neighbours
           if number of effective bond (nb) < NB the last 
           neighbour are = '*_'
        */
	NB = ((hight == 0) ? VALENCE : VALENCE-1);
        for (i = nb; i < NB; i++) {
            *s = '*'; s++; *s = '_'; s++;
            } *s = '\0';
	return(s);
	}

LOCAL	t_pointer atom_signature(atom,arg,distance)
/**********************************************************
 Compute the signature for
 - an atom.
 If the signature exist the atom signature is add in value
 if not a new signature is create with value 1
***********************************************************/
atom_p		atom;
arg_p		arg;
t_integer	distance;
{
LOCAL	t_string		S[MAXSTRING];
LOCAL	struct_signature_t     	s[2];
LOCAL	atom_signature_t	as;

	string(0,NIL,atom,&S[0]);
	as.dim = sisi_inq_dimension(NIL); 
 
	/* dis = distance(atom,root) */
        if (arg->root == NIL) as.dis = 0; else as.dis = distance;

        as.s = (t_string)(&S[0]); as.min = 1; as.max = 1;
        s[0].value = 1; s[0].as = &as; s[1].value = 0; s[1].as = NIL;

	/* add the new signature s */
	sisi_add_signature(arg->HEAP,arg->signature,(t_real)1,s,EQ);
        return((t_pointer)arg);	/* PSU Jan 92 */
	}

DEFINE	signature_p	sisi_bonding_site_signature(HEAP,bond,atom)
/***********************************************************
 Compute the signature for
 - a bonding site.
 The bonding site is a missing bond carried by an atom.
 If bond == NIL the atom doesn't carry the bond, if
 bond != NIL the bond is removed from the atom.
***********************************************************/
t_pointer	HEAP;
bond_p		bond;
atom_p		atom;
{
	t_integer	dim = sisi_inq_dimension(NIL);	
	molecule_p	molecule;
	signature_p	s;
	arg_t		arg;
	t_integer	size;
	atom_p		a1,a2;
EXTERN	signature_p	sisi_atom_signature();

	if (atom == NIL) return(NIL);
	if (bond) if ((bond->atom1 != atom) && (bond->atom2 != atom)) return(NIL);	
	molecule = (molecule_p)((group_p)atom->group)->molecule;

	/* computation of the size of the signature */
        dast_set_marq(molecule,NIL,NIL,FALSE);
        size = dast_connectivity_atom(atom,1,bond,NIL,NIL,dim-1);

	/* s = signature of the atom not containing the bond */
        if (bond) {
           a1 = bond->atom1; a2 = bond->atom2;
	   if (bond->order == 1) {bond->atom1 = NIL; bond->atom2 = NIL;}
           }
        arg.signature = (signature_p)sial_create_signature(HEAP,size); arg.HEAP = HEAP;
        dast_set_marq(molecule,NIL,NIL,FALSE);
        arg.root = atom; dast_connectivity_atom(atom,1,NIL,atom_signature,&arg,dim-1);
        dast_set_marq(molecule,NIL,NIL,FALSE);
        if (bond) { bond->atom1 = a1; bond->atom2 = a2; }
	return(arg.signature);
        }
	
DEFINE	signature_p	sisi_bond_signature(HEAP,bond)
/**********************************************************
 Compute the signature for
 - a bond.
***********************************************************/
t_pointer	HEAP;
bond_p		bond;
{
	t_integer	dim = sisi_inq_dimension(NIL);	
	molecule_p	molecule;
	signature_p	s,s1,s2;
	arg_t		arg;
	atom_p		a1,a2;
	t_integer	size;

	if ((bond->atom1 == NIL) || (bond->atom2 == NIL)) return(NIL);	
	molecule = (molecule_p)((group_p)bond->atom1->group)->molecule;

	/* computation of the size of the signature */
        dast_set_marq(molecule,NIL,NIL,FALSE);
        size = dast_connectivity_atom(bond->atom1,1,bond,NIL,NIL,dim-1) + 
               dast_connectivity_atom(bond->atom2,2,bond,NIL,NIL,dim-1);


	/* s1 = signature of the molecule containing the bond */
	s1 = (signature_p)sial_create_signature(HEAP,2 * size);
	arg.HEAP = HEAP; arg.signature = s1; 
        dast_set_marq(molecule,NIL,NIL,FALSE);
        arg.root = bond->atom1; dast_connectivity_atom(bond->atom1,1,bond,atom_signature,&arg,dim-1);
        arg.root = bond->atom2; dast_connectivity_atom(bond->atom2,2,bond,atom_signature,&arg,dim-1);
        s1 = arg.signature;

	/* s2 = signature of the molecule not containing the bond */

	a1 = bond->atom1; a2 = bond->atom2;
	bond->atom1 = NIL; bond->atom2 = NIL;
        s2 = (signature_p)sial_create_signature(HEAP,size);
       	arg.HEAP = HEAP; arg.signature = s2;
        dast_set_marq(molecule,NIL,NIL,FALSE);
        arg.root = a1; dast_connectivity_atom(a1,1,NIL,atom_signature,&arg,dim-1);
        arg.root = a2; dast_connectivity_atom(a2,2,NIL,atom_signature,&arg,dim-1);
        dast_set_marq(molecule,NIL,NIL,FALSE);
	bond->atom1 = a1; bond->atom2 = a2;
        s2 = arg.signature;
	s = sisi_add_signature(HEAP,s1,(t_real)-1,s2,EQ);
	return(s);
        }
		
DEFINE	signature_p	sisi_atom_signature(HEAP,atom)
/**********************************************************
 Compute the signature for
 - an atom.
***********************************************************/
t_pointer	HEAP;
atom_p		atom;
{
	t_integer	dim = sisi_inq_dimension(NIL);
	signature_p	signature;
	arg_t		arg;

	if (atom == NIL) return(NIL);
	signature = (signature_p)sial_create_signature(HEAP,1);
	arg.HEAP = HEAP; arg.signature = signature;
     	atom_signature(atom,&arg);
	return(arg.signature);
	}

DEFINE	signature_p	sisi_group_signature(HEAP,group)
/**********************************************************
 Compute the signature for
 - a group.
***********************************************************/
t_pointer	HEAP;
group_p		group;
{
	t_integer	dim = sisi_inq_dimension(NIL),n;
	signature_p	signature;
	arg_t		arg;

	if (group == NIL) return(NIL);
	signature = (signature_p)sial_create_signature(HEAP,dast_number_atom(NIL,group,NIL));
	arg.HEAP = HEAP; arg.signature = signature; arg.root = NIL;
	dast_all_atom(NIL,group,atom_signature,FALSE,&arg);
	return(arg.signature);
	}

DEFINE	signature_p	sisi_molecule_signature(HEAP,
					        molecule,atom)
/**********************************************************
 Compute the signature for
 - a molecule.
***********************************************************/
t_pointer	HEAP;
molecule_p	molecule;
atom_p		atom;
{
	t_integer	dim = sisi_inq_dimension(NIL), size = 0;
	signature_p	signature;
	arg_t		arg;

	if (molecule == NIL) return(NIL);
	if (atom) size = dast_number_atom(molecule,NIL,atom);
	else	  size = molecule->size;
	signature = (signature_p)sial_create_signature(HEAP,size);
	arg.HEAP = HEAP; arg.signature = signature; arg.root = NIL;
        dast_set_marq(molecule,NIL,NIL,FALSE);
	if (atom) {
           dast_set_marq(molecule,NIL,NIL,FALSE);
	   dast_connectivity_atom(atom,TRUE,NIL,NIL,molecule,MAXINTEGER);
	   dast_all_atom(molecule,NIL,atom_signature,TRUE,&arg);
           dast_set_marq(molecule,NIL,NIL,FALSE);
	   }
	else dast_all_atom(molecule,NIL,atom_signature,FALSE,&arg);     
	return(arg.signature);
	}

DEFINE	signature_p sisi_file_signature(HEAP,file)
/**********************************************************
 Read a signature from a crtc or scan file
***********************************************************/
t_pointer	HEAP;
t_bufstring	file;
{
	t_integer	dim = sisi_inq_dimension(NIL), size = 0, i = 0;
	float		v;
	char		s[MAXSIGSTRING];
	signature_p	signature;
	FILE		*FD;

	if ((FD = fopen(file,"r")) == NULL) return(NIL);
	while (fscanf(FD," %f ",&v) != EOF) { 
	  size++;
	  if (v < 1) break;
	  if (fscanf(FD,"%s",s) == EOF) {fclose(FD); return(NIL);}
	  }
	rewind(FD);
	signature = (signature_p)sial_create_signature(HEAP,size+1);

	while (fscanf(FD," %f ",&v) != EOF) { 
	  if (v < 1) break;
	  if (fscanf(FD,"%s",s) == EOF) {fclose(FD); return(NIL);}
	  signature[i].value = v;
	  signature[i].as = (atom_signature_p)sial_create_atom_signature(HEAP,dim,0,s,(t_real)0.0,(t_real)0.0);
	  i++;	   
	  }
	signature[i].value = -1;
	signature[i].as = NIL;
	fclose(FD);
	return(signature);
        }


DEFINE	signature_p sisi_signature(HEAP,molecule,group,atom,bond)
/**********************************************************
 Compute the signature for
 - a molecule,
 - a group,
 - an atom,
 - a bond.
***********************************************************/
t_pointer	HEAP;
molecule_p	molecule;
group_p		group;
atom_p		atom;
bond_p		bond;
{
	signature_p	signature = NIL;

	if (molecule)   signature = sisi_molecule_signature(HEAP,molecule,atom);
	else if (group) signature = sisi_group_signature(HEAP,group);
        else if (atom)  signature = sisi_atom_signature(HEAP,atom);
	else if (bond)  signature = sisi_bond_signature(HEAP,bond);
	return(signature);
        }

LOCAL	t_err	number_bonding_site(atom,arg)
/***********************************************************
 Return the number of bonding site of the atom having a
 signature equal to arg->signature.
 Actualy by equal I want to say 
 arg->signature is EQ to the atom signature
 arg->signature is LT to the atom signature.
***********************************************************/
atom_p		atom;
arg_p		arg;
{
	t_integer	n = 0;
	signature_p	s;
	t_err		relation;
	group_p		g;

        g = (group_p)atom->group;
	if (arg->group) {
	   if (strcmp(arg->group->comment,((group_p)atom->group)->comment)) return(OK);
           }
	if ((n = dast_number_bonding_site(NIL,NIL,atom)) == 0) return(OK);
        s = sisi_bonding_site_signature(arg->HEAP,NIL,atom);
        relation = sisi_relation_signature(arg->signature,s);
        if ((relation == EQ) || (relation == LT)) arg->number += n;
	return(OK);
	}

DEFINE	t_err	sisi_number_bonding_site(HEAP,molecule,group,
                                         atom,signature)
/***********************************************************
 Determine the number of bonding-site of a molecule a 
 molecular group !!!!! (molecular group can be formed
 by several groups)
 or an atom having a signature equal to signature
***********************************************************/
t_pointer	HEAP;
molecule_p	molecule;
group_p		group;
atom_p		atom;
signature_p	signature;
{
	arg_t	arg;

	arg.HEAP = HEAP; arg.signature = signature; arg.number = 0;
	if (group) molecule = (molecule_p)group->molecule; arg.group = group;
	if (atom) number_bonding_site(atom,&arg);
	else dast_all_atom(molecule,NIL,number_bonding_site,FALSE,&arg);
	return(arg.number);
	}

