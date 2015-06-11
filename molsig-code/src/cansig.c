/***********************************************************
 Compute canonized signature for an atom and a molecule
***********************************************************/
#include <general.h>
#include <eps.h>
#include "signature.h"
EXTERN	t_pointer	sial_create_signature();
EXTERN	t_integer	sisi_inq_dimension();

typedef	char	element_t[MAXSTRING/10];

typedef	struct	ARG {
		t_pointer	HEAP;
		group_p		group;
		atom_p		root;
		signature_p	signature;
		t_integer	number;
		} arg_t, *arg_p;

LOCAL	t_integer	MINC,MAXC,MINXL,MAXXL;
LOCAL	t_err		GRAPH = 0; /* 0 bounded valence acyclic */
                                   /* 1 bounded valence cyclic */
                                   /* 2 protein networks */

LOCAL	signature_p create_signature(HEAP,size)
/**********************************************************
 Create a signature. Initialize to NIL.
 No more than MAXSIG elements
***********************************************************/
t_pointer	HEAP; 
t_integer	size;
{
	signature_p	signature;
	t_integer	i,s;

        s = (size+1) * sizeof(struct_signature_t);
        signature = (signature_p)malloc(s);
        for (i = 0; i < size+1; i++) siut_set_signature(signature,i,ZERO,NIL);
		return(signature);
        }


LOCAL	atom_signature_p create_atom_signature(HEAP,dim,dis,s,min,max)
/**********************************************************
 Create an atomic signature.
***********************************************************/
t_pointer	HEAP; 
t_integer	dim,dis;
p_string	s;
t_real		min,max;
{
	atom_signature_p	signature = NIL;

	if (HEAP  == ((t_pointer)ERROR)) return((atom_signature_p)ERROR);
	if (s == NIL) return(NIL);
	signature = (atom_signature_p)malloc(sizeof(atom_signature_t));
	signature->dim = dim; signature->dis = dis;
	signature->s = (p_string)malloc(strlen(s)+1);
	strcpy(signature->s,s);
	signature->min = min; signature->max = max;
	return(signature);
	}

LOCAL	signature_p add_signature(HEAP,s1,k,s2,relation)
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
            else r = NE;
            if (r == relation) {
               if (HEAP) break; else siut_set_signature(s1,i1,v1 + k * v2,as1);
               }
            i1++;
	    }
          if (HEAP) {
             if (as1 == NIL) {
                as1 =(atom_signature_p)create_atom_signature (HEAP,as2->dim,as2->dis,as2->s,(t_real)1,(t_real)1); 
               }
	     if (v1 + k * v2) siut_set_signature(s1,i1,v1 + k * v2,as1);
             }
          i2++;
	  }
	return(s1);
	}

LOCAL	t_pointer atom_signature(atom,arg,distance)
/**********************************************************
 Compute the signature for
 - an atom.
 If the signature exist the atom signature is add in value
 if not a new signature is create with value 1
 If dimension id ODD then the signature is computed
 for all bond that originated @ atom
***********************************************************/
atom_p		atom;
arg_p		arg;
t_integer	distance;
{
EXTERN	t_err				dast_set_ID_atom();
LOCAL	p_string			S;
LOCAL	struct_signature_t  s[2];
LOCAL	atom_signature_t	as;
		t_name			pt;
		molecule_p		molecule = (molecule_p)((group_p)atom->group)->molecule;
		t_integer		ID = 0, dim = sisi_inq_dimension(NIL);
		set_bond_p		SB;  
	

	as.dim = dim; 	
	if (arg->root == NIL) as.dis = 0; else as.dis = distance;
	if (EVEN(dim)) {
		S = (p_string)sicd_signature_atom(atom,NIL,dim/2);
		if (S!=NIL) { // PC 05/12
			as.s = (t_string)(&S[0]); as.min = 1; as.max = 1;
			s[0].value = 1; s[0].as = &as; s[1].value = 0; s[1].as = NIL;
			/* add the new signature s */
			add_signature(arg->HEAP,arg->signature,(t_real)1,s,EQ);
			free(as.s);
		}
		return((t_pointer)arg);	/* PSU Jan 92 */
	}
	
	SB = (set_bond_p)atom->SB;
	while (SB) {
	    bond_p	b = SB->bond;
	    atom_p	aa = ((b->atom1 == atom)?b->atom2:NIL);
		if (aa == NIL) { SB = SB->succ; continue; }
		S = (p_string)sicd_signature_atom(atom,aa,dim/2);
		if (S!=NIL) { // PC 05/12
			as.s = (t_string)(&S[0]); as.min = 1; as.max = 1;
			s[0].value = 1; s[0].as = &as; s[1].value = 0; s[1].as = NIL;
			/* add the new signature s */
			add_signature(arg->HEAP,arg->signature,(t_real)1,s,EQ);
			free(as.s);
		}
	    SB = SB->succ;
	}
	return((t_pointer)arg);	/* PSU Jan 92 */	
}

	
DEFINE	signature_p	sisc_bond_signature(HEAP,bond)
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
	s = add_signature(HEAP,s1,(t_real)-1,s2,EQ);
	return(s);
	}
		
DEFINE	signature_p	sisc_atom_signature(HEAP,atom)
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

DEFINE	signature_p	sisc_group_signature(HEAP,group)
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

DEFINE	signature_p	sisc_molecule_signature(HEAP,
					        molecule,atom)
/**********************************************************
 Compute the signature for
 - a molecule 
***********************************************************/
t_pointer	HEAP;
molecule_p	molecule;
atom_p		atom;
{
	t_integer	dim = sisi_inq_dimension(NIL), size = 0;
	signature_p	signature;
	arg_t		arg;

	if (molecule == NIL) return(NIL);

	if (molecule->size < 1) return(NIL);
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

LOCAL	int	sigcmp(s1,s2)
/***********************************************************
 compare s1 s2 
***********************************************************/
signature_p	s1,s2;
{
	if ((s1 == NIL) && (s2 == NIL)) return(0);
	if (s2 == NIL) return(-1);
	if (s1 == NIL) return( 1);
	if (s2->as == NIL) return(-1);
	if (s1->as == NIL) return( 1);
	return(-strcmp(s1->as->s,s2->as->s));
	}

DEFINE	signature_p sisc_signature(HEAP,molecule,group,atom,bond)
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
	GRAPH = 1;
	
	if (molecule)   signature = sisc_molecule_signature(HEAP,molecule,atom);
	else if (group) signature = sisc_group_signature(HEAP,group);
	else if (atom)  signature = sisc_atom_signature(HEAP,atom);
	else if (bond)  signature = sisc_bond_signature(HEAP,bond);
	else return(signature);
		
	/* sort signature */
	qsort(signature,sial_size_signature(signature),\
              sizeof(struct_signature_t),sigcmp);
	return(signature);
	}

#define	MAXPT		200
LOCAL	t_name		CANPT;
LOCAL	t_name *	ATOMPT;
typedef	struct		STRUCT_CLASS {
		t_integer	n,c,l;
		p_string	s;
		} class_t, *class_p;
LOCAL	p_string	SMAX = NIL;
LOCAL   int	cmp_atom_potential_type(apt1,apt2)
/***********************************************************
 ***********************************************************/
t_name apt1, apt2;
{
	return(strcmp(apt2,apt1));
}

LOCAL	t_err 		atom_potential_type(atom)
/**********************************************************
 Compute and sort the pt of atom and its neighbors
 ***********************************************************/
atom_p		atom;
{	
	set_bond_p	sb;
	t_integer	n = 0,i;
	t_name		APT[MAXPT];
	sb = (set_bond_p)atom->SB; strcpy(ATOMPT[atom->ID],"");
	if (sb != NIL && sb != (set_bond_p)ERROR) strcpy(ATOMPT[atom->ID],atom->potential_type);
	while(sb) { 
		if (sb == (set_bond_p)ERROR) break;
		bond_p	b = sb->bond;
		atom_p	a = ((b->atom1 == atom) ? b->atom2 : b->atom1);
		strcpy(APT[n],"2"); strcat(APT[n++],a->potential_type);
		set_bond_p	ssb = (set_bond_p)a->SB;
		while (ssb) {
			if (ssb == (set_bond_p)ERROR) break;
			bond_p	bb = ssb->bond;
			atom_p	aa = ((bb->atom1 == a) ? bb->atom2 : bb->atom1);
			if (aa != atom) { strcpy(APT[n],"1"); strcat(APT[n++],aa->potential_type); }
			set_bond_p	sssb = (set_bond_p)aa->SB;
			while (sssb) {
				if (sssb == (set_bond_p)ERROR) break;
				bond_p	bbb = sssb->bond;
				atom_p	aaa = ((bbb->atom1 == aa) ? bbb->atom2 : bbb->atom1);
				if (aaa != a) { strcpy(APT[n],"0"); strcat(APT[n++],aaa->potential_type); }
				sssb = sssb->succ; 
			}
			ssb = ssb->succ; 
		}
		sb = sb->succ; 
	}
	strcpy(APT[n],"");
	qsort(APT,n,sizeof(t_name),cmp_atom_potential_type);
	for (i = 0; i < n; i++) strcat(ATOMPT[atom->ID],APT[i]);
	return(OK);
}

LOCAL	t_err 		print_atom_signature(atom,CLASS)
/**********************************************************
 Recursive routine to compute the signature string for
 - an atom.
***********************************************************/
atom_p		atom;
class_p		CLASS;
{
	p_string	s = NIL;
	t_integer	n,i;
	t_integer *	label;

	n = ((molecule_p)(((group_p)atom->group)->molecule))->size;
	CLASS[atom->ID].n = atom->ID;
	if (strcmp(ATOMPT[atom->ID],CANPT) != 0) return(OK);
	label = (t_integer *)malloc((n+1)*sizeof(t_integer));
 	s = (p_string)sicd_signature_atom_label(atom,sisi_inq_dimension(NIL)/2,label);
	CLASS[atom->ID].s = s;
	if (SMAX) if (strcmp(s,SMAX) < 0) {
	   free(label); return(OK);
	   }
	SMAX = s; for (i = 0; i < n; i++) CLASS[i].l = label[i];
	free(label);
	return(OK);
	}

LOCAL	t_err 		label_atom_signature(atom,CLASS)
/**********************************************************
 Recursive routine to label an atom
***********************************************************/
atom_p		atom;
class_p		CLASS;
{
	atom->ID = CLASS[atom->ID].l;
	return(OK);
	}

LOCAL	int	classcmp(c1,c2)
/***********************************************************
 compare c1 c2 
***********************************************************/
class_p	c1,c2;
{
	if ((c1 == NIL) && (c2 == NIL)) return(0);
	if (c2 == NIL) return(-1);
	if (c1 == NIL) return( 1);
	if (c2->s == NIL) return(-1);
	if (c1->s == NIL) return( 1);
	return(-strcmp(c1->s,c2->s));
	}

LOCAL	t_err 		canonical_atom_potential_type(atom)
/**********************************************************
 Search the canonical potential type
 ***********************************************************/
atom_p		atom;
{
	t_name	pt;
	atom_potential_type(atom,pt);
	if (strcmp(pt,CANPT) > 0) strcpy(CANPT,pt);
	return(OK);
}

DEFINE	t_err sisc_canonize(file,molecule)
/**********************************************************
 Compute a canonical representation for a molecule
 and relabel the atom of that molecule
***********************************************************/
t_name		file;
molecule_p	molecule;
{
	signature_p	signature = NIL;
	t_integer	i,n = molecule->size;
	class_p		CLASS;
	FILE		*FD;
	
	GRAPH = 1; SMAX = NIL;
	CLASS = (class_p)malloc((n+1)*sizeof(class_t));
	ATOMPT = (t_name *)malloc((n+1)*sizeof(t_name));
	t_name * atompt = (t_name *)malloc((n+1)*sizeof(t_name));
	for (i = 0; i < n; i++) {
	    CLASS[i].s = NIL;
	    CLASS[i].n = CLASS[i].c = CLASS[i].l = -1;
	    }
	dast_all_atom(molecule,NIL,canonical_atom_potential_type,FALSE,NIL);
	for (i = 0; i < n; i++) strcpy(atompt[i],ATOMPT[i]);
	qsort(atompt,n,sizeof(t_name),cmp_atom_potential_type);
	t_integer ni = 1, nmin = n, imin = 0;
	for (i = 1; i < n; i++) 
//		if (strcmp(ATOMPT[i],ATOMPT[i-1]) == 0) ni++;  // PC Changed to atompt
		if (strcmp(atompt[i],atompt[i-1]) == 0) ni++;
		else {
			if (ni < nmin) { nmin = ni; imin = i-1; }
			ni = 1;
			}
	strcpy(CANPT,atompt[imin]);free(atompt);
    if (VERBOSE) printf("%d %s : ",nmin,CANPT);
	dast_all_atom(molecule,NIL,print_atom_signature,FALSE,CLASS);
	dast_all_atom(molecule,NIL,label_atom_signature,FALSE,CLASS);
	/*
	qsort(CLASS,n,sizeof(class_t),classcmp);
	CLASS[0].c = 0;
	for (i = 1; i < n ; i++) {
	  if (classcmp(&CLASS[i-1],&CLASS[i]) == 0) 
	     CLASS[i].c = CLASS[i-1].c;
	  else CLASS[i].c = CLASS[i-1].c + 1;
	  }
	if (VERBOSE) {
	for (i = 0; i < n; i++)
	   printf("CLASS[%d] = %d (%s)\n",CLASS[i].n,CLASS[i].c,CLASS[i].s);
	   }
	printf("(%d",CLASS[0].n);
	for (i = 1; i < n; i++) 
	  if (CLASS[i].c == CLASS[i-1].c) printf(" %d",CLASS[i].n);
	  else                            printf(")\n(%d",CLASS[i].n);
	printf(")\n");
	*/
	if ((FD = fopen(file, "w")) == NULL) return (ERROR);
	fprintf(FD,"1.0 %s\n0.0\n",SMAX);
	fclose(FD);
	/* 
	printf("( ",molecule->name);
	for (i = 0; i < n; i++) printf("%d ",CLASS[i].l);
	printf(")\n");
	 */
    
	for (i = 0; i < n; i++) free(CLASS[i].s);
	free(CLASS); free(ATOMPT);
	return(OK);
	}

DEFINE	t_err sisc_orbit(molecule,ORBIT)
/**********************************************************
 Compute the atoms orbits of molecule
***********************************************************/
molecule_p	molecule;
t_integer	ORBIT[];
{
	signature_p	signature = NIL;
	t_integer	i,n = molecule->size;
	class_p		CLASS;

	if (ORBIT == NIL) return(ERROR);
	GRAPH = 1; SMAX = NIL;
	CLASS = (class_p)malloc((n+1)*sizeof(class_t));
	for (i = 0; i < n; i++) {
	    CLASS[i].s = NIL;
	    CLASS[i].n = CLASS[i].c = CLASS[i].l = -1;
	    }
	dast_all_atom(molecule,NIL,print_atom_signature,FALSE,CLASS);
	qsort(CLASS,n,sizeof(class_t),classcmp);
	CLASS[0].c = 0;
	for (i = 1; i < n ; i++) {
	  if (classcmp(&CLASS[i-1],&CLASS[i]) == 0) 
	     CLASS[i].c = CLASS[i-1].c;
	  else CLASS[i].c = CLASS[i-1].c + 1;
	  }
	/*
	for (i = 0; i < n; i++)
	   printf("CLASS[%d] = %d (%s)\n",CLASS[i].n,CLASS[i].c,CLASS[i].s);
	*/
	for (i = 0; i < n; i++) ORBIT[CLASS[i].n] = CLASS[i].c;
	for (i = 0; i < n; i++) free(CLASS[i].s);
	free(CLASS);
	return(OK);
	}
