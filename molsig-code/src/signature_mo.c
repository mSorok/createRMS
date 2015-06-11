 /***********************************************************
  Create a mol file our of a signature
***********************************************************/
#include <general.h>
#include <eps.h>
#include "signature.h"


LOCAL	t_integer signature_diameter(s)
/**********************************************************
 return signature diameter
***********************************************************/
p_string	s;
{
	t_integer	i, h = 0, hmax = 0, dim = 0, root = 0;
	for (i = 0; s[i] != '\0'; i++) {
	    if (s[i] == '[' && h == 0) root++;
	    if (s[i] == '(') h = h+1;
	    if (s[i] == ')') h = h-1;
	    if (h > hmax) hmax = h;
	    }
	dim = 2*hmax; if (root > 1) dim++;
	return(dim);
	}

LOCAL	signature_p create_signature(size)
/**********************************************************
 Create a signature. Initialize to NIL.
 of size s
***********************************************************/
t_integer	size;
{
	signature_p	signature;
	t_integer	i;

        signature = (signature_p)malloc((size+1)*sizeof(struct_signature_t));
        for (i = 0; i < size+1; i++) {
            signature[i].as = NIL;
			signature[i].value = 0;
	    }
	return(signature);
        }

LOCAL	atom_signature_p create_atom_signature(dim,dis,s,min,max)
/**********************************************************
 Create an atomic signature.
***********************************************************/
t_integer	dim,dis;
p_string	s;
t_real		min,max;
{
	atom_signature_p	signature = NIL;

	if (s == NIL) return(NIL);
	signature = (atom_signature_p)malloc(sizeof(atom_signature_t));
	signature->dim = dim; signature->dis = dis;
	signature->s = (p_string)malloc(strlen(s)+1);
	strcpy(signature->s,s);
	signature->min = min; signature->max = max;
	return(signature);
	}
    
DEFINE	signature_p simo_read_signature_file(file,dim)
/**********************************************************
 Read a signature from a scan file
 ***********************************************************/
t_bufstring	file;
t_integer	*dim;
{
	t_integer	size = 0, i = 0;
	float		v;
	char		s[MAXSIGSTRING];
	signature_p	signature;
	FILE		*FD;
	
	*dim = -1;
	if ((FD = fopen(file,"r")) == NULL) return(NIL);
	while (fscanf(FD," %f ",&v) != EOF) { 
		size++;
		if (v < 1) break;
		if (fscanf(FD,"%s",s) == EOF) {fclose(FD); return(NIL);}
		if (signature_diameter(s) > *dim) *dim = signature_diameter(s);
	}
	rewind(FD);
	if (size < 1) { fclose(FD); return(NIL); }
	signature = (signature_p)create_signature(size+1);
	
	while (fscanf(FD," %f ",&v) != EOF) { 
		if (v < 1) break;
		if (fscanf(FD,"%s",s) == EOF) {fclose(FD); return(NIL);}
		signature[i].value = v;
		signature[i].as = (atom_signature_p)create_atom_signature(\
				*dim/2,0,s,(t_real)0.0,(t_real)0.0);
		i++;	   
	}
	signature[i].value = -1;
	signature[i].as = NIL;
	fclose(FD);
	return(signature);
}


LOCAL	t_integer	SIZE = 0;
LOCAL	t_integer	ID = 0;
typedef	struct	ARGNAME {
                atom_p		atom;
		t_name		name;
		} argname_t, *argname_p;

LOCAL	t_err	atom_name(atom,arg)
/***********************************************************
 Return the number of bonding site of the atom having a
 signature equal to arg->signature.
***********************************************************/
atom_p		atom;
argname_p	arg;
{
	if (atom->ID < SIZE) return(OK);
	if (strcmp(atom->name,arg->name) == OK) {
	   arg->atom = atom;
	   return(ERROR);
	   }
	return(OK);
	}

LOCAL	atom_p	signature_atom(molecule,s,I,ap)
/***********************************************************
 Build a molecule from an atomic signature string
***********************************************************/
molecule_p	molecule;
p_string	s;
t_integer	*I;
atom_p		ap;
{
	t_integer	i,j,k;
	t_bool		order = 1;
	t_name		name,element;
	atom_p		a = NIL;
	t_coord	    ZEROC = {0,0,0};
	t_bufstring	comment;
	argname_t	arg;
	t_real		charge = 0;
	
	if (s == NIL) return(NIL); i = 0;	
	if (s[0] == '=') { order = 2; i++; }
	if (s[0] == 't') { order = 3; i++; }
	if (s[0] == 'p') { order = 4; i++; }
	if (s[i] != '[') return(NIL); i++; j = 0; 
	while (s[i] != ']') {element[j] = s[i]; i++; j++;}
	element[j] = '\0'; strcpy(name,element);
	for (j = 0; element[j] != '\0'; j++) if (element[j] == ',') break;
	if (element[j] == ',') {
	   strcpy(arg.name,name); arg.atom = NIL;
	   dast_all_atom(molecule,NIL,atom_name,FALSE,&arg);
	   a = arg.atom; element[j] = '\0';
	   }
	// charge
	for (j = 0; element[j] != '\0'; j++) 
		if ((element[j] == '+') || (element[j] == '-') || (element[j] == ':')) break;
	if (j != strlen(element)) {
		t_integer	c;
		if (element[j] == ':') c == 4;
		else c = atoi(&element[j+1]);
		if (element[j] == '-') c = -c;
		// printf("before charge = %d\n",c);
		if (c ==  2)      charge = 2; 
		else if (c ==  1) charge = 3; 
		else if (c == -1) charge = 5; 
		else if (c == -2) charge = 6; 
		else if (c == -3) charge = 7; 
		// printf("after name=%s element=%s charge = %f\n",name,element,charge);
	}
	strcpy(comment,name);
	if (a == NIL) a = (atom_p)daal_create_atom(molecule->HEAP,name,ID++,\
	molecule->SG->group,NIL,NIL,NIL,ZEROC,element,charge,0,comment,FALSE);
	molecule->size = ID;
	if (ap) {
		daal_create_bond(molecule->HEAP,a,ap,NIL,NIL,ZEROC,ZEROC,order,"",FALSE);
		ap->degre += order; a->degre += order;
	}
	i++;
	if (s[i] == '\0') {*I += i;   return(a);}
	if (s[i] == '[')  {*I += i;   return(a);}
	if (s[i] == '=')  {*I += i;   return(a);}
	if (s[i] == 't')  {*I += i;   return(a);}
	if (s[i] == 'p')  {*I += i;   return(a);}
	if (s[i] == ')')  {*I += i;   return(a);}
	if (s[i] != '(')  {xerr_c("reading signature string\n");}
	i++; k = 0;
	LOOP {
          if ((s[i] != '[') && (s[i] != '=') &&
              (s[i] != 't') && (s[i] != 'p')) break;
	  signature_atom(molecule,&s[i],&i,a);
	  }
	*I += i+1;
	return(a);
	}

DEFINE	void simo_reset_molecule(molecule)
/**********************************************************
 ***********************************************************/
molecule_p		molecule;
{
	molecule->SG->group->SA = NIL;
	molecule->SG->group->SB = NIL;
	molecule->size = SIZE = ID = 0;
}

DEFINE	atom_p simo_signature_to_atom(molecule,as)
/**********************************************************
 Add an atom to a molecule according to its signature
***********************************************************/
molecule_p		molecule;
atom_signature_p	as;
{
	atom_p		atom = NIL;
	p_string	s,s2;
	t_integer	I = 0;

	if (as == NIL) return(atom);
	if ((s = as->s) == NIL) return(atom);
	atom = signature_atom(molecule,s,&I,NIL);
	return(atom);
}

DEFINE	molecule_p	simo_signature_to_molecule(file)
/******************************************************************************
 Return a molecule from a signature
 ******************************************************************************/
t_bufstring	file;
{
	LOCAL	t_pointer		HEAP = NIL;
	t_integer				Na = 0,i,j,hi = -1;
	signature_p				signature;
	molecule_p				MF;
	
	signature = (signature_p)simo_read_signature_file(file,&hi);
	if (signature == NIL) xerr_c("no input signature");
	
	/* compute the number of atoms */
	for (i = 0; signature[i].as; i++) Na += signature[i].value;

	/* create file and solution molecules  */
	if ((HEAP = (t_pointer)heat_creer_hea(MAXHEAP)) == NIL) 
	   xerr_c("no more space");
	MF = (molecule_p)daal_create_molecule(HEAP,file,0,0,NIL,NIL,NIL,NIL,"",0);
	daal_create_group(MF->HEAP,MF->name,1,MF,0,NIL,NIL,MF->name,FALSE,NIL); 
	for (i = 0; signature[i].as; i++) {
	  if (signature[i].as == NIL || signature[i].value < 1) break;
	  for (j = 0; j < (t_integer)signature[i].value; j++) simo_signature_to_atom(MF,signature[i].as);
	}
	return(MF);
}







