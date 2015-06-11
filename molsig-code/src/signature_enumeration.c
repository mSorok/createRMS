 /***********************************************************
  Labeled enumeration from signatures
  can enumerate signature of any height
  from signature of any height
  work with scan and sig
  does not work with odd diameter
  the number of desired connected components
  is given by the user (0 for any)
  co. Jean-Loup Faulon
  Jan-March 2011
  Jan-March 2012
***********************************************************/
#include <general.h>
#include <eps.h>
#include "signature.h"

EXTERN	t_bool		COUNT;
LOCAL	t_integer	NBCC; // nbr of connected components

typedef struct		BONDINGSITE_STR {
			t_integer	x,y,index;
			t_bool		o;
			p_string	s;
			} bonding_site_t, *bonding_site_p;

LOCAL	t_integer	m_ID = 1,g_ID,a_ID = 1;
LOCAL	bond_p	create_bond(a1,a2,order)
/**********************************************************
  Create a bond
***********************************************************/
atom_p		a1,a2;
t_bool		order;
{
	bond_p		bond = NIL;
	set_bond_p	sb;
	t_coord     ZEROC = {0,0,0};

	if ((a1 == NIL) && (a2 == NIL)) return((bond_p)ERROR);

	/* the bond is created only if it doesn't exist already */
	if ((bond = (bond_p)dast_inq_bond(a1,a2))) return(bond);
	bond = (bond_p)malloc(sizeof(bond_t));

	bond->atom1 = a1; bond->atom2 = a2;
	bond->order = order;
	if (a1) { 
	   sb = (set_bond_p)malloc(sizeof(set_bond_t));
	   sb->bond = bond; sb->succ = a1->SB; a1->SB = sb;
	}
	if (a2) {
	   sb = (set_bond_p)malloc(sizeof(set_bond_t));
	   sb->bond = bond; sb->succ = a2->SB; a2->SB = sb;
	}
	strcpy(bond->comment,""); 
	bond->marq = FALSE; bond->visit = FALSE;
	bond->geom1 = ZEROC; bond->geom2 = ZEROC; 
	return(bond);
}

LOCAL	t_err	delete_bond(bond)
/**********************************************************
  Delete a bond, and create the bonding sites.
***********************************************************/
bond_p	bond;
{
	set_bond_p	sbs,sbp;

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
		    free(sbs);
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
			free(sbs);
		}
	   bond->atom2 = NIL;
	}
	free(bond);
 	return(OK);
	}

/*--------------------------------------------------------------------
  Enumeration algorithm 
---------------------------------------------------------------------*/
LOCAL  int		str_cmp(p_string *s1,p_string *s2) {return(strcmp(*s2,*s1));}
LOCAL  p_string	compute_signature_atom_target(a,SAF)
/**********************************************************************
 Compute the target signature for atom x
 It's the target signature of a + the sorted target signature
 of all atoms bonded to a
 **********************************************************************/
p_string *		SAF;
atom_p			a;
{
	t_integer	x = a->ID;
	set_bond_p	SB = a->SB;
	t_integer	i = 0,n = 0;
	p_string	s,*sx;
	
	SB = a->SB; while (SB) { i++; SB = SB->succ; } SB = a->SB; 
	sx = (p_string *)malloc((i)*sizeof(p_string)); 
	n = strlen(SAF[x]); i = 0;
	while (SB) {
		atom_p	y = (SB->bond->atom1 == a)?SB->bond->atom2:SB->bond->atom1;
		if (y) { 
			sx[i] = malloc((strlen(SAF[y->ID])+4)*sizeof(char)); 
			sprintf(sx[i]," %1d %s",SB->bond->order,SAF[y->ID]);
			n += strlen(sx[i]); 
			i++; 
		}
		SB = SB->succ;
	}
	s = (p_string)malloc((n+1)*sizeof(char)); strcpy(s,"");
	n = i; 	qsort(sx,n,sizeof(p_string),str_cmp);
	strcat(s,SAF[x]); for (i = 0; i < n; i++) { 
		strcat(s,sx[i]); free(sx[i]);
	} free(sx);
	return(s);
}

LOCAL  p_string	compute_signature_atom_target_2(a,SAF)
/**********************************************************************
 Compute the target signature for atom x
 It's the target signature of a + the sorted target signature
 of all atoms bonded to a where a is the atom1 of the bond
 (originated from a)
 **********************************************************************/
p_string *		SAF;
atom_p			a;
{
	t_integer	x = a->ID;
	set_bond_p	SB = a->SB;
	t_integer	i = 0,n = 0;
	p_string	s,*sx;
	
	SB = a->SB; while (SB) { i++; SB = SB->succ; } SB = a->SB; 
	sx = (p_string *)malloc((i)*sizeof(p_string)); 
	n = strlen(SAF[x]); i = 0;
	while (SB) {
		atom_p	y = SB->bond->atom2;
		if (y && y != a) { 
			sx[i] = malloc((strlen(SAF[y->ID])+4)*sizeof(char)); 
			sprintf(sx[i]," %1d %s",SB->bond->order,SAF[y->ID]);
			n += strlen(sx[i]); 
			i++; 
		}
		SB = SB->succ;
	}
	s = (p_string)malloc((n+1)*sizeof(char)); strcpy(s,"");
	n = i; 	qsort(sx,n,sizeof(p_string),str_cmp);
	strcat(s,SAF[x]); for (i = 0; i < n; i++) { 
		strcat(s,sx[i]); free(sx[i]);
	} free(sx);
	return(s);
}

LOCAL  p_string	compute_signature_atom(a,h,scan)
/**********************************************************************
 return the signature height h of atom a
 simple form if scan == FALSE
 **********************************************************************/
atom_p		a;
t_integer	h;
t_bool		scan;
{
EXTERN	t_bufstring	SCANTYPE;
	t_bufstring	scantype;
	p_string	s;
	strcpy(scantype,SCANTYPE);
	sicd_signature_reset();
	if (scan) strcpy(SCANTYPE,"scan"); else strcpy(SCANTYPE,"sig");
	s = (p_string)sicd_signature_atom(a,NIL,h);
	strcpy(SCANTYPE,scantype);
	return(s);
}

typedef	struct	ARGCAN {
	// the upper and lower signature of b
	t_integer	h;
	t_bool		scan;
	p_string	s; 
} arg_can_t, *arg_can_p;

LOCAL  arg_can_p	signature_atom_canonical(a,arg)
/**********************************************************************
 return the signature height h of atom a
 simple form if scan == FALSE
 **********************************************************************/
atom_p		a;
arg_can_p	arg;
{
	p_string	s = compute_signature_atom(a,arg->h,arg->scan);
	if (strcmp(s,arg->s) < 0) {
		free(arg->s); arg->s = s;
		return(arg);
	}
	free(s);
	return(arg);
}

LOCAL  p_string	compute_signature_atom_canonical(a,h,scan)
/**********************************************************************
 compute the canonical signature for all atoms attached to a
 **********************************************************************/
atom_p		a;
t_integer	h;
t_bool		scan;
{
	arg_can_t	arg;
	if (a == NIL) return(NIL);
	arg.h = h; arg.scan = scan;
	arg.s = compute_signature_atom(a,arg.h,arg.scan);
	dast_connectivity_atom(a,TRUE,NIL,signature_atom_canonical,&arg,MAXATOM);
	return(arg.s);
}

LOCAL  p_string	compute_signature_molecule(A,M,h,scan)
/**********************************************************************
 compute the canonical signature for all atoms
 the molecular signature is spliced into connected components
 **********************************************************************/
atom_p	*	A;
molecule_p	M;
t_integer	h; 
t_bool		scan;
{
	p_string *	CAN;
	p_string *	CANCC;
	p_string *	CANcc;
	p_string	can,CANv;
	t_bufstring	nb; 
	t_integer	n,N = M->size, *NCC, nbcc, ncc,x,xcc,nxcc;
	t_real		v;
	
	if (N < 1) return(NIL);
	NCC = (t_integer *)malloc((N+1)*sizeof(t_integer));
	nbcc = dast_connectivity_distribution(M,NIL,FALSE,NCC);
	if (NBCC > 0 && nbcc != NBCC) { free(NCC); return(NIL); }
	CAN = (p_string *)malloc((N+1)*sizeof(p_string));
	CANcc = (p_string *)malloc((N+1)*sizeof(p_string));
	CANCC = (p_string *)malloc((nbcc+1)*sizeof(p_string));
	for (x = 0; x < N; x++) CAN[x] = compute_signature_atom(A[x],h,scan);

	/* For each CC */
	for (ncc=1; ncc <= nbcc; ncc++) {
		n = 0; nxcc = 0; 
		for (x = 0; x < N; x++) {
			if (NCC[x] != ncc) continue;
			CANcc[nxcc++] = CAN[x];
			n+= (strlen(CAN[x])+10);
		}
		qsort(CANcc,nxcc,sizeof(p_string),str_cmp);
		can = (p_string)malloc((n+1)*sizeof(char)); can[0] = '\0';
		v = 1; CANv = CANcc[0];
		for (x = 1; x < nxcc; x++) {
			if (strcmp(CANcc[x],CANv) == 0) { v++; continue; }
			sprintf(nb,"%.1f",v);
			strcat(can,nb); strcat(can," "); strcat(can,CANv); strcat(can,"\n");
			v = 1; CANv = CANcc[x];
		}
		sprintf(nb,"%.1f",v);
		strcat(can,nb); strcat(can," "); strcat(can,CANv); strcat(can,"\n");
		CANCC[ncc-1] = can;
	}
	for (x = 0; x < N; x++) free(CAN[x]); free(CAN); free(CANcc); free(NCC);
	qsort(CANCC,nbcc,sizeof(p_string),str_cmp);
	n = 0; for (ncc = 0; ncc < nbcc; ncc++) n = n+strlen(CANCC[ncc])+2;
	can = (p_string)malloc((n+1)*sizeof(char)); can[0] = '\0';
	for (ncc = 0; ncc < nbcc; ncc++) { 
		strcat(can,CANCC[ncc]); 
		if (nbcc > 1 && ncc != nbcc-1) strcat(can,"\n"); 
		free(CANCC[ncc]); 
	}
	free(CANCC);
	return(can);
}

LOCAL	t_integer 	index(A,x)
/**********************************************************
 Return the first bonding site number for atom x
***********************************************************/
atom_p		*A;
t_integer	x;
{
	t_integer	i,I = 0;
	for (i = 0; i < x; i++) I += A[i]->degre;
	return(I);
	}

LOCAL  int	compare_bonding_site(b1,b2)
/**********************************************************************
 Compare two signature
**********************************************************************/
bonding_site_p		b1,b2;
{
	p_string	s1,s2;
	t_integer	c;

	if (b1 == NIL && b2 == NIL) return(0);
	if (b1 == NIL) return( 1);
	if (b2 == NIL) return(-1);
	if (b1->s == NIL && b2->s == NIL) return(0);
	if (b1->s == NIL) return( 1);
	if (b2->s == NIL) return(-1);
	return(str_cmp(&(b1->s),&(b2->s)));
	}

EXTERN	t_integer	MAXSOL;
LOCAL	t_bufstring	INPUTFILE;
LOCAL	t_bool		END_ENUMERATE = FALSE;	
LOCAL	t_bool		HIT = FALSE;
LOCAL	t_integer	Nsol = -1;
LOCAL	t_integer	Nenum = -1;
LOCAL	t_integer	Ntry = 0;
LOCAL	t_err		end_enumerate(A,SF,hi,scani,ho,scano)
/**********************************************************************
 Check that molecule match the provided signature
**********************************************************************/
atom_p		*A;
p_string	*SF;
t_integer	hi,ho;
t_bool		scani,scano;
{
LOCAL	p_string *	SOL;
EXTERN	t_real	dast_weight_element_atom();
	p_string	s = NIL;
	molecule_p	M = ((molecule_p)(((group_p)(A[0]->group))->molecule));
	t_integer	x;
	t_bufstring	file,nb;
	set_molecule_t	SM;
	FILE			*FD;
	t_real			mass = 0;
	
	if (Nsol == -1) {
	   SOL = (p_string *)malloc((MAXSOL+1)*sizeof(p_string)); Nsol = 0;
	   }
	Ntry++;

	if (Nsol >= MAXSOL || Ntry >= MAXSOL) { 
	   END_ENUMERATE = TRUE;
	   return(ERROR);
	   }
	
//	if (dast_connectivity(M,NIL,FALSE) > 1) JLF 20/2012
//	   xerr_c("bug");

	/* canonize structure and compute the mass */
	for (x = 0; x < M->size; x++) {
		mass +=  dast_weight_element_atom(A[x]);
	    if (A[x]->degre > 1) {
	       s = (p_string)compute_signature_atom(A[x],hi,scani);
	       if (strcmp(s,SF[x]) != OK) { 
			   free(s); return(ERROR);
		   }
	       free(s);
	       }
	    }
	if (Ntry % 10000 == 0) {
		fprintf(stderr,"%s Ntry = %d Nsol = %d MAXSOL = %d HIT = %d\n",M->name,Ntry,Nsol,MAXSOL,HIT);
	}
	if ((s = (p_string)compute_signature_molecule(A,M,ho,scano)) == NIL) return(ERROR); 
	for (x = 0; x < Nsol; x++) if (strcmp(s,SOL[x]) == OK) break;
	if (x < Nsol) { free(s); return(ERROR); }
	SOL[Nsol++] = s;
	if (HIT == FALSE) if (hit(s)) HIT = TRUE;
	if (COUNT == TRUE) return(OK);

	strcpy(file,M->name); strcat(file,"_out"); 
	sprintf(nb,"%d",Nsol); strcat(file,nb);
	if (scano) strcat(file,".scan"); else strcat(file,".sig");
	sprintf(nb,"%d",2*ho); strcat(file,nb);
	if ((FD = fopen(file,"w")) == NULL) return(ERROR);  
	fprintf(FD,"%s0.0\nMASS %.6f\n",s,mass);
	fclose(FD);

	// writting mol file
	if (ho >= MAXATOM/2) {
		SM.molecule = M; SM.succ = NIL;
		strcpy(file,M->name); strcat(file,"_out"); sprintf(nb,"%d",Nsol);
		strcat(file,nb);
		moou_write(file,&SM,FALSE);
	}
	return(OK);
	}

LOCAL	t_integer number_bond(atom)
/**********************************************************
 count the number of bonding sites for a atom
 ***********************************************************/
atom_p	atom;
{
	set_bond_p	SB;
	t_integer	n = 0;
	if (atom == NIL) return(0);
	SB = atom->SB;
	while (SB) { n++; SB = SB->succ;}
	return(n);
}

LOCAL  t_pointer	saturated_atom(a,sat,h)
/**********************************************************************
 Return TRUE is a is in a saturated subgraph
 and there are still insaturated atoms
**********************************************************************/
atom_p		a;
t_bool		*sat;
t_integer	h;
{
	if (a->degre != number_bond(a)) *sat = FALSE;
	return((t_pointer)sat);
	}

LOCAL  t_bool		saturated_connectivity(A,a)
/**********************************************************************
 Return TRUE is a is in a saturated subgraph
 and there are still unsaturated atoms
**********************************************************************/
atom_p	*A,a;
{
	molecule_p	M = ((molecule_p)(((group_p)\
                            (a->group))->molecule));
	t_bool		sat = TRUE;
	t_integer	i;
	
	if (NBCC != 1) return(FALSE); // JLF 02/2012
	for (i = 0; i < M->size; i++) 
	    if (A[i]->degre != number_bond(A[i])) break;
	if (i == M->size) return(FALSE);
	dast_set_marq(M,NIL,NIL,FALSE); 
	dast_connectivity_atom(a,TRUE,NIL,saturated_atom,&sat,MAXATOM);
	dast_set_marq(M,NIL,NIL,FALSE);
	return(sat);
	}

LOCAL  t_err	enumerate(n,N,BOND,BS,A,SF,hi,scani,ho,scano)
/**********************************************************************
 Enumerate all structures matching BOND
	BOND is an array of all possible bonds
	BS is an array of active bonding sites
**********************************************************************/
t_integer		n,N,hi,ho;
t_bool			scani,scano;
bonding_site_p  *BOND;
t_bool			*BS;
atom_p			*A;
p_string		*SF;
{
	t_coord     ZEROC = {0,0,0};
	molecule_p	M = ((molecule_p)(((group_p)(A[0]->group))->molecule));
	t_integer	i,k,ncan = 0;
	p_string	sp = NIL,*CAN;
	t_err		status;

	if (Nenum++ > MAXSOL) { END_ENUMERATE = TRUE; Ntry = MAXSOL; return(ERROR); }
	if (n >= N) return(end_enumerate(A,SF,hi,scani,ho,scano));

	if (BS[n] == FALSE) return(enumerate(n+1,N,BOND,BS,A,SF,hi,scani,ho,scano));
	BS[n] = FALSE;
	for (k = 0; BOND[n][k].o > 0; k++) if (BS[BOND[n][k].index]) ncan++;
	CAN = (p_string *)malloc((ncan+1)*sizeof(p_string));
	for (k = 0; k <= ncan; k++) CAN[k] = NIL;
	for (i = 0; BOND[n][i].o > 0; i++) {
	  bond_p	b;
	  t_integer	j = BOND[n][i].index;
	  t_bool	o = BOND[n][i].o;
	  atom_p	ax = A[BOND[n][i].x];
	  atom_p	ay = A[BOND[n][i].y];
	  if (BS[j] == FALSE) continue;
	  if (dast_inq_bond(ax,ay)) continue;
	  BS[j] = FALSE;
	  b = create_bond(ax,ay,o); 
	  if (saturated_connectivity(A,ax) == FALSE) { 
	     p_string 	s = (p_string)compute_signature_atom_canonical(A[BOND[n][i].x],M->size+1,TRUE);
	     t_bool	dp = FALSE;
	     if (sp == NIL) dp = TRUE;
	     else if (strcmp(BOND[n][i].s,sp) != OK) dp = TRUE; /* different than previous */
	     for (k = 0; CAN[k]; k++) if (strcmp(s,CAN[k]) == OK) break;
	     /* recursivity if does not already exist with same initial signature */
	     if (CAN[k] == NIL ||  dp) {
		    /* printf("n = %d N = %d add bond [%d %d]\n",n,N,ax->ID,ay->ID); */
		    enumerate(n+1,N,BOND,BS,A,SF,hi,scani,ho,scano);
			}
	     if (CAN[k] == NIL) CAN[k] = strdup(s);
	     free(s);
	     }
	  delete_bond(b);
	  BS[j] = TRUE; sp = BOND[n][i].s;
	  if (END_ENUMERATE) break;
	  }

	for (k = 0; CAN[k]; k++) free(CAN[k]); free(CAN);
	BS[n] = TRUE;
	return(OK);
	}

LOCAL  t_err	enumerate_plusone(n,px,py,N,BOND,BS,A,SF,hi,scani,ho,scano)
/**********************************************************************
 Enumerate all structures matching BOND
 BOND is an array of all possible bonds
 BS is an array of active bonding sites
 The routine is on average 10 times faster than enumerate
 **********************************************************************/
t_integer		n,N,hi,ho;
t_bool			scani,scano;
bonding_site_p  *BOND;
t_bool			*BS;
atom_p			*A, px,py;
p_string		*SF;
{
	t_coord     ZEROC = {0,0,0};
	molecule_p	M = ((molecule_p)(((group_p)(A[0]->group))->molecule));
	t_integer	i,k,ncan = 0;
	p_string	sp = NIL,spp = NIL, spx = NIL,spy = NIL,spi = NIL, *CAN;
	t_err		status;
	
	if (Nenum++ > 100*MAXSOL) { END_ENUMERATE = TRUE; Ntry = MAXSOL; return(ERROR); }
	if (n >= N) return(end_enumerate(A,SF,hi,scani,ho,scano));
	
	if (BS[n] == FALSE) return(enumerate_plusone(n+1,px,py,N,BOND,BS,A,SF,hi,scani,ho,scano));
	BS[n] = FALSE;	
	for (k = 0; BOND[n][k].o > 0; k++) if (BS[BOND[n][k].index]) ncan++;
	CAN = (p_string *)malloc((ncan+1)*sizeof(p_string));
	for (k = 0; k <= ncan; k++) CAN[k] = NIL;
	for (i = 0; BOND[n][i].o > 0; i++) {
		bond_p	b;
		t_integer	j = BOND[n][i].index,cmp;
		t_bool	o = BOND[n][i].o;
		atom_p	ax = A[BOND[n][i].x];
		atom_p	ay = A[BOND[n][i].y];
		if (BS[j] == FALSE) continue;
		if (dast_inq_bond(ax,ay)) continue;
		BS[j] = FALSE;
		b = create_bond(ax,ay,o); 
		if (VERBOSE)  { printf("----------create bond [%d %d] o=%d n,i = %d %d\n",BOND[n][i].x,BOND[n][i].y,o,n,i);
		for (k = 0; k < M->size; k++) { set_bond_p	SBx = A[k]->SB; 
			while (SBx) { 
				atom_p	xy = (SBx->bond->atom1 == A[k])?SBx->bond->atom2:SBx->bond->atom1; 
				if (k<xy->ID) printf("[%d %d (%d)] ",k,xy->ID,SBx->bond->order); SBx = SBx->succ; } } 
		printf("\n"); }
	

		if (saturated_connectivity(A,ax) == FALSE) { 
			t_bool		can = FALSE;
			p_string 	sax = (p_string)compute_signature_atom_target(ax,SF);
			p_string 	say = (p_string)compute_signature_atom_target(ay,SF);
			p_string 	sa  = (p_string)malloc((strlen(sax)+strlen(say)+6)*sizeof(char));
			if (strcmp(sax,say) < 0) sprintf(sa,"%s-%1d-%s",say,o,sax);
			else sprintf(sa,"%s-%1d-%s",sax,o,say);
			for (k = 0; CAN[k]; k++) if (strcmp(sa,CAN[k]) == OK) break;
			if (CAN[k] == NIL) { can = TRUE; CAN[k] = strdup(sa); } else can = FALSE; 
			free(sax); free(say); free(sa);
//			sax = (p_string)compute_signature_atom_target_2(ax,SF);
//			if (px == NIL || px == ax) spx = sax; else spx = (p_string)compute_signature_atom_target_2(px,SF);
//			cmp = strcmp(sax,spx);  if (spx != sax) free(spx); free(sax);
			if (VERBOSE)  printf("----------[%d %d] n,i = %d %d can=%d cmp = %d\n",BOND[n][i].x,BOND[n][i].y,n,i,can,cmp);
			cmp=0;
			/* recursivity if does not already exist with same initial signature */
			if (cmp <= 0 && can) enumerate_plusone(n+1,ax,ay,N,BOND,BS,A,SF,hi,scani,ho,scano);
		}
		delete_bond(b);
		if (VERBOSE)  printf("----------delete bond [%d %d] n,i = %d %d\n",BOND[n][i].x,BOND[n][i].y,n,i);
		BS[j] = TRUE; 
		if (END_ENUMERATE) break;
	}
	for (k = 0; CAN[k]; k++) free(CAN[k]); free(CAN); 	
	BS[n] = TRUE;
	return(OK);
}

#define MAXSIGENUM	300
DEFINE	t_err	signature_enumeration(in_name,out_name,nbcc)
/******************************************************************************
 Enumerate all labeled molecule matching the signature in in_name
 ******************************************************************************/
t_bufstring	in_name,out_name;
t_integer	nbcc;
{
	EXTERN	t_real	dara_random_real();
	LOCAL	t_pointer		HEAP = NIL;
	bonding_site_p 		*BOND;
	t_bufstring		file;
	t_integer		Na = 0,Nb = 0,i,j,x,y,hi = -1,ho = -1;
	t_integer		t1,t2,Nc = 0;
	t_real			Nest = 1,nsol;
	atom_p			*AF,*A;
	p_string 		*SF,*SF1,*SB;
	t_bool			*BS;
	signature_p		signature;
	molecule_p		MF,M;
	group_p			GF,G;
	t_coord         ZEROC = {0,0,0};
	t_bool			scani = TRUE, scano = TRUE;
	
	NBCC = nbcc;
	if (strstr(in_name,"sig"))	scani = FALSE;
	if (strstr(out_name,"sig")) scano = FALSE;
	for (i = 0; i < strlen(out_name); i++) {
		if (out_name[i] >= '0' && out_name[i] <= '9') break;
	}
	if (i == strlen(out_name)) ho = MAXATOM;
	else ho = atoi(&out_name[i]);
	if (ODD(ho)) xerr_c("output signature diameter must be even");
	else ho = ho/2;
	strcpy(file,in_name);
	strcpy(INPUTFILE,in_name);
	signature = (signature_p)simo_read_signature_file(file,&hi);
	if (signature == NIL) xerr_c("no input signature");
	if (signature[1].value < 0) xerr_c("need more than one signature in file");
	if (hi < 2) xerr_c("input signature diameter must be greater than 2");
	if (ODD(hi)) xerr_c("input signature diameter must be even");
	else hi = hi/2;
	
	/* compute the number of atoms */
	for (i = 0; signature[i].as; i++) Na += signature[i].value;

	/* create file and solution molecules  */
	if ((HEAP = (t_pointer)heat_creer_hea(MAXHEAP)) == NIL) 
	   xerr_c("no more space");
	if (Na > MAXSIGENUM) 
	    xerr_c("maximum number of atoms exceeded");
	MF = (molecule_p)daal_create_molecule(HEAP,in_name,0,0,NIL,NIL,NIL,NIL,"",0);
	daal_create_group(MF->HEAP,MF->name,1,MF,0,NIL,NIL,MF->name,FALSE,NIL); 
	M  = (molecule_p)daal_create_molecule(HEAP,in_name,0,0,NIL,NIL,NIL,NIL,"",0);
	daal_create_group(M->HEAP,M->name,1,M,0,NIL,NIL,M->name,FALSE,NIL); 

	AF =  (atom_p *)malloc(Na*sizeof(atom_p));
	SF  = (p_string *)malloc(Na*sizeof(p_string));
	SF1 = (p_string *)malloc(Na*sizeof(p_string));
	for (x = 0; x < Na; x++) { AF[x] = NIL; SF[x] = NIL; }

	/* assign a signature to each atom
	  create signature SF and SF1 and SB for bonding sites
	 */
	x = 0; Nb = 0; 	SB = (p_string *)malloc((VALENCE+1)*Na*sizeof(p_string)); 
	for (i = 0; signature[i].as; i++) {
	  if (signature[i].as == NIL || signature[i].value < 1) break;
	  for (j = 0; j < (t_integer)signature[i].value; j++) {
		set_bond_p	SBx;
		simo_reset_molecule(MF); // each signature is a molecule
		AF[x] =  (atom_p)simo_signature_to_atom(MF,signature[i].as);
		AF[x]->degre = number_bond(AF[x]);
		if (/*AF[x]->degre < 1 ||*/ AF[x]->degre > VALENCE)
			xerr_c("atom with valence > 4");
		SF[x] =  (p_string)compute_signature_atom(AF[x],hi,scani);   
		SF1[x] = (p_string)compute_signature_atom(AF[x],hi-1,scani);
		if (VERBOSE) { printf("SF%d[%d]=%s SF%d[%d]=%s\n",hi,x,SF[x],hi-1,x,SF1[x]); }
		SBx = AF[x]->SB; while (SBx) {
			  atom_p	xy = (SBx->bond->atom1 == AF[x])?SBx->bond->atom2:SBx->bond->atom1;
			  SB[Nb] = (p_string)compute_signature_atom(xy,hi-1,scani);
			  if (VERBOSE) { printf("x=%d SB[%d]=%s\n",x,Nb,SB[Nb]); }
			  SBx = SBx->succ; Nb++;
		  }
	    if (++x > Na) xerr_c("in .scan file");
	    }
	  }

	/* create initial solution molecule */
	A =  (atom_p *)malloc(Na*sizeof(atom_p));
	for (x = 0; x < Na; x++) {
		// ZEROC.x = dara_random_real((t_real)100);
		// ZEROC.y = dara_random_real((t_real)100);
		// ZEROC.z = dara_random_real((t_real)100);
		A[x] = (atom_p)daal_create_atom(M->HEAP,AF[x]->name,x,\
		M->SG->group,NIL,NIL,NIL,ZEROC,AF[x]->potential_type,\
	    AF[x]->charge,AF[x]->degre,AF[x]->comment,FALSE);
	}
	M->size = Na;

	/* create array of possible bonds */
	BOND = (bonding_site_p *)malloc(Nb*sizeof(bonding_site_p));
	BS = (t_bool *)malloc(Nb*sizeof(t_bool));
	for (i = 0; i < Nb; i++) {
	    BOND[i] = (bonding_site_p)malloc(Nb*sizeof(bonding_site_t));
	    BS[i] = FALSE;
	    }
	i = 0;
	for (x = 0; x < Na; x++) {
	  t_integer	k = number_bond(AF[x]),l;
	  for (l = 0; l < k; l++) {
	    for (j = 0; j < Nb; j++) {
	        BOND[i][j].x = x; BOND[i][j].y = -1; 
	        BOND[i][j].index = -1; BOND[i][j].s = NIL; BOND[i][j].o = 0; 
	        } i++; 
	    }
	  }

	for (x = 0; x < Na; x++) {
	  set_bond_p	SBx = AF[x]->SB;
	  t_integer	I = index(AF,x);
	  i = 0; while (SBx) {
	    for (y = x+1; y < Na; y++) {
	      set_bond_p	SBy = AF[y]->SB; 
	      t_integer		J = index(AF,y);
	      if (strcmp(SF1[y],SB[I+i]) != OK) continue;  
		  if (VERBOSE) { printf("x=%d y=%d I+i=%d J=%d match SF1[y] = %s SB[I+i] = %s\n",x,y,I+i,J,SF1[y],SB[I+i]); }
	      j = 0; while (SBy) {
	        t_integer	order = 0;
	        if (strcmp(SF1[x],SB[J+j]) == OK) 
	           order = ((SBx->bond->order == SBy->bond->order)?SBx->bond->order:0);
			if (order < 1) {SBy = SBy->succ; j++; continue;}
			if (VERBOSE) { printf("x=%d y=%d I+i=%d J+j=%d match SF1[x] = %s SB[J+j] = %s\n\n",x,y,I+i,J+j,SF1[x],SB[J+j]); }
			if (BS[I+i] == -1 || BS[J+j] == -1) {SBy = SBy->succ; j++; continue;}
			if (dast_inq_bond(A[x],A[y])) {SBy = SBy->succ; j++; continue;}
	        /* create a bond if x or y is of degree 1 */ 
	        if (hi > 1 && (A[x]->degre < 2 || A[y]->degre < 2)) {
				create_bond(A[x],A[y],order);
				if (VERBOSE) printf("bond %d %d created BS %d %d\n",x,y,I+i,J+j);
				BS[I+i] = -1; BS[J+j] = -1;
				SBy = SBy->succ; j++; continue;
			}
			if (strcmp(A[x]->potential_type,"H") == OK || 
				strcmp(A[x]->potential_type,"F") == OK ||
				strcmp(A[x]->potential_type,"Cl") == OK) {
				create_bond(A[x],A[y],order);
				if (VERBOSE) printf("bond %d %d created BS %d %d\n",x,y,I+i,J+j);
				BS[I+i] = -1; BS[J+j] = -1;
				SBy = SBy->succ; j++; continue;
			}  
			if (strcmp(A[y]->potential_type,"H") == OK || 
				strcmp(A[y]->potential_type,"F") == OK ||
				strcmp(A[y]->potential_type,"Cl") == OK) {
				create_bond(A[x],A[y],order);
				if (VERBOSE) printf("bond %d %d created BS %d %d\n",x,y,I+i,J+j);
				BS[I+i] = -1; BS[J+j] = -1;
				SBy = SBy->succ; j++; continue;
			}
			BS[I+i] = TRUE; BS[J+j] = TRUE;
			BOND[I+i][J+j].index = J+j;
			BOND[I+i][J+j].y = y; 
			BOND[I+i][J+j].s = malloc((strlen(SF[y])+4)*sizeof(char));
			strcpy(BOND[I+i][J+j].s,"");
			sprintf(BOND[I+i][J+j].s,"%1d %s",order,SF[y]);
			BOND[I+i][J+j].o = order;
			if (VERBOSE) { printf("----BOND[%d][%d] index = %d y = %d s = %s\n",I+i,J+j,BOND[I+i][J+j].index,BOND[I+i][J+j].y,BOND[I+i][J+j].s); }
	        SBy = SBy->succ; j++;
	        } }
	    SBx = SBx->succ; i++;
	    } }

	/* sort BOND according to signature */
	for (i = 0; i < Nb; i++)
		qsort(&BOND[i][0],Nb,sizeof(bonding_site_t),compare_bonding_site);

	/* check and initial count */
	for (i = 0; i < Nb; i++)
	  if (BS[i] == FALSE) {
	     printf("ERROR: atom %d: %s bonding site %d: %s cannot be linked\n",\
                        BOND[i][0].x,SF[BOND[i][0].x],i,SB[i]);
	     xerr_c("signature inconsistent");
	     }
	  else if (BS[i] == -1) BS[i] = FALSE;

	for (i = 0; i < Nb; i++) {
		t_integer	n = 0;
		if (BS[i] == FALSE) continue;
		for (j = 0; j < Nb; j++) if (BOND[i][j].o > 0) n++;
		if (n < 1) continue;
		if (VERBOSE) printf("bs %d [%d] possibilities %d \n",i,BOND[i][0].x,n); 
		Nest = Nest * n;
	}
	
	/* free some memory */
	for (i = 0; i < Na; i++) free(SF1[i]);
	free(SF1); free(AF);

	/* enumeration */
	END_ENUMERATE = FALSE; 
	t1 = clock();
	if (ho <= hi+1 && scano == FALSE) 
		enumerate_plusone(0,NIL,NIL,Nb,BOND,BS,A,SF,hi,scani,ho,scano);
	else 
		enumerate(0,Nb,BOND,BS,A,SF,hi,scani,ho,scano);
	t2 = clock();
	if (Ntry >= MAXSOL && HIT == 0) nsol = Nest;
	else nsol = Nsol;
	for (x = 0; x < Na; x++) if (A[x]->potential_type[0] != 'H') Nc++;
	fprintf(stdout,"File = %s out = %s dim-out = %d size = %d %d dim-in = %d Nb = %d Nsol = %d NBCC = %d Ntry/Nsol= %.1f HIT = %d nsol = %E CPU-time = %.6f\n",\
						  in_name,out_name,2*ho,M->size,Nc,2*hi,Nb,Nsol,NBCC,(float)Ntry/(float)Nsol,HIT,nsol,(float)(t2-t1)/(float)1000000L);
	return(OK);
	}







