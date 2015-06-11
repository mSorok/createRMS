/***********************************************************
 Routines to canonize DAG (Directed Acyclic Graph)
 using signature 
 Jean-Loup Faulon & Pablo Carbonell, February 2012
 The code was modified to 
 1) Handle simple signature without ring closure (here the
 signature is an extended conectivity and is identical
 to the tree used when developping the CIP rules
 2) Treat the CIP rules R/S and Z/E
 The code does not deal with 
    - relative configuration (R* and S*)
	- pseudo asymmetric center (r and s)
 Beware before changing invariant calculation any increment must
 be greater than e-7
 ***********************************************************/
#include <general.h>
#include <eps.h>
#include "signature.h"
#include <time.h>

extern t_err compute_cip();

#define ENDINVLIST -1

typedef	char	element_t[MAXSTRING/10];

typedef	struct	VERTEX {
	atom_p			atom;
	element_t		element;
	// invlist is populated by arrays of 5 items for
	//	the vertex and its neigbors 
	// Each array comprise
	// - Z number + stereo
	// - maximum bond order to the parent
	// - charge
	// - invariant CIP 
	// - bond order ij where i is the 
	//       node and j the clid/parent,
	// invlist terminates with an array of n+1 node invariant
	// where we find in decreasing order the invariant of the n neigbors
	// followed by the invariant of the node itself
	// inlist terminates with the charactere ENDINVLIST
	t_integer *		invlist;
	t_integer		invariant;
	struct VERTEX * *child;
	struct VERTEX * *parent;
} vertex_t, *vertex_p;

typedef	struct	EDGE { /* [x,y] edge */
	t_integer	*y;
	t_integer	*l;
} edge_t, *edge_p;


// a special structure for compute invariant occurance
typedef	struct	INOCCSTR {
	t_integer	inv,occ,ID;
} inocc_t, *inocc_p;

LOCAL	t_bool		SCAN = TRUE; /* SCAN mode */
LOCAL	t_bool		CIP = TRUE; /* CIP mode */
LOCAL	t_bool		SIG = TRUE; /* SIG mode */
LOCAL	t_integer *	INVAR = NIL;
LOCAL	t_integer **LAYER = NIL;
LOCAL	t_integer *	LABEL = NIL; /* for printing only */
LOCAL	t_integer *	OCCUR = NIL;
LOCAL	t_integer *	COLOR = NIL;
LOCAL	p_string	SMAX = NIL;
LOCAL	t_integer	NSTR = 0; /* number of strings */
LOCAL	t_integer	NITR = 0; /* max ITER */
LOCAL	t_integer	SIZE = 0;
LOCAL	t_integer	LL	= -1; /* used for ring closer */

/* the following variables
	are used when the structure
	is relabeled i.e. canonizing or finding
    the orbits */
LOCAL	t_integer *	LACAN = NIL;
LOCAL	t_integer 	NBCUR = 0;
LOCAL	t_integer *	LACUR = NIL;

/*----------------------------------------------------------
 COMPUTE EDGE and VERTEX invariant
 -----------------------------------------------------------*/

LOCAL   t_integer	order(v1,v2)
/***********************************************************
 return the order of the directed bond v1<-v2
 ***********************************************************/
vertex_p	v1,v2;
{
	bond_p		b;
	t_integer	order = -1;
	
	if ((v1 == NIL) || (v2 == NIL)) return(order);
	b = (bond_p)dast_inq_bond_directed(v1->atom,v2->atom); 
	if (b) order = b->order;
	
	if (v1->atom->potential_type[strlen(v1->atom->potential_type)-1] == 'p' &&
		v2->atom->potential_type[strlen(v2->atom->potential_type)-1] == 'p') 
		order = 4;
	
	return(order);
}

LOCAL   t_integer	edge_invariant(v1,v2)
/***********************************************************
 return the invariant of the directed bond v1<-v2
 10 single
 15 aromatic
 20 double
 24 seqtrans \
 28 seqcis  / seqcis > seqtrans CIP rule 3
 30 triple
 ***********************************************************/
vertex_p	v1,v2;
{
	bond_p			b;
	t_integer		inv = order(v1,v2)*10;
	if (inv < 0) return(0);
	if (inv == 40) inv = 15;
	if (inv != 20) return(inv);
	

	if (v1->atom->potential_type[strlen(v1->atom->potential_type)-1] == '\\' &&
		v2->atom->potential_type[strlen(v2->atom->potential_type)-1] == '\\') 
		inv = 24;
	if (v1->atom->potential_type[strlen(v1->atom->potential_type)-1] == '/' &&
		v2->atom->potential_type[strlen(v2->atom->potential_type)-1] == '/') 
		inv = 28;
	return(inv);
}

LOCAL	t_integer	Z_vertex(v)
/**********************************************************
 Return the Z number of a vertex
 ***********************************************************/
vertex_p	v;
{
	EXTERN		t_real	dast_Z_element_atom();
	if (v == NIL) return(0);
	if (v->atom == NIL) return(0);
	return((t_integer)dast_Z_element_atom(v->atom));
}

LOCAL	t_integer	Stereo_vertex(v)
/**********************************************************
Precedence of stereocenters
IUPAC CIP Rule 5. R has priority over S that has priority over non-chiral
 ***********************************************************/
vertex_p		v;
{
	t_name		element;
	t_integer		Z;
	t_integer	i;

	if (v == NIL) return(0);
	if (v->atom == NIL) return(0);

	strcpy(element,v->atom->potential_type);
	for (i = 0; element[i]; i++)
		if ((element[i] < 'A') || (element[i] > 'z')) break;
	if (element[i] != '\0') {
		if (element[i++] == '@')
			if (element[i] == '@') //R
				Z = 6;
			else Z= 3; //S
		else Z=0;
	}
	else Z=0;
	return(Z);
}


LOCAL   t_integer	charge_invariant(v)
/***********************************************************
 1 +3, 2 +2, 3 +1, 4 :, 5 -1, 6 -2, 7 -4
 It has to be added with the lowest priority (i.e. only is used if branches
 are the same)
 ***********************************************************/
vertex_p	v;
{
	if (v == NIL) return(0);
	if (v->atom == NIL) return(0);	
	return((t_integer)v->atom->charge);
}

LOCAL	t_integer	sig_parent(n,l,P)
/***********************************************************
 Search if n has a parent *P having the same atom
 l is the layer of n
 The routine return *P and the layer of *P
 ***********************************************************/
vertex_p	n,*P;
t_integer	l; 
{
	vertex_p	m = NIL;
	t_integer	lm = l+1,i;
	element_t	element;
	
	return(0); // PC: Cip invariant is now computed in the CIP section!!!
	if (n == NIL) return(0);
	if (CIP == FALSE) return(0); // only in CIP mode
	if (n->parent == NIL) return(0);
	for (m = n->parent[0]; m; m = m->parent[0]) {
		if (m->atom == n->atom) break;
		if (m->parent == NIL) break;
		lm++;
	}
	if (m == NIL) return(0);
	if (m->atom != n->atom) return(0);
	*P = m;
	return(lm);
}

LOCAL	t_integer	cip_invariant(n,l)
/***********************************************************
 To follow the CIP rule 1.b one check if vertex n
 does not already occur in the list of its parents
 Let lm be the layer of the parent m
 First a label is added to both n and m
 Second lm is added (returned) to the the invariant of n 
 ***********************************************************/
vertex_p	n;
t_integer	l; 
{
	vertex_p	m;
	return(sig_parent(n,l,&m));
}

LOCAL	void	vertex_invariant(v,l)
/***********************************************************
  Populate invlist for vertex v in layer l
  - Z number including stereo
  - the max order of the bonds to its parents
  - the charge
  - the CIP invariant (in CIP mode)
  - 0
 ***********************************************************/
vertex_p	v;
t_integer	l;
{
	t_integer	size = 0, i = 0;
	t_integer	max = 0, m;
	if (v == NIL) return;
	if (v->invlist) return;
	
	// first allocate invlist 
	if (v->parent) for (i = 0;v->parent[i];i++); 
	size = i; i = 0;
	if (v->child) for (i = 0;v->child[i];i++); 
	if (i > size) size = i;
	v->invlist = (t_integer *)calloc((6*(size+1)+1),sizeof(t_integer));

	// populate invlist for the vertex 
	v->invlist[0] = 10*Z_vertex(v)+Stereo_vertex(v);
		if (v->parent)
			for (i = 0; v->parent[i]; i++) {
				if ((m = edge_invariant(v,v->parent[i])) > max) max = m;
		}
	v->invlist[1] = max;
	v->invlist[2] = charge_invariant(v);
	if (CIP) v->invlist[3] = cip_invariant(v,l); // PC: not in use, cip invariants computed in the CIP section
	// 4 is empty
}

/*----------------------------------------------------------
 COMPARE INVARIANTS
 ALL SORT ARE IN INCREASING ORDER
 -----------------------------------------------------------*/

LOCAL   int	cmp_invariant(i1,i2)
/***********************************************************
 i1 < i2  -1 ; i1 > i2 +1 ; i1 = i2 0
 ***********************************************************/
t_integer	*i1,*i2;
{
	if ((i1 == NIL) || (i2 == NIL)) return(0);
	if (i1 == NIL) return(-1);
	if (i2 == NIL) return( 1);
	if (*i1 < *i2) return(-1);
	if (*i1 > *i2) return( 1);
	return(0);
}

LOCAL   int	cmp_6_invariant(i1,i2)
/***********************************************************
 i1[] < i2[]  -1 ; i1[] > i2[] +1 ; i1[] = i2[] 0
 ***********************************************************/
t_integer	i1[],i2[];
{
	int	r,j;
	if ((i1 == NIL) || (i2 == NIL)) return(0);
	if (i1 == NIL) return(-1);
	if (i2 == NIL) return( 1);
	for (j = 0; j < 6; j++)
		if ((r=cmp_invariant(&i1[j],&i2[j]))) return(r);
	return(0);
}


LOCAL   int	cmp_vertex_invariant(n1,n2)
/***********************************************************
 i1 < i2 -1 ; i1 > i2  1 ; i1 = i2 0
 ***********************************************************/
vertex_p	*n1,*n2;
{
	t_integer	r;
	if ((n1 == NIL) || (n2 == NIL))   return(0);
	if ((*n1 == NIL) && (*n2 == NIL)) return(0);
	if (*n1 == NIL) return(-1);
	if (*n2 == NIL) return( 1);
	return(cmp_invariant(&(*n1)->invariant,&(*n2)->invariant));
}

LOCAL   int	cmp_vertex_invlist(n1,n2)
/***********************************************************
 i1 < i2 -1 ; i1 > i2  1 ; i1 = i2 0
 ***********************************************************/
vertex_p	*n1,*n2;
{
	t_integer	r,j;
	if ((n1 == NIL) || (n2 == NIL))   return(0);
	if ((*n1 == NIL) && (*n2 == NIL)) return(0);
	if (*n1 == NIL) return(-1);
	if (*n2 == NIL) return( 1);
		for (j = 0; (*n1)->invlist[j] != ENDINVLIST && (*n2)->invlist[j] != ENDINVLIST ; j++)
			if ((r=cmp_invariant(&(*n1)->invlist[j],&(*n2)->invlist[j]))) return(r);  
		if ((*n1)->invlist[j] == ENDINVLIST) // PC Vertex1 is at the end
			if ((*n2)->invlist[j] == ENDINVLIST) // PC Vertex2 is at the end
				return (0); // PC Both vertices have same number of children, return who is heavier
			else return(-1); // PC Vertex1 is smaller
		else return(1); // PC Vertex2 is smaller
}

LOCAL   int		cmp_compute_invariant_occurence(inocc1,inocc2)
/**********************************************************************
 i1 < i2  -1 ; i1 > i2 +1 ; i1 = i2 0
 This routine is called by compute_invariant_occurence
 and used the special data str inocc_p work like the previous one
 **********************************************************************/
inocc_p	inocc1,inocc2;
{
	int r,j;
	if ((inocc1 == NIL) || (inocc2 == NIL)) return(0);
	if (inocc1 == NIL) return(-1);
	if (inocc2 == NIL) return( 1);
	for (j = 0; LAYER[inocc1->ID][j] != ENDINVLIST && LAYER[inocc2->ID][j] != ENDINVLIST ; j++)
		if ((r=cmp_invariant(&(LAYER[inocc1->ID][j]),&(LAYER[inocc2->ID][j])))) return(r);
	return(0);
}

/*----------------------------------------------------------
  UTILITIES
-----------------------------------------------------------*/

DEFINE	p_string sicd_signature_reset()
/***********************************************************
 Reset all LOCAL variables
 ***********************************************************/
{
	t_integer	i;
	if (LABEL) free(LABEL); LABEL = NIL; 
	if (INVAR) free(INVAR); INVAR = NIL;
	if (LAYER) {
		for (i = 0; i < SIZE; i++) free(LAYER[i]);
		free(LAYER); LAYER = NIL;
	}
	if (OCCUR) free(OCCUR); OCCUR = NIL;
	if (COLOR) free(COLOR); COLOR = NIL;
	SMAX = NIL; NSTR = 0;  NITR = 0; SIZE = 0;
	if (LACAN) free(LACAN); LACAN = NIL;
	NBCUR = 0;
	if (LACUR) free(LACUR); LACUR = NIL;
	return(OK);
}

LOCAL	t_err		print_layer(L,h,INV)
/***********************************************************
 ***********************************************************/
vertex_p *	L[];
t_integer	h;
t_integer *	INV;
{
	t_integer	l,i,j;
	vertex_p *	N;
	
	if (L == NIL) return(ERROR);
	for (l = h; l >= 0; l--) if ((N = (vertex_p *)L[l])) {
		for (i = 0; N[i]; i++) {
			vertex_p	vertex = N[i];
			if (INV) printf("layer = %d vertex = %d %s inv = %d (%d) : ",\
				   l,vertex->atom->ID+1, vertex->element, vertex->invariant,INV[vertex->atom->ID]);
			else printf("layer = %d vertex = %d %s inv = %d : ",\
						l,vertex->atom->ID+1, vertex->element, vertex->invariant);
			for (j = 0; vertex->child[j]; j++)
                printf(" %d %s %d order=%d",\
					   vertex->child[j]->atom->ID+1,vertex->child[j]->element,\
					   vertex->child[j]->invariant,order(vertex,vertex->child[j]));
			printf("\n");
		}
	}
	
	for (l = h; l >= 0; l--) {
		printf("l = %d: ",l);
		if ((N = (vertex_p *)L[l])) for (i = 0; N[i]; i++) {
			vertex_p	vertex = N[i];
			printf("%d ",vertex->atom->ID+1);
		}
		printf("\n");
	}
	
	for (l = h; l >= 0; l--) {
			if ((N = (vertex_p *)L[l]))
				for (i = 0; N[i]; i++) {
					vertex_p	vertex = N[i];
					printf("%d Atom %d: \t",l,vertex->atom->ID+1);
					t_integer *invlist = vertex->invlist;
					int c = 0;
					while (*invlist != ENDINVLIST) {
						if (!((c++)%5))  printf("| ");
						printf("%2d ",*invlist++);
					}
					printf("| ");
					printf("\n");
			}
	}

	if (SCAN) {

		for (i=0; i < SIZE; i++) {
			printf("ATOM %d [",i+1);
			for (l = 0; l <= h; l++) {
				printf(" %d",LAYER[i][l]);
			}
			printf(" ]\n");
		}
	}
	return(OK);
}

LOCAL	t_err		delete_dag(L,h)
/***********************************************************
 Delete L[] data str
***********************************************************/
vertex_p *	L[];
t_integer	h;
{
	t_integer	l,i;
	vertex_p *	N;

	for (l = h; l >= 0; l--) if ((N = (vertex_p *)L[l])) {
	  for (i = 0; N[i]; i++) {
	      vertex_p	vertex = N[i];
 	      free(vertex->invlist); 
 	      free(vertex->parent); 
		  free(vertex->child);
	      free(vertex);
	      }
	  free(N);
	  }
	free(L);
	return(OK);
	}

/*----------------------------------------------------------
  TREE BUILDER
-----------------------------------------------------------*/

LOCAL	t_bool	inq_edge(E,x,y,l)
/***********************************************************
 TRUE if edge already created at previous layer
***********************************************************/
atom_p		x,y;
t_integer	l;
edge_p		E;
{
	t_integer	i;

	if (E[x->ID].y == NIL) return(FALSE);
	for (i = 0; E[x->ID].y[i] != -1; i++) 
	    if (E[x->ID].y[i] == y->ID) {
	       if (l < 0) return(TRUE);
               if (E[x->ID].l[i] > l) return(TRUE);
	       }
	return(FALSE);
	}

LOCAL	t_err	set_edge(E,x,y,l)
/***********************************************************
 insert edge [x,y] in E
***********************************************************/
atom_p		x,y;
t_integer	l;
edge_p		E;
{
	t_integer	i;

	if (E[x->ID].y == NIL) {
	   E[x->ID].y = (t_integer *)\
            malloc((x->degre+1)*sizeof(t_integer));
	   E[x->ID].l = (t_integer *)\
            malloc((x->degre+1)*sizeof(t_integer));
	   E[x->ID].y[0] = -1;
	   }
	if (E[y->ID].y == NIL) {
	   E[y->ID].y = (t_integer *)\
            malloc((y->degre+1)*sizeof(t_integer));
	   E[y->ID].l = (t_integer *)\
            malloc((y->degre+1)*sizeof(t_integer));
	   E[y->ID].y[0] = -1;
	   }
	for (i = 0; E[x->ID].y[i] != -1; i++) 
	    if (E[x->ID].y[i] == y->ID) return(OK);
	E[x->ID].y[i] = y->ID; E[x->ID].l[i] = l; E[x->ID].y[i+1] = -1; 

	for (i = 0; E[y->ID].y[i] != -1; i++) 
	    if (E[y->ID].y[i] == x->ID) return(OK);
	E[y->ID].y[i] = x->ID; E[y->ID].l[i] = l; E[y->ID].y[i+1] = -1; 
	return(OK);
	}

LOCAL	t_err	add_child(v,c)
/***********************************************************
 add child c to v
***********************************************************/
vertex_p	v,c;
{
	t_integer	i;
	/* create children */
	if (v->child == NIL) {
	   atom_p	a = v->atom;
	   set_bond_p	SB;
	   t_integer	n = 0;
	   SB = (set_bond_p)a->SB;
	   while (SB) { n++; SB = SB->succ;}
	   v->child = (vertex_p *)malloc((n+1)*sizeof(vertex_p));
	   for (i = 0; i < n+1; i++) v->child[i] = NIL;
	   }
	if (c == NIL) return(OK);

	/* add child */
	for (i = 0; v->child[i]; i++);
	v->child[i] = c;
	return(OK);
	}

LOCAL	t_err	add_parent(v,p)
/***********************************************************
 add parent p to v
***********************************************************/
vertex_p	v,p;
{
	t_integer	i;
	/* create children */
	if (v->parent == NIL) {
	   atom_p	a = v->atom;
	   set_bond_p	SB;
	   t_integer	n = 0;
	   SB = (set_bond_p)a->SB;
	   while (SB) { n++; SB = SB->succ;}
	   v->parent = (vertex_p *)malloc((n+1)*sizeof(vertex_p));
	   for (i = 0; i < n+1; i++) v->parent[i] = NIL;
	   }
	if (p == NIL) return(OK);

	/* add parent */
	for (i = 0; v->parent[i]; i++);
	v->parent[i] = p;
	return(OK);
	}

LOCAL	t_err	add_vertex(n,a,E,l,N)
/***********************************************************
 Create a new vertex in N corresponding to atom a
 if a is not already in the layer.
***********************************************************/
vertex_p	n,N[];
atom_p		a;
t_integer	l;
edge_p		E;
{
	t_integer	i;

	/* check if a is a parent of n  */
	if (SCAN) {
		if (inq_edge(E,n->atom,a,l) == TRUE) return(ERROR);
	} else if (n) if (n->parent) { // both CIP and SIG modes
		vertex_p	m = n->parent[0];
		if (m->atom == a) return(ERROR); // a - n - a
		// CIP mode only 
		if (CIP) for (m = n->parent[0]; m; m = m->parent[0]) {
			if (m->atom == n->atom) return(ERROR); 
				// cycle n-n closed no need to add atom a
			if (m->parent == NIL) break;
		}
	}
	/* check if a is already in layer */
	if (SCAN) for (i = 0; N[i]; i++) {
	  if (N[i]->atom == a) break;
	  }
	else for (i = 0; N[i]; i++);
	
	if (N[i] == NIL) {
		N[i] = (vertex_p)malloc(sizeof(vertex_t));
		N[i]->atom = a; 
		strcpy(N[i]->element,"[");
		strcat(N[i]->element,a->potential_type);
		strcat(N[i]->element,"]");
	    N[i]->invariant = 1; N[i]->invlist = NIL;
		N[i]->child = N[i]->parent = NIL; 
	    add_child(N[i],NIL); add_parent(N[i],NIL);
	    N[i+1] = NIL;
	   }
	/* add vertex N[i] to child of n */
	add_child(n,N[i]);
	/* add vertex n to parent of N[i]  */
	add_parent(N[i],n);
	
	/* update edge array */
	if (SCAN) return(set_edge(E,n->atom,a,l));
	return(OK);
	}

LOCAL	t_err	build_layer(N,E,L,h,a1,a2,ha12)
/***********************************************************
 Build a signature tree layer by layer
 of efficient space size of height h, 
 N[] is the list of vertex in the previous layer
 E[] is a list of edges already inserted
***********************************************************/
vertex_p	N[];
t_integer	h,ha12;
edge_p		E;
vertex_p * 	L[];
atom_p		a1,a2; // the roots
{
	vertex_p *	NN;
	t_integer	i,nb = 0;
	
	if (h < 0) return(OK);	
	if (N[0] == NIL) return(OK);
	
	/* compute the number of vertexs to be created 
           in the current layer */
	for (i = 0; N[i]; i++) {
	  atom_p	a = N[i]->atom,aa;
	  set_bond_p	SB;  
	  if (a == NIL) break;
	  SB = (set_bond_p)a->SB;
	  while (SB) {
	    bond_p	b = SB->bond;
	    aa = ((b->atom1 == a)?b->atom2:b->atom1);
		if (aa == NIL) { SB = SB->succ; continue; }
		if (h == ha12-1) if (aa == a1 || aa == a2) { SB = SB->succ; continue; }
		if (SCAN) { if (inq_edge(E,a,aa,h) == FALSE) nb++; } // JLF Feb. 2011
		else nb++;
	    SB = SB->succ;
	    }
	  }

	/* create vertex for current layer  */
	NN = (vertex_p *)malloc((nb+2)*sizeof(vertex_p)); NN[0] = NIL; 
	for (i = 0; N[i]; i++) {
	  vertex_p	n = N[i];
	  atom_p	a = N[i]->atom,aa;
	  set_bond_p	SB;  
	  if (a == NIL) break;
	  SB = (set_bond_p)a->SB;
	  while (SB) {
	    bond_p	b = SB->bond;
	    aa = ((b->atom1 == a)?b->atom2:b->atom1);
		if (aa == NIL) { SB = SB->succ; continue; }
		if (h == ha12-1) if (aa == a1 || aa == a2) { SB = SB->succ; continue; }
		add_vertex(n,aa,E,h,NN);
	    SB = SB->succ;
	    }
	  }
	L[h] = &NN[0];
	return(build_layer(NN,E,L,h-1,a1,a2,ha12));
	}

LOCAL	t_err	build_dag(a1,a2,L,h)
/***********************************************************
 Build a dag (signature tree) of height h for the root
 atoms a1 and a2.
***********************************************************/
atom_p		a1,a2;
t_integer	h;
vertex_p * 	L[];
{
	edge_p		E = NIL; /* array of directed edges */
	vertex_p	root1 = NIL,root2 = NIL;
	vertex_p	*N;
	t_integer	i,j;
	t_err		err;

	if (a1 == NIL && a2 == NIL) return(NIL);
	if (h < 0) return(NIL);
	if (a1) {
		root1 = (vertex_p)malloc(sizeof(vertex_t));
		root1->atom = a1; 
		strcpy(root1->element,"[");
		strcat(root1->element,a1->potential_type);
		strcat(root1->element,"]");
		root1->invariant = 1; root1->invlist = NIL;
		root1->child = root1->parent = NIL; 
		add_child(root1,NIL);
	}
	if (a2) {
		root2 = (vertex_p)malloc(sizeof(vertex_t));
		root2->atom = a2; 
		strcpy(root2->element,"[");
		strcat(root2->element,a2->potential_type);
		strcat(root2->element,"]");
		root2->invariant = 1; root2->invlist = NIL;
		root2->child = root2->parent = NIL; 
		add_child(root2,NIL);
	}

	N = (vertex_p *)malloc(3*sizeof(vertex_p));
	N[0] = root1; N[1] = root2; N[2] = NIL;
	L[h] = &N[0]; 
	if (h < 1) return(OK);
	
	if (SCAN == FALSE) return(build_layer(N,NIL,L,h-1,a1,a2,h));

	E = (edge_p)malloc((SIZE+1)*sizeof(edge_t));
	for (i = 0; i < SIZE; i++) E[i].y = E[i].l = NIL;
	err = build_layer(N,E,L,h-1,a1,a2,h);

	/* delete edges */
	for (i = 0; i < SIZE; i++) {free(E[i].y); free(E[i].l);}
	free(E); 

	return(err);
	}

/*----------------------------------------------------------
  PRINT SIGNATURES IN SCAN MODE
-----------------------------------------------------------*/
LOCAL	t_err	reset_print_string(L,h)
/***********************************************************
 Remove all element labels
 ***********************************************************/
vertex_p *	L[];
t_integer	h;
{
	vertex_p	*N;
	t_integer	n = 0,i,j,l;
	
	N = L[h]; if (N[0] == NIL) return(OK);
	
	for (l = h; l >= 0; l--) if ((N = (vertex_p *)L[l])) {
		for (i = 0; N[i]; i++) {
			vertex_p	r = N[i];
			for (j = 0; j < strlen(r->element); j++) 
				if (r->element[j] == ',') {
					r->element[j] = ']';
					r->element[j+1] = '\0';
					break;
				}
		}
	}
	return(OK);
}

LOCAL	t_err	print_bond(s,p,r)
/***********************************************************
 Print tree r in string s
 p is the parent r is the current vertex
 ***********************************************************/
vertex_p	p,r;
p_string	s;
{
	t_integer	i = 0,k,ii,oder;
	t_bufstring	tagAtom;
		
	if (r == NIL || p == NIL) return(OK);
	if (order(p,r) < 2) return(OK);  
	if (order(p,r) == 2) { 
		if (r->element[strlen(r->element)-2] == '/' && r->element[strlen(r->element)-2] == '/') 
			 strcat(&s[0],"\\=/"); // seqcis
		else if (r->element[strlen(r->element)-2] == '\\' && r->element[strlen(r->element)-2] == '\\' ) 
			 strcat(&s[0],"\\=\\"); // seqtrans
		else strcat(&s[0],"=");
		}
	else if (order(p,r) == 3) strcat(&s[0],"t");
	else if (order(p,r) == 4) strcat(&s[0],"p"); /*  aromatics JLF */
	else sprintf(&s[0],"-%d-",order(p,r));
	return(OK);
}


LOCAL	t_integer	print_string(s,p,r,E,LAB,OCC)
/***********************************************************
 Print tree r in string s
  p is the parent r is the current vertex
***********************************************************/
vertex_p	p,r;
p_string	s;
edge_p		E;
t_integer	LAB[],OCC[];
{
	t_integer	i = 0,j = 0, k,ii,oder;	
	element_t	element;
	if (r == NIL) return(0);

	/* add label for atoms occuring more than once */
	if (OCC[r->atom->ID] >  1) {
	   strcpy(element,r->element);
	   for (i = 0; i < strlen(element); i++) 
	       if (element[i] == ',') break;
	   if (i == strlen(element)) {
		  if (LAB[r->atom->ID] <  0) {
			   LAB[r->atom->ID] = (++LL);
		  }
	      strcpy(r->element,"[");
	      strcat(r->element,r->atom->potential_type);
	      strcpy(element,""); sprintf(element,"%d",LAB[r->atom->ID]);
	      strcat(r->element,",");
	      strcat(r->element,element);
	      strcat(r->element,"]");
	      }
	   }
	print_bond(s,p,r); 
	// remove / and \ signs
	for (i = 0; i < strlen(r->element); i++) 
		if ((r->element[i] != '/') && (r->element[i] != '\\'))
			element[j++] = r->element[i];
	element[j++] = '\0';
	strcat(s,element);	
	if (LACUR[r->atom->ID] < 0) LACUR[r->atom->ID] = NBCUR++;
	i = strlen(&s[0]);
	if (r->child == NIL) return(i);
	for (k = 0; r->child[k]; k++); if (k < 1) return(i);

	/* recursion */
	ii = i; sprintf(&s[i++],"("); 
	for (k = 0; r->child[k]; k++) {
	  if (inq_edge(E,r->atom,\
					 r->child[k]->atom,-1)) continue;
	  set_edge(E,r->atom,r->child[k]->atom,1);
	  i += print_string(&s[i],r,r->child[k],E,LAB,OCC); 
	  }
	sprintf(&s[i++],")");
	if (i == ii+2) i = ii;
	return(i);
	}

LOCAL	t_err		occur_string(r,E,OCC)
/***********************************************************
 Count the number of occurence for vertex r
***********************************************************/
vertex_p	r;
edge_p		E;
t_integer	OCC[];
{
EXTERN	int 		cmp_vertex_invariant();
	t_integer	k;

	if (r == NIL) return(OK);
	if (OCCUR[r->atom->ID] >  1) OCC[r->atom->ID] += 1;

	/* sort children by invariant */
	for (k = 0; r->child[k]; k++);
	if (k == 0) return(OK);
	qsort(r->child,k,sizeof(vertex_p),cmp_vertex_invariant);

	/* recursion */
	for (k = 0; r->child[k]; k++) {
	  if (inq_edge(E,r->atom,\
					 r->child[k]->atom,-1)) continue;
	  set_edge(E,r->atom,r->child[k]->atom,1);
	  occur_string(r->child[k],E,OCC); 
	  }
	return(OK);
	}

LOCAL	t_integer	delete_edge(r,E)
/***********************************************************
 Delete the edges of r
***********************************************************/
vertex_p	r;
edge_p		E;
{
	t_integer	k;

	if (r == NIL) return(OK);
	for (k = 0; r->child[k]; k++) 
	  set_edge(E,r->atom,r->child[k]->atom,-1);
	return(OK);
	}

LOCAL	t_integer	layer_size(L,h)
/***********************************************************
 Return the number of vertices in all layers
***********************************************************/
vertex_p *	L[];
t_integer	h;
{
	t_integer	i,l,n=0;
	vertex_p	*N;
	for (l = 0; l <= h; l++) {
	  if ((N = L[l])) {
             for (i = 0; N[i]; i++) n++;
	     }
	  }
	return(n);
	}

LOCAL	p_string	layer_print_string(L,h,LAB)
/***********************************************************
 Print a signature
 ***********************************************************/
vertex_p *	L[];
t_integer	h,LAB[];
{
	edge_p		E = NIL;
	t_integer *	OCC = NIL;
	vertex_p	*N;
	p_string	s;
	t_integer	n = 0,i,j,l;
	
	N = L[h]; if (N[0] == NIL) return(NIL);
	
	/* OCCUR is reset here to its initial value */
	for (i = 0; i < SIZE; i++)  OCCUR[i] = 0;
	for (l = h; l >= 0; l--) if ((N = (vertex_p *)L[l])) {
		for (i = 0; N[i]; i++) {
			vertex_p	vertex = N[i];
			if (vertex->atom->degre < 2) OCCUR[vertex->atom->ID] = 1; 
			else {
				j = 0; if (vertex->parent) for (j = 0; vertex->parent[j]; j++);
				OCCUR[vertex->atom->ID] += j;
			}
		}
	}

	n = layer_size(L,h);
	s = (p_string)malloc(10*(n+1)*sizeof(element_t)); s[0] = '\0';

	/* reset arrays */
	N = L[h]; LL = 0;
	E = (edge_p)malloc((SIZE+1)*sizeof(edge_t));
	OCC = (t_integer *)malloc((SIZE+1)*sizeof(t_integer));
	for (i = 0; i < SIZE; i++) { E[i].y = E[i].l = NIL; OCC[i] = 0; }	
	for (n = 0; N[n]; n++) occur_string(N[n],E,OCC);
	for (i = 0; i < SIZE; i++) if (E[i].y) E[i].y[0] = -1;
	for (n = 0; N[n]; n++) {
		print_string(&s[strlen(s)],NIL,N[n],E,LAB,OCC);
		print_bond(&s[strlen(s)],N[n],N[n+1]); 
	}
	
	/* free arrays */
	for (i = 0; i < SIZE; i++) {
	    free(E[i].y); free(E[i].l); 
	}
	free(E); free(OCC);
	reset_print_string(L,h);
	return(s);
}

/*----------------------------------------------------------
 PRINT SIGNATURES IN SIG MODE
 -----------------------------------------------------------*/

LOCAL	t_integer	print_string_sig(s,p,r)
/***********************************************************
 Print tree r in string s
 p is the parent r is the current vertex
 ***********************************************************/
vertex_p	p,r;
p_string	s;
{
	t_integer	i = 0,j = 0,k,ii,oder;
	vertex_p	m;
	element_t	element;
	
	if (r == NIL) return(0);
	print_bond(s,p,r); 
	
	/* change element for both r and m */
	if (sig_parent(r,0,&m)) {
		strcpy(element,m->element);
		for (i = 0; i < strlen(element); i++) 
			if (element[i] == ',') break;
		if (i == strlen(element)) {
			strcpy(m->element,"[");
			strcat(m->element,m->atom->potential_type);
			strcpy(element,""); sprintf(element,"%d",++LL);
			strcat(m->element,",");
			strcat(m->element,element);
			strcat(m->element,"]");
		}
		strcpy(r->element,m->element); 
	}
	for (i = 0; i < strlen(r->element); i++) 
		if ((r->element[i] != '/') && (r->element[i] != '\\'))
			element[j++] = r->element[i];
	element[j++] = '\0';
	strcat(s,element);	
	i = strlen(&s[0]);
	if (r->child == NIL) return(i);
	
	/* recursion */
	for (k = 0; r->child[k]; k++); if (k < 1) return(i);
	ii = i; sprintf(&s[i++],"("); 
	
	/* sort edges according to stereochemistry */
	for (k = 0; r->child[k]; k++) {
		i += print_string_sig(&s[i],r,r->child[k]); 
	}
	sprintf(&s[i++],")");
	if (i == ii+2) i = ii;
	return(i);
}

LOCAL	t_err		sort_string_sig(r)
/***********************************************************
 sort the children of r
 ***********************************************************/
vertex_p	r;
{
	EXTERN	int 		cmp_vertex_invariant();
	t_integer	k;
	
	if (r == NIL) return(OK);
	
	/* sort children by invariant */
	for (k = 0; r->child[k]; k++);
	if (k == 0) return(OK);
	qsort(r->child,k,sizeof(vertex_p),cmp_vertex_invariant);
	
	/* recursion */
	for (k = 0; r->child[k]; k++) sort_string_sig(r->child[k]); 
	return(OK);
}

LOCAL	p_string	layer_print_string_sig(L,h)
/***********************************************************
 Print a signature
 ***********************************************************/
vertex_p *	L[];
t_integer	h;
{
	vertex_p	*N;
	p_string	s;
	t_integer	n = 0;
	
	N = L[h]; 
	if (N == NIL) return(NIL);
	if (N[0] == NIL) return(NIL);
	n = layer_size(L,h); LL = 0;
	s = (p_string)malloc(10*(n+1)*sizeof(element_t)); s[0] = '\0';
	for (n = 0; N[n]; n++) { 
		sort_string_sig(N[n]);
		print_string_sig(&s[strlen(s)],NIL,N[n]); 
		print_bond(&s[strlen(s)],N[n],N[n+1]); 
	}
	reset_print_string(L,h);
	return(s);
}

/*----------------------------------------------------------
 PRINT SIGNATURES IN CIP MODE (just the 2 first layer)
 -----------------------------------------------------------*/

LOCAL	t_integer	print_string_cip(s,r,h)
/***********************************************************
 Print tree r in string s
 up to h == 2
 ***********************************************************/
vertex_p	r;
p_string	s;
t_integer	h;
{
	t_integer	i = 0,k,ii,oder;
	t_bufstring	tagAtom;
	
	if (r == NIL) return(0);
	sprintf(tagAtom,"[%d,%d]",r->atom->ID,r->invariant); // To print invariants 
	sprintf(&s[0],"%s",tagAtom);
	i = strlen(&s[0]);
	if (r->child == NIL) return(i);
	if (h == 2) return(i); // !!!!!
	
	/* recursion */
	for (k = 0; r->child[k]; k++); if (k < 1) return(i);
	ii = i; sprintf(&s[i++],"("); 
	for (k = 0; r->child[k]; k++) {
		i += print_string_cip(&s[i],r->child[k],h+1); 
	}
	sprintf(&s[i++],")");
	if (i == ii+2) i = ii;
	return(i);
}

LOCAL	p_string	layer_print_string_cip(L,h)
/***********************************************************
 Print a signature
 ***********************************************************/
vertex_p *	L[];
t_integer	h;
{
	vertex_p	*N;
	p_string	s;
	t_integer	n = 0;
	
	N = L[h]; 
	if (N == NIL) return(NIL);
	if (N[0] == NIL) return(NIL);
	n = layer_size(L,h);
	s = (p_string)malloc(3*(n+1)*sizeof(element_t)); s[0] = '\0';
	for (n = 0; N[n]; n++) { 
		print_string_cip(&s[strlen(s)],N[n],0);
		print_bond(&s[strlen(s)],N[n],N[n+1]); 
	}
	return(s);
}

/*----------------------------------------------------------
  INVARIANTS
-----------------------------------------------------------*/

LOCAL	t_err	compute_vertex_invariant(vertex,INV,relation,l)
/***********************************************************
 Compute the invariant for vertex
 INV the atom invariant, relation is "" parent or child
 l is the current layer
 This function was changed 07/2012 by JLF
 ***********************************************************/
vertex_p	vertex;
t_integer	INV[],l;
t_bufstring	relation;
{
	t_integer *	invar;
	t_integer	n=0,i=0,j,inv=0;
	vertex_p *	neighbor = NIL;
	
	if (vertex == NIL) return(ERROR);
	// run only when called by init_invariant
	// it populates vertex->invlist[0...5]
	vertex_invariant(vertex,l);

	// compile and sort the invariant of the neigbors
	// the data structure invar comprises the 5 elements used
	// for CIP sorting 5Z, max bond, charge, cip, bons order ij
	// and the neighbor invariant itself (as the 6th element)
	if (strcmp(relation,"child") == OK) { neighbor = vertex->child; }
	else if (strcmp(relation,"parent") == OK) { neighbor = vertex->parent; }
	else { neighbor = NIL; }
	if (neighbor) for (n = 0; neighbor[n]; n++);
	if (n) invar = (t_integer *)calloc(6*n,sizeof(t_integer));  
	for (i = 0; i < n; i++) {
		// First, we add the invariant of the neighbor vertex
		for (j = 0; j <= 3; j++) invar[6*i+j] = neighbor[i]->invlist[j];
		// We add the term for the edge type
		invar[6*i+4] = edge_invariant(vertex,neighbor[i]);
		// Next, the invariant of its branch scaled 
		invar[6*i+5] = neighbor[i]->invariant;
	}
	// sort the 6 element-invar in increasing order
	if (n) qsort(invar,n,6*sizeof(t_integer),cmp_6_invariant);

	// append invar to vertex->invlist 
	// First we append the 5 elements used in cip sorting
	// next we add the invariants of the n neighbors and finally
	// the node invariant itself 
	 for (i = 0; i < n; i++)  
		for (j = 0; j < 5; j++)
			vertex->invlist[5*(i+1)+j] = invar[6*(n-1-i)+j]; // PC Changed i*6+j (n-1-i)*6+j (Heavier children first)
	 for (i = 0; i < n; i++)
			vertex->invlist[5*(n+1)+i] = invar[6*(n-1-i)+5]; // PC Changed i*6+j (n-1-i)*6+j (Heavier children first)
		
	// add the previous invariant (the least important)
	if (SCAN && INV) inv = INV[vertex->atom->ID]; else inv = vertex->invariant;
	vertex->invlist[5*(n+1)+n] = inv;
	vertex->invlist[6*(n+1)] = ENDINVLIST; // end of string
	if (n) free(invar);
	return(OK);
}

LOCAL	t_err	compute_layer_invariant(N,INV,relation,L,h,l)
/***********************************************************
 Compute invariant for all vertices in N[]
 L and h are passed for debuging prupose only
 l is the layer
 ***********************************************************/
vertex_p	N[];
t_integer	INV[];
t_bufstring	relation;
vertex_p *	L[];
t_integer	h,l;
{
	t_integer	i,n;
	t_integer	*invar,inv;
	
	for (n = 0; N[n]; n++)
	    compute_vertex_invariant(N[n],INV,relation,l);

	if (VERBOSE) { printf("Layer Invariants before rescale layer\n");print_layer(L,h,INV); } 
	qsort(N,n,sizeof(vertex_p),cmp_vertex_invlist); // JLF 02/2012
	invar = (t_integer *)malloc((n+1)*sizeof(t_integer));
	invar[0] = inv = 1; i = 1;
	while (i < n) {
		while (cmp_vertex_invlist(&N[i],&N[i-1]) == 0) { // JLF 02/2012
	        invar[i] = inv; i++;
		}
		inv++; invar[i] = inv; i++;
	}
	for (i = 0; i < n; i++) N[i]->invariant = invar[i];
	free(invar);
	if (VERBOSE) { printf("Layer Invariants after rescale layer \n");print_layer(L,h,INV); }
	return(OK);
}

LOCAL	t_err	compute_invariant_occurence(L,h,INV,OCC)
/***********************************************************
 Compute invariant and occurence for all atoms
 and return the number of invariants
 At the end of the process all invariants are rescaled
 1..N and are integers
***********************************************************/
vertex_p *	L[];
t_integer	h, INV[], OCC[];
{
	vertex_p *	N;
	t_integer	l,nb,i0,inv,i;
	inocc_p inocc;

	inocc = (inocc_p) calloc(SIZE,sizeof(inocc_t));

	for (i=0;i<SIZE;i++) inocc[i].ID = i;
	for (l = 0; l < h+1 ; l++) 
		if ((N = (vertex_p *)L[l]))
	     for (i = 0; N[i]; i++)
			 LAYER[N[i]->atom->ID][l] = N[i]->invariant;
	qsort(inocc,SIZE,sizeof(inocc_t),cmp_compute_invariant_occurence);

	// if (VERBOSE) { printf("before rescale init invariant\n"); print_layer(L,h,INV); }
	inocc[0].inv = inv = 1;
	for (i = 1; i < SIZE; i++) 
	  if (cmp_compute_invariant_occurence(&inocc[i],&inocc[i-1]) == 0)
		  inocc[i].inv = inocc[i-1].inv;
	  else inocc[i].inv = (++inv);
	for (i = 0; i < SIZE; i++) INV[inocc[i].ID] = inocc[i].inv;
	// if (VERBOSE) { printf("after rescale init invariant\n"); print_layer(L,h,INV); }

	/* populate INOCC[].occ
	   occ the nbr of atoms having the same given invariant
	   inv the invariant
	   ID the atom number
	*/
	nb = 1; i0 = 0;
	for (i = 1; i < SIZE; i++) {
	  t_integer	k;
	  if (inocc[i].inv == inocc[i-1].inv) { nb++; continue;}
	  for (k = i0; k < i0+nb; k++) inocc[k].occ = nb;
	  nb = 1; i0 = i; 
	  }
	for (i = i0; i < SIZE; i++) inocc[i].occ = nb;
	for (i = 0; i < SIZE; i++) 
	  if (OCC[inocc[i].ID] > 1) OCC[inocc[i].ID] = inocc[i].occ;
	free(inocc);
	return(inv);
	}

/*----------------------------------------------------------
  INVARIANTS
-----------------------------------------------------------*/

LOCAL	t_err	end_invariant(L,h)
/***********************************************************
 print the signature string 
***********************************************************/
vertex_p *	L[];
t_integer	h;
{
	p_string	s;
	t_integer	i;

	for (i = 0; i < SIZE; i++) LABEL[i] = LACUR[i] = -1; NBCUR = 0; 
	s = layer_print_string(L,h,LABEL); NSTR++;
	if (SMAX) if (strcmp(s,SMAX) < 0) {free(s); return(OK);}
	if (SMAX) free(SMAX); SMAX = s; 
	for (i = 0; i < SIZE; i++) LACAN[i] = LACUR[i];
	return(OK);
	}


LOCAL	t_err	compute_invariant(L,h,OCC,INV,ITER)
/***********************************************************
 Compute label and invariant for all vertices.
***********************************************************/
vertex_p *	L[];
t_integer	h;
t_integer 	OCC[];
t_integer	INV[];
t_integer	ITER;
{
	t_integer	l,i;
	t_integer	omax = 1,imax = -1, inv = -1;
	vertex_p *	N; 
	t_integer	*occur = NIL;
	t_integer	*invar = NIL; // PC: Added initialization

	if (L == NIL) return(ERROR);
	if (SCAN == FALSE) {
		/* SIG and CIP compute invariants from children to parent */
		t_integer	l0 = 0, ln = h+1, li = 1;
		if (CIP) { // CIP MODE
			/*********************************/
			/* Special algorithm for CIP mode */
			/*********************************/
			compute_cip(L,h);
			//for (l = l0; l != ln; l = l+li)
			//	if ((N = (vertex_p *)L[l]))
			//		compute_layer_invariant(N,invar,"child",L,h,l);
			/*********************************/
			SMAX = layer_print_string_cip(L,h);
		}
		else { // SIG MODE
			for (l = l0; l != ln; l = l+li)
				if ((N = (vertex_p *)L[l]))
					compute_layer_invariant(N,invar,"child",L,h,l);
			if (VERBOSE) { printf("CIP child\n");print_layer(L,h,invar); }
			SMAX = layer_print_string_sig(L,h);
		}
		return(OK);
	}
	
	/* SCAN mode initialization */
	occur = (t_integer *)malloc((SIZE+1)*sizeof(t_integer));
	invar = (t_integer *)malloc((SIZE+1)*sizeof(t_integer));
	for (i = 0; i < SIZE; i++) { invar[i] = INV[i]; occur[i] = OCC[i]; }
	
	LOOP {
		/* compute invariants from parent to children */
		t_integer l0 = h, ln = -1, li = -1;
		for (l = l0; l != ln; l = l+li) 
			if ((N = (vertex_p *)L[l])) {
				// vertex_invariant(*N,l);
				compute_layer_invariant(N,invar,"parent",L,h,l);
			}
		if (VERBOSE) { printf("before rescale parents\n");print_layer(L,h,invar); }
		i = compute_invariant_occurence(L,h,invar,occur);
		if (VERBOSE) { printf("after rescale parents\n");print_layer(L,h,invar); }
		
		/* compute invariants from children to parent */
		l0 = 0; ln = h+1; li = 1;
		for (l = l0; l != ln; l = l+li) 
			if ((N = (vertex_p *)L[l])) {
				// vertex_invariant(*N,l);
				compute_layer_invariant(N,invar,"child",L,h,l);
			}
		if (VERBOSE) { printf("before rescale children\n");print_layer(L,h,invar); }
		i = compute_invariant_occurence(L,h,invar,occur);
		if (VERBOSE) { printf("after rescale children\n"); print_layer(L,h,invar); }
	    if (i <= inv) break; inv = i;
	} inv = -1;
	if (VERBOSE) { printf("after loop parent children children\n"); print_layer(L,h,invar); }
	/* Process ends here in CIP mode */
	if (CIP) { // this is the CIP in SCAN mode
		SMAX = layer_print_string_cip(L,h);
		free(invar); free(occur);
		return(OK);
	}
	/* At this stage computation of invariant is completed
	   and all invariants are integers.
	   If several atoms have the same invariants
	   the way to break up tie is to add a weigt (< 1) 
	   to each atoms (one after another) with the same 
	   invariant and call again compute invarinats
	 */

	/* last step */
	for (i = 0; i < SIZE; i++) if (occur[i] > 1) break;
	if (i == SIZE) {
	   end_invariant(L,h); free(invar); free(occur);
	   return(OK);
	   }
	/* find the orbit to singularize from leaves to root 
	   this orbit has the maximum number of elements
           and the max invariant */
	for (l = 0; l <= h; l++) if ((N = (vertex_p *)L[l])) 
	  for (i = 0; N[i]; i++) 
             if (COLOR[N[i]->atom->ID] >  1)
             if ((occur[N[i]->atom->ID] > omax) && (invar[N[i]->atom->ID] > inv)) {
	            imax = N[i]->atom->ID; inv = invar[imax]; omax = occur[imax]; 
	        }

	/* recursion */
	if (ITER > NITR) NITR = ITER;
	if (omax > 1) for (i = 0; i < SIZE; i++) 
	   if (invar[i] == invar[imax]) {
		  invar[i] += SIZE;
//		  printf("ITER = %d root = %d labeled atom = %d imax = %d omax = %d invar %d MAX ITER = %d\n",\
//								 ITER,(*L[h])->atom->ID,i,imax,omax,invar[i],NITR);  // PC Commented lines
  	      compute_invariant(L,h,occur,invar,ITER+1);
		  invar[i] -= SIZE;
	      /* when the number of iteration is greater than the
			maximum number of color allowed only one atom is
			colored. sisi_inq_color is set to a value lower
			than MAXINTEGER when the option scanx is used
			with x < MAXINTEGER */
	      if (ITER >= sisi_inq_color()) break;

	      }
	/* else omax > 1 */
	if (omax == 1) compute_invariant(L,h,occur,invar,ITER+1);
	free(invar); free(occur);
	return(OK);
	}
	
LOCAL	t_err	init_invariant(L,h)
/***********************************************************
 Compute initial invariant for all vertices which is
 the number of parents in the DAG data str.
 OCCUR is the occurence number in the signature-tree
 COLOR is the occurence number in the signature-tree
 for bactracking/coloring purposes
 INVAR is the array of initial invariants
 LABEL is for printing only not invariant computation
 it is used for ring closer
 ***********************************************************/
vertex_p *	L[];
t_integer	h;
{
	t_integer	l,i,j;
	vertex_p *	N;
	
	if (L == NIL) return(ERROR);	
	if (SCAN) {
		if (LABEL == NIL) {
			LABEL = (t_integer *)malloc((SIZE+1)*sizeof(t_integer));
			OCCUR = (t_integer *)malloc((SIZE+1)*sizeof(t_integer));
			COLOR = (t_integer *)malloc((SIZE+1)*sizeof(t_integer));
			INVAR = (t_integer *)malloc((SIZE+1)*sizeof(t_integer));
			LAYER = (t_integer **)malloc((SIZE+1)*sizeof(t_integer *));
			LACAN = (t_integer *)malloc((SIZE+1)*sizeof(t_integer));
			LACUR = (t_integer *)malloc((SIZE+1)*sizeof(t_integer));
		}
		for (i = 0; i < SIZE; i++) {
			OCCUR[i] = COLOR[i] = INVAR[i] = 0; LABEL[i] = LACAN[i] = LACUR[i] = -1;
			LAYER[i] = (t_integer *)malloc((h+2)*sizeof(t_integer));
			for (j = 0; j < h+1; j++) LAYER[i][j] = 0;
			LAYER[i][h+1] = ENDINVLIST; // end of string
		}
		
		/* vertices with degree 1 have OCCUR = 1 
		 JLF 10/03 vertices occuring alone 
		 with more than one parent have OCC += 1
		 all other vertices have OCC += nbr of parent 
		 each time they occur 
		 INVAR are initialize corresponding to the 
		 */
		
		for (l = h; l >= 0; l--) if ((N = (vertex_p *)L[l])) {
			t_integer	parent = 0;
			
			/* count the number of vertices with more
             than one parent in layer l */
			for (i = 0; N[i]; i++) {
				vertex_p	vertex = N[i];
				j = 0; if (vertex->parent) for (j = 0; vertex->parent[j]; j++);
				if (j > 1) parent += 1;                      
			}
			for (i = 0; N[i]; i++) {
				vertex_p	vertex = N[i];
				if (vertex->atom->degre < 2) {
					OCCUR[vertex->atom->ID] = COLOR[vertex->atom->ID] = 1; 
				}
				else {
					j = 0; if (vertex->parent) for (j = 0; vertex->parent[j]; j++);
					OCCUR[vertex->atom->ID] += j;
					if (parent < 2)  COLOR[vertex->atom->ID] += 1; 
					else             COLOR[vertex->atom->ID] += j; 
				}
			}
		}
	}
	// all modes 
	for (l = 0; l != h+1; l++) 
		if ((N = (vertex_p *)L[l])) 
			compute_layer_invariant(N,NIL,"",NIL,h,l);			
	if (SCAN == FALSE) return(OK);
	if (VERBOSE) { printf("init before rescale\n");print_layer(L,h,INVAR); }
	compute_invariant_occurence(L,h,INVAR,OCCUR);
	if (VERBOSE) { printf("init after rescale\n");print_layer(L,h,INVAR); }
	return(OK);
}

DEFINE	p_string sicd_signature_atom(a1,a2,h)
/***********************************************************
 Build a signature of height h for atoms a1 and a2 
***********************************************************/
atom_p		a1,a2;
t_integer	h;
{
	vertex_p *	*L;
	t_integer	i,t1,t2;
	p_string	scan_type;
	
	SIZE = ((molecule_p)((group_p)a1->group)->molecule)->size;
	scan_type = (p_string)xsio_scan_type();
	if (strstr(scan_type,"scan")) SCAN = TRUE; else SCAN = FALSE;
	if (strstr(scan_type,"sig"))  SIG = TRUE;  else SIG = FALSE;
	if (strstr(scan_type,"cip"))  CIP = TRUE;  else CIP = FALSE;
	if (h > SIZE+1) h = SIZE+1;

	if (CIP)  {// PC skip in CIP mode unnecessary atoms
 		if (!strcmp(a1->comment,"0")) {return(NIL);}
 		else if (a2 != NIL) {
 			if (!strcmp(a2->comment,"0")) {return(NIL);}
 		}
 	}


	/* build dag */
	L = (vertex_p * *)malloc((h+2)*sizeof(vertex_p *));
	for (i = 0; i < h+2; i++) L[i] = NIL;
	build_dag(a1,a2,L,h);
	if (VERBOSE) printf("building dag completed for atom %s %d\n",a1->name,a1->ID+1);

	/* compute label invariant */
	init_invariant(L,h); 
	if (VERBOSE) printf("init invariant completed\n");
	SMAX = NIL; NSTR = 0; NITR = 0;
	if (VERBOSE) {
	   t1 = clock();
	   }
	compute_invariant(L,h,OCCUR,INVAR,0);
	if (VERBOSE) {
	   t2 = clock();
	   printf("vertex %d iteration = %d permutations %d CPU TIME %.2f s\n",a1->ID+1,NITR,NSTR,(float)(t2-t1)/(float)1000000L);
	   }

	/* delete layer data structures */
	delete_dag(L,h);
	if (VERBOSE) printf("atom = %d: %s\n",a1->ID+1,SMAX);
	return(SMAX);
	}

DEFINE	p_string sicd_signature_atom_label(a,h,label)
/***********************************************************
 Build a signature of height h for atom a
***********************************************************/
atom_p		a;
t_integer	h;
t_integer	label[];
{
	vertex_p *	*L;
	t_integer	i,t1,t2;
	atom_p		a1 = a, a2 = NIL;
	
	SIZE = ((molecule_p)(((group_p)a->group)->molecule))->size;
	SCAN = TRUE; SIG = FALSE; CIP = FALSE;
	if (h > SIZE+1) h = SIZE+1;
	
	/* build dag */
	L = (vertex_p * *)malloc((h+2)*sizeof(vertex_p *));
	for (i = 0; i < h+2; i++) L[i] = NIL;
	build_dag(a1,a2,L,h);
	if (VERBOSE) printf("building dag completed for atom %s %d\n",a1->name,a1->ID+1);
	
	/* compute label invariant */
	init_invariant(L,h); 
	if (VERBOSE) printf("init invariant completed\n");
	SMAX = NIL; NSTR = 0; NITR = 0;
	if (VERBOSE) {
		t1 = clock();
	}
	compute_invariant(L,h,OCCUR,INVAR,0);
	if (VERBOSE) {
		t2 = clock();
		printf("vertex %d iteration = %d permutations %d CPU TIME %.2f s\n",a1->ID+1,NITR,NSTR,(float)(t2-t1)/(float)1000000L);
	}
	
	/* delete layer data structures */
	delete_dag(L,h);
	if (VERBOSE) printf("atom = %d: %s\n",a1->ID+1,SMAX);
	for (i = 0; i < SIZE; i++) label[i] = LACAN[i];
	return(SMAX);
	}

