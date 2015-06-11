/*
 * PCcip.c
 *
 *  Created on: 9 August 2012
 *      Author: carbonell
 *      Special CIP algorithm
 */

#include "PCcip.h"
// A negative value means that nodes are sorted from heavier to lighter
#define CMPS -1

// Print CIP comparisons
t_bool VERBOSECIP = FALSE;

// If TRUE, apply cip invariants as in Rassat et al.
// Otherwise, apply IUPAC rules
t_bool CIPTEST = FALSE;

t_integer	sig_parent(n,l)
/***********************************************************
 Search if n has a parent having the same atom
 l is the layer of n
 The routine return * the layer of the parent
 ***********************************************************/
vertex_p	n;
t_integer	l;
{
	vertex_p	m = NIL;
	t_integer	lm = l+2,i;
	element_t	element;

	if (n == NIL) return(0);
	if (n->parent == NIL) return(0);
	for (m = n->parent[0]; m; m = m->parent[0]) {
		if (m->atom == n->atom) break;
		if (m->parent == NIL) break;
		lm--;
	}
	if (m == NIL) return(0);
	if (m->atom != n->atom)
	return(0);
	return(lm);
}


/******************************************
 * Init branch structure
 *****************************************/
branch_t init_branch(branch_p branch,vertex_p N) {
	branch_t init;
	init.parent = N; // Store the initial parents of the branches
	init.cip = (cipnode_p )calloc(2,sizeof(cipnode_t));
	init.cip[0].vertex = N;
	init.cip[0].nvertex = 1;
	init.cip[0].cipinvariant = 0;
	init.ncip = 1;
	init.ncipinv = 0;  // PC TEST
	return(init);
}

/******************************************
 * Deallocate memory of branch structure
 *****************************************/
void 	free_branch(branch_p* branch,vertex_p* N,int h) {
	int k,l,m;
	for(k=0;N[k];k++) {
		for(l=0;l<h;l++) {
			if (branch[k][l].parent!=NIL) {
				for(m=0;branch[k][l].cip[m].children!=NIL;m++)
					free(branch[k][l].cip[m].children);
				free(branch[k][l].cip);
			}
		}
		free(branch[k]);
	}
	free(branch);
}

/******************************************
 * If any of the vertex appeared before in a higher layer
 * the one that appeared at the highest is greater
 *****************************************/
int cip_cmp(cipnode_p n1, cipnode_p n2) {
	int sign = 1;
	if (n1->vertex==NIL) {
		return 0;
	}
	if (n2->vertex==NIL) {
		return 0;
	}
	if (CIPTEST) return 0;
	if (n1->cipinvariant==0)
		if (n2->cipinvariant==0)
			return 0;
		else return -1*sign;
	else if (n2->cipinvariant==0)
		return 1*sign;

	if (n1->cipinvariant < n2->cipinvariant) // Appears in a higher layer
			return 1*sign;
		else if (n1->cipinvariant > n2->cipinvariant)
			return -1*sign;
	return 0;
}

// Criteria in Rassat et al.
int cip_branch_cmp(branch_p b1,branch_p b2) {
	int sign = -1;
	if (b1->ncipinv > b2->ncipinv) return 1*sign;
	else if  (b1->ncipinv < b2->ncipinv) return -1*sign;
	return 0;
}


/******************************************
 * Compare invariants of a cip node
 * This is the most basic comparison
 *****************************************/
int invvertex_cmp(cipnode_p n1, cipnode_p n2) {
	int sign = CMPS; //Invert order from heavier to lighter
	if ((n1== NIL)||(n1->vertex==NIL))
		if ((n2== NIL)||(n2->vertex==NIL))
			return 0;
		else return -1*sign;
	else if ((n2==NIL)||(n2->vertex==NIL))
		return 1*sign;
	if ((n1)->vertex->invariant>(n2)->vertex->invariant) {
		return 1*sign;
	} else {
		if ((n1)->vertex->invariant<(n2)->vertex->invariant)
			return -1*sign;
	}
	return CMPS*cip_cmp(n1,n2);
}

/******************************************
 * Compares two nodes based on their children
 * in order to sort them from heavier to lighter
 * Useful when we are moving to the next layer
 *****************************************/
int children_cmp(cipnode_p n1,cipnode_p n2)
{
	int n;
	int  sign= CMPS;
	if (n1 == NIL)
		if (n2 ==NIL)
			return 0;
		else return -1*sign;
	else if (n2 ==NIL)
		return 1*sign;
	for(n=0;(n1)->children[n].vertex!=NIL;n++) {
		if ((n2)->children[n].vertex!=NIL) {
			if (((n1)->children[n].vertex)->invariant > ((n2)->children[n].vertex)->invariant)
				return 1*sign;
			else if ((n1)->children[n].vertex->invariant<(n2)->children[n].vertex->invariant)
				return -1*sign;
		} else return 1*sign;
		int cip;
		if ( (cip=cip_cmp( &((n1)->children[n]) ,&((n2)->children[n]) ) )!=0) return sign*cip;
	}
	if ((n2)->children[n].vertex!=NIL) return -1*sign;
	return 0;
}



//
/******************************************
 * Adds a new layer to the branches (if it no already done)
 * Returns the number of children in the current layer,
 * i.e., the vertices of the new layer
 *****************************************/
int add_branch_layer(branch_p thisbranch,int l) {
	branch_t nextbranch;
	cipnode_p vsort;
	int n,m,nv=0;

	// Check if we are at the end of the branch
	for(n=0;n<thisbranch[l].ncip;n++)
		if (thisbranch[l].cip[n].vertex!=NIL)
			for(m=0;thisbranch[l].cip[n].vertex->child[m]!=NIL;m++) nv++;
	if (nv==0) return(nv); // No more children

	// Check if layer was already stored
	if (thisbranch[l+1].parent!= NIL)
		return(thisbranch[l+1].ncip);

	// Count how many children,then store and sort them
	for(n=0;n<thisbranch[l].ncip;n++) {
		if (thisbranch[l].cip[n].vertex==NIL) continue;
		for(m=0;thisbranch[l].cip[n].vertex->child[m]!=NIL;m++) nv++;
		thisbranch[l].cip[n].children = calloc(m+1,sizeof(cipnode_t));
		for(m=0;thisbranch[l].cip[n].vertex->child[m]!=NIL;m++) {
			thisbranch[l].cip[n].children[m].vertex =  thisbranch[l].cip[n].vertex->child[m];
			thisbranch[l].cip[n].children[m].cipinvariant = sig_parent(thisbranch[l].cip[n].vertex->child[m],l+1);
		}
		qsort(thisbranch[l].cip[n].children,m,sizeof(cipnode_t),invvertex_cmp);
	}

	// Copy in vsort the vertices
	vsort = (cipnode_p)calloc(n+1,sizeof(cipnode_t));
	for(n=0;n<thisbranch[l].ncip;n++)
		vsort[n] = thisbranch[l].cip[n];
	// Initially, we sort the parents according to their children
	if (l==0)
		qsort(vsort,n,sizeof(cipnode_t),children_cmp);
	else {
	// Sort now the vertices of the layer according to their children
	// in case that they have the same priorities
	if (n>0) {
		for(n=1;n<thisbranch[l].ncip;n++) {
			int dup = 1;
			if (vsort[n-1].vertex==NIL) continue;
			if (vsort[n].vertex==NIL) continue;
			while(vsort[n-1].vertex->invariant == vsort[n].vertex->invariant) {
				dup++;
				n++;
				if ((n==thisbranch[l].ncip)||(vsort[n].vertex==NIL)) break;
			}
			if (dup>1) qsort(&(vsort[n-dup]),dup,sizeof(cipnode_t),children_cmp);  // PC CHECK
		}
	}
	}


	// Define the next branch
	nextbranch.parent = thisbranch[l].parent;
	// Allocate a pointer to each children vertex
	nextbranch.cip = (cipnode_p)calloc(nv+n+1,sizeof(cipnode_t));
	nv = 0;
	int nvertex = 0;
	int ncipinv = 0; // PC TEST
	for(n=0;n<thisbranch[l].ncip;n++) {
		if (vsort[n].children==NIL) continue;
		for(nvertex=0;vsort[n].children[nvertex].vertex!=NIL;nvertex++);
		for(m=0;vsort[n].children[m].vertex!=NIL;m++) {
			nextbranch.cip[nv].vertex = vsort[n].children[m].vertex;
			nextbranch.cip[nv].cipinvariant = vsort[n].children[m].cipinvariant; // PC ASSIGN FROM NEWCHILDREN
			nextbranch.cip[nv].nvertex = nvertex;
			if (nextbranch.cip[nv].cipinvariant!=0) { // PC TEST
				ncipinv++;
			}
			nv++;
		}
		// If there are no more children, store a null pointer
		if (m==0) {
			nextbranch.cip[nv].nvertex = 0;
			nextbranch.cip[nv++].vertex = NIL;
		}
	}
	nextbranch.ncip = nv;
	nextbranch.ncipinv = ncipinv;
	free(vsort);
	thisbranch[l+1] = nextbranch;
	return(thisbranch[l+1].ncip);
}

/******************************************
 * Print a comparison between the two branches
 *****************************************/
void print_branch_cmp(branch_t n1,branch_t n2) {
	int i;
	if (n1.parent!=NIL) {
		if (n1.parent->parent!=NIL)
			printf("{%d} ",n1.parent->parent[0]->atom->ID+1);
		 else printf("{} ");
		printf("%d:",n1.parent->atom->ID+1);
		for(i=0;i<n1.ncip;i++)
			if (n1.cip[i].vertex!=NIL)
				printf(" %d[%d,%d,%d]",n1.cip[i].vertex->atom->ID+1,n1.cip[i].vertex->invariant,n1.cip[i].nvertex,n1.cip[i].cipinvariant);
			else
				printf(" [0]");
	}
	printf(" [[%d]] ",n1.ncipinv);
	printf(" <> ");
	if (n2.parent!=NIL) {
		if (n2.parent->parent!=NIL)
			printf("{%d} ",n2.parent->parent[0]->atom->ID+1);
		else printf("{} ");
		printf("%d:",n2.parent->atom->ID+1);
		for(i=0;i<n2.ncip;i++) {
			if (n2.cip[i].vertex!=NIL)
				printf(" %d[%d,%d,%d]",n2.cip[i].vertex->atom->ID+1,n2.cip[i].vertex->invariant,n2.cip[i].nvertex,n2.cip[i].cipinvariant);
			else
				printf(" [0]");
		}
	}
	printf(" [[%d]] ",n2.ncipinv);
	printf("\n");
}

/******************************************
 * Compares two branches (at the current layer)
 *****************************************/
int branch_cmp(branch_p n1,branch_p n2)
{
	int n;
	int  sign= 1;
	if (n1 == NIL)
		if (n2 ==NIL)
			return 0;
		else return -1*sign;
	else if (n2 ==NIL)
		return 1*sign;




	// While both branches are not empty
	// compare the vertex invariants
	for(n=0;n<n1->ncip;n++) {


		if (n>=n2->ncip) {
			if (n1->cip[n].vertex!= NIL)
				return 1*sign;
			else
				continue;
		}



		if (n1->cip[n].vertex== NIL) {
			if (n2->cip[n].vertex!= NIL)
				return -1*sign;
		} else if (n2->cip[n].vertex== NIL)
				return 1*sign;
		else {
			if (n1->cip[n].vertex->invariant > n2->cip[n].vertex->invariant)
				return 1*sign;
			else if (n1->cip[n].vertex->invariant < n2->cip[n].vertex->invariant)
				return -1*sign;
		}
		// If there is a tie, check for the duplicate vertex rule
		int cip;
		if ((cip = cip_cmp(&((n1)->cip[n]),&((n2)->cip[n])))!=0) return cip;

		// Finally, make sure that each vertex has the same number of children
		//if (n1->cip[n].nvertex > n2->cip[n].nvertex) return 1*sign;
		//else if  (n1->cip[n].nvertex < n2->cip[n].nvertex) return -1*sign;
	}
	// Continue with n2 to see if there are no empty branches
	for(;n<n2->ncip;n++) {
		if (n2->cip[n].vertex!= NIL)
			return -1*sign;
	}
	// This is the Rassat et al criteria, only for comparison purposes
	if (CIPTEST) return cip_branch_cmp(n1,n2);
	return 0;
}

/******************************************
 * Compares two branches (at the current layer)
 * starting from the top layer until they differ
 *****************************************/
int full_branch_cmp(branch_p *n1, branch_p *n2) {
	int l = 0;
	int r;
	int sign = CMPS;
	int b1,b2;
	// while branches are the same, compare at next layer
	if (VERBOSECIP==TRUE){
		printf("Layer %d ",l);
		print_branch_cmp(*n1[0],*n2[0]);
	}
	// While two branches are equal, descend to the next layer
	while((r = branch_cmp( &((*n1)[l]),&((*n2)[l]) ) )==0) {
		// Add branch layer if not already computed
		b1 = add_branch_layer(*n1,l);
		b2 = add_branch_layer(*n2,l);
		if (VERBOSECIP==TRUE){
			if ( (b1>0)&&(b2>0) ) {
				printf("Layer %d ",l+1);
				print_branch_cmp((*n1)[l+1],(*n2)[l+1]);
			}
		}
		if (b1==0)
			if (b2==0)
				return 0;
			else
				return -1*sign;
		else if (b2==0)
			return 1*sign;
		l += 1;
	}
	return(r*sign);
}

/******************************************
 * Apply CIP rules from the given layer
 *****************************************/
t_err cip_layer(vertex_p *N,int h) {

	t_integer	k,l,m;
	branch_p *branch;

	// For each children
	for (k=0;N[k];k++);
	branch = (branch_p *)calloc(k,sizeof(branch_p));
	for (k=0;N[k];k++)
			branch[k] = (branch_p)calloc(h,sizeof(branch_t));

	// Initially, each branch contains one of the children of the root
	// They are stored at 0
	for (k=0;N[k];k++)
		branch[k][0] = init_branch(branch[k],N[k]);

	// Sort branches
	qsort(branch,k,sizeof(branch_p*),full_branch_cmp);
	//Assign invariant
	int *invlist = calloc(k,sizeof(int));
	int inv = k;
	invlist[0] = inv;
	for (k=1;N[k];k++) {
		if (full_branch_cmp(&branch[k-1],&branch[k])!=0) inv--;
		invlist[k] = inv;
	}
	// If not all invariants are different, we shift them
	for (k=0;N[k];k++) {
			branch[k][0].parent->invariant = invlist[k]-(inv-1);
	}
	free(invlist);
	if (VERBOSECIP==TRUE) {
		printf("INV:");
		for (k=0;N[k];k++) {
			printf(" %d[%d]",branch[k][0].parent->atom->ID+1,branch[k][0].parent->invariant);
		}
		printf("\n");
	}
	free_branch(branch,N,h);
	return(OK);

}

/***********************************************************
 This routine applies the CIP rule in order to determine priorities
 of children. The algorithm goes from parent to children and stops
 when children have different configurations.
 ***********************************************************/
t_err compute_cip(L,h)
vertex_p *	L[];
t_integer	h;
{
	vertex_p *	N;

	// Compute cip rules for the children branches of the atom or the bond
	if ((N = (vertex_p *)L[h-1]))  cip_layer(N,h);

	return(OK);
}

