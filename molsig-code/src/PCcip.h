/*
 * PCcip.h
 *
 *  Created on: 9 August 2012
 *      Author: carbonell
 */

#ifndef PCCIP_H_
#define PCCIP_H_

#include "general.h"
#include "signature.h"

// PC: These definitions are in candag.c, they should be moved to a common header file (or make extern)

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


// This structure contains the information needed in order to order vertices
// based on cip priorities, other terms can be added if necessary
typedef struct cipnode {
	vertex_p vertex; // The list of vertices of the branch at the current layer
	struct cipnode* children;
	int cipinvariant;
	int nchildren; // Number of children of the vertex
	int nvertex; // Number of vertices in the list
} cipnode_t, *cipnode_p;

// This structure is used in order to associate a branch with the list of vertices
// at the current layer and its corresponding invlist for the CIP comparison
typedef struct BRANCH {
	vertex_p parent;  // The initial parent at top of the branch
	cipnode_p cip; // The list of subbranches of the branch at the current layer
	int ncip; //Number of cip nodes
	int ncipinv; //Number of cip invariants in the layer (Rassat et al, criteria)
} branch_t,*branch_p;



t_err compute_cip(vertex_p * L[],t_integer h);

#endif /* PCCIP_H_ */
