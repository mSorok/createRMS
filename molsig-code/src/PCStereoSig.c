/*
 * PCStereoSig.c
 *
 *  Created on: March 2012
 *      Author: Pablo Carbonell
 *      Routines that are used in order to iteratively assign the parity from the stereochemistry.
 *      It makes use of the information from:
 *      	a) The geometry of the molecular structure
 *      	b) Branch priorities determined through the application of the CIP rules
 */

#include "general.h"
#include "signature.h"
#include "DFgeneral.h"
#include "PCStereoSig.h"

t_bool PCINFO= FALSE;
// Structure containing the info about the signature and the atoms
struct SigGraph{
	int size;
	int nbonds;
	struct SigInfo* sig;
	struct SigInfo* bsig;
	struct AtomInfo* atom;
	struct BondInfo* bond;
};

// Stores the info about each atom and its associated invariants
struct AtomInfo {
	t_real* rootpriorities;
	int* children;
	int branches;
};

// Stores the info about each bond and its associated invariants
struct BondInfo {
	t_real* rootpriorities;
	int head;
	int* children_head;
	int branches_head;
	int tail;
	int* children_tail;
	int branches_tail;
};

// Stores the info about the signature
struct SigInfo {
	char *s;
	int nlabels;
	struct LabelPosition* labels;
};

// Store the position of atom labels in the signature and their corresponding weights
struct LabelPosition {
	char* start;
	char* end;
	int color;
	int id;
	t_real weight;
	int depth;
	int dbond;
};

// Structure creation and freed
struct SigGraph* pcss_NewSigGraph(int size, int nbonds);
void pcss_AddAtomSigGraph(struct SigGraph *sg,signature_p s,t_bool SIGSTRING);
void pcss_AddBondSigGraph(struct SigGraph *sg, signature_p s,t_bool SIGSTRING);
void pcss_FreeSigGraph(struct SigGraph *sg);
void pcss_OrderChildrenAsInGraph(struct SigGraph *sg,struct GRAPH Graph);
void pcss_OrderChildrenBondAsInGraph(struct SigGraph *sg,struct GRAPH Graph);
t_err 	pcss_AtomList(molecule_p 	m, struct GRAPH Graph);
void pcss_InitialParities(int* Parity,int* BondParity,struct GRAPH Graph);

void pcss_FlagStereoCenters(int size, int* newparities, struct GRAPH Graph) ;
void pcss_FlagStereoBonds(int nbonds, int* newbondparities, struct GRAPH Graph);

// Diameter for computing the atom or the bond signature
int pcss_DiameterAtomSignature(int init_dim);
int pcss_DiameterBondSignature(int init_dim);

// Reading and writing signatures
signature_p	pcss_ReadSignature(set_molecule_p SM);

// Check about parity update
t_bool pcss_ParitiesChange(int size, int* oldparities,int* newparities);
t_bool pcss_AtomParitiesChange(int size,int* currentparities,struct GRAPH Graph);
t_bool pcss_BondParitiesChange(int nbonds, int* currentbondpriorities,struct GRAPH Graph);

// Global stereo assignment based on given priorities
void pcss_GlobalStereoBonds(struct SigGraph *sg,struct GRAPH Graph,int *BondParity);
void pcss_GlobalStereoCenters(struct SigGraph *sg,struct GRAPH Graph,int *Parity);

void pcss_GetAtomInfo(int size,struct LabelPosition *labels,int nlabels, struct AtomInfo *atom);
void pcss_GetBondInfo(int size,struct LabelPosition *labels,int nlabels, struct BondInfo *bond);
void pcss_PrintStereobondInfo(struct SigGraph *sg,	struct GRAPH Graph,int i,int state);
void pcss_PrintStereocenterInfo(struct SigGraph *sg, int i,int state,struct GRAPH Graph);
void pcss_UpdateAtomParity(int* Parity,int* OldParity,struct GRAPH Graph);
void pcss_UpdateBondParity(int* BondParity,int* OldBondParity,struct GRAPH Graph);

int HeadBond(char *sigcursor,t_bool SIGSTRING);
int TailBond(char *sigcursor,t_bool SIGSTRING);
int HowManyLabels(char *sigcursor);
int RootAtom(char *sigcursor,t_bool SIGSTRING);
int BranchPriorities(t_real* Priorities,int size,int nlabels,struct LabelPosition* AtomLabels,int depth,t_bool RESET);

struct LabelPosition* pcss_ParseSignature(int nlabels,char* s,t_bool SIGSTRING);

// External calls not following header prototype convention
EXTERN t_err	dast_all_atom(molecule_p molecule,group_p group,f_traitement _function,t_bool marq,t_pointer arg);
EXTERN signature_p sisc_signature(t_pointer HEAP,molecule_p molecule,group_p group,atom_p atom,bond_p bond);
EXTERN t_err  sisi_inq_dimension_min();
EXTERN t_err  sisi_inq_dimension_max();
EXTERN t_err  sisi_set_dimension(t_integer d);


//  Main routine for computing and assigning the parities of stereocenters
void pcss_StereocenterParities(int* Parity,int* OldParity,struct GRAPH Graph, set_molecule_p SM) {
	signature_p s = NIL;
	int i;
	struct SigGraph* pcsg;

	pcss_AtomList(SM->molecule,Graph);
	pcsg = pcss_NewSigGraph(Graph.size, Graph.NumbEdges);
	// CALL THE CIP ATOM CANONIZER
	pcss_FlagStereoCenters(pcsg->size,OldParity,Graph);
	s = pcss_ReadSignature(SM);
	pcss_ResetFlags(Graph,SM);
	pcss_AddAtomSigGraph(pcsg, s, FALSE);
	// Sort children as in data structure Graph
	pcss_OrderChildrenAsInGraph(pcsg, Graph);
	// MODIFY ATOM TYPE FOR STEREOCENTERS
	pcss_GlobalStereoCenters(pcsg, Graph, Parity);
	for (i = 0; i < Graph.size; i++) {
		if (Parity[i] != 0) {
			if (OldParity[i] == 0) {
				if (Parity[i] == -1) // S
					strcat(Graph.MyAtom[i][0].potential_type, "@");
				if (Parity[i] == 1) // R
					strcat(Graph.MyAtom[i][0].potential_type, "@@");
			}
		}
	}
	pcss_FreeSigGraph(pcsg);
}

//  Main routine for computing and assigning the parities of stereobonds
void pcss_StereobondParities(int* BondParity,int* OldBondParity,struct GRAPH Graph, set_molecule_p SM) {
	signature_p s = NIL;
	int i;
	struct SigGraph* pcsg;

	pcss_AtomList(SM->molecule,Graph);
	pcsg = pcss_NewSigGraph(Graph.size, Graph.NumbEdges);
	// CALL THE CIP BOND CANONIZER
	pcss_FlagStereoBonds(pcsg->nbonds,OldBondParity,Graph);
	s = pcss_ReadSignature(SM);
	pcss_ResetFlags(Graph,SM);
	pcss_AddBondSigGraph(pcsg, s, FALSE);
	pcss_OrderChildrenBondAsInGraph(pcsg, Graph);
	pcss_GlobalStereoBonds(pcsg, Graph, BondParity);
	// Sort children as in data structure Graph
	// MODIFY ATOM TYPE FOR STEREOBONDS
	for (i = 0; i < Graph.NumbEdges; i++) {
		if (BondParity[i] != 0) {
			if (OldBondParity[i] == 0) {
				if (BondParity[i] == 1) { // cis
					strcat(Graph.MyAtom[Graph.head[i]][0].potential_type,"/");
					strcat(Graph.MyAtom[Graph.end[i]][0].potential_type,"/");
				} else if (BondParity[i] == -1) { // trans
					strcat(Graph.MyAtom[Graph.head[i]][0].potential_type,"\\");
					strcat(Graph.MyAtom[Graph.end[i]][0].potential_type,"\\");
				}
			}
		}
	}
	pcss_FreeSigGraph(pcsg);
}

// Allocate and init parities arrays
void pcss_InitParities(int **Parity,int **BondParity, int **OldParity, int **OldBondParity,struct GRAPH Graph) {
	*Parity = calloc(Graph.size, sizeof(int));
	*BondParity = calloc(Graph.NumbEdges, sizeof(int));
	*OldParity = calloc(Graph.size, sizeof(int));
	*OldBondParity = calloc(Graph.NumbEdges, sizeof(int));
	pcss_InitialParities(*Parity,*BondParity,Graph);
}

// Deallocate parities arrays
void pcss_FreeParities(int *Parity,int *BondParity, int *OldParity, int *OldBondParity) {
	free(Parity);
	free(BondParity);
	free(OldParity);
	free(OldBondParity);
}

// Initialize the parities
void pcss_InitialParities(int* Parity,int* BondParity,struct GRAPH Graph) {
	int i;
	for (i = 0; i < Graph.size; i++)
		if ((Graph.atom_identity[i].candidate) && (Graph.atom_identity[i].order!=5)) // Skip planar non-flagged candidates
			Parity[i] = 1;
	for (i = 0; i < Graph.NumbEdges; i++)
		if (Graph.bond_identity[i].align)
			BondParity[i] = 1;
}

// Update atom parities arrays
void pcss_UpdateAtomParity(int* Parity,int* OldParity,struct GRAPH Graph) {
	int i;
for (i = 0; i < Graph.size; i++)
	if (OldParity[i] == 0)
		OldParity[i] = Parity[i];
}

// Update bond parities arrays
void pcss_UpdateBondParity(int* BondParity,int* OldBondParity,struct GRAPH Graph) {
	int i;
	for (i = 0; i < Graph.NumbEdges; i++)
		if (OldBondParity[i] == 0)
			OldBondParity[i] = BondParity[i];
}


// Data structure constructor
// The flag is TRUE if the signature contains first the string otherwise only the invariants
struct SigGraph* pcss_NewSigGraph(int size, int nbonds) {

	struct SigGraph *pcsg = calloc(1,sizeof(struct SigGraph));
	pcsg->size = size;
	pcsg->nbonds = nbonds;
	pcsg->atom = NULL;
	pcsg->sig = NULL;
	pcsg->bond = NULL;
	pcsg->bsig = NULL;
	return(pcsg);
}

// Add atoms to the data structure from the signature
void pcss_AddAtomSigGraph(struct SigGraph *sg,signature_p s,t_bool SIGSTRING) {

	if (PCINFO) printf("ADDING ATOM PRIORITIES\n");

	sg->sig = calloc(sg->size,sizeof(struct SigInfo));
	sg->atom = calloc(sg->size,sizeof(struct AtomInfo));

	int i;
	for(i=0;i<sg->size;i++) {
		if (s[i].value>0) { // If a signature actually exists
			sg->sig[i].s = s[i].as[0].s;
			sg->sig[i].nlabels = HowManyLabels(sg->sig[i].s);
			sg->sig[i].labels = pcss_ParseSignature(sg->sig[i].nlabels,sg->sig[i].s,SIGSTRING);
			int root = RootAtom(sg->sig[i].s,SIGSTRING);
			pcss_GetAtomInfo(sg->size,sg->sig[i].labels,sg->sig[i].nlabels,&sg->atom[root]);
		}
	}
}
// Add bonds to the data structure from the signature
void pcss_AddBondSigGraph(struct SigGraph *sg, signature_p s,t_bool SIGSTRING) {

	if (PCINFO) printf("ADDING BOND PRIORITIES\n");

	sg->bsig = calloc(sg->nbonds,sizeof(struct SigInfo));
	sg->bond = calloc(sg->nbonds,sizeof(struct BondInfo));

	int i;
	for(i=0;i<sg->nbonds;i++) {
		if (s[i].value>0) { // If a signature actually exists
			sg->bsig[i].s = s[i].as[0].s;
			sg->bsig[i].nlabels = HowManyLabels(sg->bsig[i].s);
			sg->bsig[i].labels = pcss_ParseSignature(sg->bsig[i].nlabels,sg->bsig[i].s,SIGSTRING);

			sg->bond[i].head = HeadBond(sg->bsig[i].s,SIGSTRING);
			sg->bond[i].tail = TailBond(sg->bsig[i].s,SIGSTRING);

			sg->bond[i].branches_head = 0;
			sg->bond[i].branches_tail = 0;
			pcss_GetBondInfo(sg->size,sg->bsig[i].labels,sg->bsig[i].nlabels,&sg->bond[i]);
		}
	}
}

// Free allocation for data structure
void pcss_FreeSigGraph(struct SigGraph *sg) {
	int i;
	if (sg->atom==NULL) return;
	if (sg->atom!=NULL) {
		for(i=0;i<sg->size;i++) {
			free(sg->sig[i].labels);
			free(sg->atom[i].rootpriorities);
			free(sg->atom[i].children);
		}
	}
	if (sg->bond!=NULL) {
		for(i=0;i<sg->nbonds;i++) {
			free(sg->bsig[i].labels);
			free(sg->bond[i].rootpriorities);
		}
	}
	if (sg->sig!=NULL) free(sg->sig);
	if (sg->bsig!=NULL) free(sg->bsig);
	if (sg->atom!=NULL) free(sg->atom);
	if (sg->bond!=NULL) free(sg->bond);
}

// Set Diameter for computing bond signature (odd)
int pcss_DiameterBondSignature(int init_dim)
{
	if (init_dim%2 == 0) { // EVEN
		return(init_dim-1);
	} else { // ODD
		return(init_dim);
	}
}

// Set Diameter for computing atom signature (even)
int pcss_DiameterAtomSignature(int init_dim)
{
	if (init_dim%2 == 0) { // EVEN
		return(init_dim);
	} else { // ODD
		return(init_dim-1);
	}
}

// Order children atoms in SigGraph structure as in Graph structure
void pcss_OrderChildrenAsInGraph(struct SigGraph *sg,struct GRAPH Graph){
	int i,j;
	for(i=0;i<sg->size;i++) {
		for(j=0;j<sg->atom[i].branches;j++)
		sg->atom[i].children[j] = Graph.Neigh[i][j];
	}
}

// Order children bonds in SigGraph structure as in Graph structure
void pcss_OrderChildrenBondAsInGraph(struct SigGraph *sg,struct GRAPH Graph){
	int i,j;

	for(i=0;i<sg->nbonds;i++) {
		for(j=0;j<sg->bond[i].branches_head;j++)
				sg->bond[i].children_head[j] = Graph.Neigh[sg->bond[i].head][j];
		for(j=0;j<sg->bond[i].branches_tail;j++)
				sg->bond[i].children_tail[j] = Graph.Neigh[sg->bond[i].tail][j];
	}
}

//////////// PARITIES AND SYMMETRY SECTION

// Check if parity vector has changed
t_bool pcss_ParitiesChange(int size, int* oldparities,int* newparities) {
	int i;
	for(i=0;i<size;i++)
		if ((oldparities[i]==0) && (newparities[i]!=0))  return(TRUE);
	return FALSE;
}

// Check if atom parities have changed
t_bool pcss_AtomParitiesChange(int size,int* currentparities,struct GRAPH Graph) {
	int i;
	for(i=0;i<size;i++) {
		if ((currentparities[i]==0) && (Graph.atom_identity[i].candidate!=0) && (Graph.atom_identity[i].order!=5))  {
			return(TRUE);
		}
	}
	return(FALSE);
}

// Check if bond parities have changed
t_bool pcss_BondParitiesChange(int nbonds, int* currentbondparities,struct GRAPH Graph) {
	int i;
	for(i=0;i<nbonds;i++)
		if ((currentbondparities[i]==0) && (Graph.bond_identity[i].align != 0))  {
			return(TRUE);
		}
	return(FALSE);
}

// Flag those stereocenter candidates whose parity still have not been assigned
void pcss_FlagStereoCenters(int size, int* newparities, struct GRAPH Graph) {
	int i;
	for(i=0;i<size;i++) {
		if ((newparities[i]==0) && (Graph.atom_identity[i].candidate!=0) && (Graph.atom_identity[i].order!=5))  {
			strcpy(Graph.MyAtom[i][0].comment, "1");
		} else {
			Graph.MyAtom[i][0].comment[0] = '0';
		}
	}
}

// Flag those stereobond candidates whose parity still have not been assigned
void pcss_FlagStereoBonds(int nbonds, int* newbondparities, struct GRAPH Graph) {
	int i;
	for(i=0;i<Graph.size;i++) {
				strcpy(Graph.MyAtom[i][0].comment, "0");
	}
	for(i=0;i<nbonds;i++)
		if ((newbondparities[i]==0) && (Graph.bond_identity[i].align != 0))  {
			strcpy(Graph.MyAtom[Graph.head[i]][0].comment, "1");
			strcpy(Graph.MyAtom[Graph.end[i]][0].comment, "1");
		}
}

// Reset flags in atoms
void pcss_ResetFlags(struct GRAPH Graph,set_molecule_p SM) {
	int i;
	for(i=0;i<Graph.size;i++)
		strcpy(Graph.MyAtom[i][0].comment, "0");
}

//Universal orderer of sequences (from larger to smaller is +1 (of course with all the cyclic permutations))
int cip_InvDetOrderVec(t_real *Vec, int size){
        int count = 1;
        int i,j, aux;
        for(i=0;i<size;i++){
                aux = Vec[i];
                for(j=i+1;j<size;j++)
                        if(Vec[j] < aux)
                                count++;
        }
        count += size; //To compare parity of different length so that 123, 0123 have the same parity but 321 and 3210 do not.
        if (size==4) count -= 1;

        count = 2 * (count % 2) -1;
        return count;
}

// Check vector parity
int CheckSameWeight(t_real* Vec,int size) {
	int i,j;
	int diff = 1;
	for (i=0;i<size;i++) {
		for (j=i+1;j<size;j++) {
			if (Vec[i]==Vec[j]) diff = 0;
		}
	}
	return(diff);
}

// Returns the parity of the stereocenter in function of the Priorities of the children and the geometry of the atoms
int pcss_SterocenterParity( struct GRAPH Graph,struct SigGraph *sg,int root,int atom){

	int i;
	int nbranches=sg->atom[atom].branches;
	if (nbranches<3) return(0);

	// PC March 2013
	if ((nbranches==4) && (Graph.atom_identity[atom].order>10)){
		int ord;
		int K=10;
		int j =0;
		t_real* auxvec = calloc(4,sizeof(t_real));
		for(i=0;i<nbranches;i++) {
			if (sg->atom[root].rootpriorities[sg->atom[atom].children[i]]==4) {
				ord = Graph.atom_identity[atom].order % K;
				ord = ord / (K/10);
				if (ord == 1) {
					ord = 1;
				} else if (ord==2){
					ord = -1;
				}
				continue;
			}
			auxvec[j] = sg->atom[root].rootpriorities[sg->atom[atom].children[i]];
			if (auxvec[j]==0) return(0); // Some branch is missing
			K *= 10;
			j++;
		}
		int State =  -1*ord * cip_InvDetOrderVec(auxvec,3) *CheckSameWeight(auxvec,3);
		free(auxvec);
		return(State);


	}

	t_real* auxvec = calloc(nbranches,sizeof(t_real));
	for(i=0;i<nbranches;i++) {
		auxvec[i] = sg->atom[root].rootpriorities[sg->atom[atom].children[i]];
		if (auxvec[i]==0) return(0); // Some branch is missing
	}

	int State =  -1*Graph.atom_identity[atom].order * cip_InvDetOrderVec(auxvec,nbranches) *CheckSameWeight(auxvec,nbranches);
	free(auxvec);
	return(State);
}



////// SIGNATURES PARSING SECTION

// Returns the signature of each atom
// Based on DF's function
signature_p	pcss_ReadSignature(SM)
/******************************************************************************
 Write scan files
 ******************************************************************************/
set_molecule_p	SM;
{
	set_molecule_p  S;
	molecule_p	molecule;
	signature_p	s;

	if (SM == NIL) return(NIL);
	if (SM->molecule == NIL) return(NIL);

	S = SM;
	while (S) {
		int	h;
		molecule = S->molecule;
		if ((molecule == NIL) || (molecule == (molecule_p)ERROR)) break;
		s = NIL;
		for (h = sisi_inq_dimension_min(); h <= sisi_inq_dimension_max(); h++) {
			sisi_set_dimension(h);
			s = sisc_signature(molecule->HEAP,molecule,NIL,NIL,NIL);
		}
		S = S->succ;
	}
	return(s);
}


// Find out how many labels does contain the signature
// for allocation purposes, works with both the signature or the invariants
int HowManyLabels(char *sigcursor){
	int natoms = 0;
	while( (sigcursor[0] != '<' ) && (sigcursor[0] != '\0'  ) ) {
		if (sigcursor[0] == '[') natoms++;
		sigcursor++;
	}
	return(natoms);
}

// Returns label of next atom
char* NextAtomLabel(char* sigcursor, struct LabelPosition* position, int depth){

	position->depth = depth;
	position->color = -1; // no color
	sigcursor++;
	position->start = sigcursor;
	while (sigcursor[0]!=']') {
		sigcursor++;
		if (sigcursor[0] == ',') {
			position->color = atoi(++sigcursor);
		}
	}
	position->end = sigcursor;
	return(sigcursor);
}

// Read atom weight in the signature
char* NextAtomWeight(char* sigcursor, struct LabelPosition* position){

	position->id = atoi(++sigcursor);
	while((sigcursor[0]!=',') && (sigcursor!='\0')) sigcursor++;
	if (sigcursor[0] == ',') position->weight = (t_real) atof(++sigcursor);
	return(sigcursor);
}

// Returns the priorities of the children of the given atom
int PrioritiesChildren(t_real* Priorities,int size,int nlabels,struct LabelPosition* AtomLabels,int atom,t_bool RESET) {
	// Find the first instance of the atom
	int i, NextLevel, depth = 0;
	for (i=0;i<size;i++) Priorities[i] = 0.0;
	while( (Priorities[atom] == 0) && NextLevel ) {
		NextLevel = BranchPriorities(Priorities,size,nlabels,AtomLabels,depth,RESET);
		depth++;
	}
	return(BranchPriorities(Priorities,size,nlabels,AtomLabels,depth,RESET));
}


// Find out which atom is the root of this atomic signature
int RootAtom(char *sigcursor,t_bool SIGSTRING) {
	if (SIGSTRING)
		while(sigcursor[0]!='>') sigcursor++;
	while(sigcursor[0]!='[') sigcursor++;
	return(atoi(++sigcursor));
}

// Find out which bond head
int HeadBond(char *sigcursor,t_bool SIGSTRING) {
	if (SIGSTRING)
		while(sigcursor[0]!='>') sigcursor++;
	while(sigcursor[0]!='[') sigcursor++;
	int head = atoi(++sigcursor);
	return(head);
}
// Find out which bond tail
int TailBond(char *sigcursor,t_bool SIGSTRING) {
	if (SIGSTRING)
		while(sigcursor[0]!='>') sigcursor++;
	int depth = 0;
	while(sigcursor[0]!='[') sigcursor++;
	sigcursor++;
	while((sigcursor[0]!='[') || (depth!=0)) {
		if (sigcursor[0]=='(') {
			depth++;
		} else if (sigcursor[0]==')') {
			depth--;
		} else if (sigcursor[0]=='\0') {
			break;
		}
		sigcursor++;
	}
	int tail = atoi(++sigcursor);
	return(tail);
}


///// GRAPH SECTION

// Returns number of children
int pcss_NumberChildren(struct LabelPosition *labels,int nlabels) {
	int i;
	int nchild = 0;
	for (i=0;i<nlabels;i++)
		if (labels[i].depth == 1) nchild++;
	return(nchild);
}

// This code is only intended for bond signatures
int pcss_NumberChildrenParent(struct LabelPosition *labels,int nlabels,int parent) {
	int i;
	t_bool BRANCH = FALSE;
	int nchild = 0;
	for (i=0;i<nlabels;i++) {
		if (labels[i].depth == 0) {
					if (labels[i].id == parent) BRANCH = TRUE;
					else {
							nchild++; // The other side of the bond
							BRANCH = FALSE;
					}
				}
		if (BRANCH &&  (labels[i].depth == 1)) nchild++;
	}
	return(nchild);
}

// Compare children
int comparechildren (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

// Get children
int* pcss_GetChildren(struct LabelPosition *labels,int nlabels,int nchild) {
	int i;
	int* children = calloc(nchild,sizeof(int));
	int n = 0;
	for (i=0;i<nlabels;i++)
		if (labels[i].depth == 1) children[n++] = labels[i].id;
	qsort(children,n,sizeof(int),comparechildren);
	return(children);
}

// This code is only intended for bond signatures
int* pcss_GetChildrenParent(struct LabelPosition *labels,int nlabels,int nchild,int parent) {
	int i;
	t_bool BRANCH = FALSE;
	int* children = calloc(nchild,sizeof(int));
	int n = 0;
	for (i=0;i<nlabels;i++) {
		if (labels[i].depth == 0) {
			if (labels[i].id == parent) BRANCH = TRUE;
			else {
				children[n++] = labels[i].id; // Add the other side of the bond
				BRANCH = FALSE;
			}
		}
		if (BRANCH && (labels[i].depth == 1)) children[n++] = labels[i].id;
	}
	qsort(children,n,sizeof(int),comparechildren);
	return(children);
}


// Get Priorities from the labels
void pcss_GetAtomInfo(int size,struct LabelPosition *labels,int nlabels, struct AtomInfo *atom) {
	int i;
	t_real* priorities = calloc(size,sizeof(t_real));
	for (i=0;i<nlabels;i++) {
		if (labels[i].depth==1) { // We want priorities of the children layer
			priorities[labels[i].id] = labels[i].weight;
		}
	}
	atom->rootpriorities = priorities;
	atom->branches = pcss_NumberChildren(labels,nlabels);
	atom->children =  pcss_GetChildren(labels,nlabels,atom->branches);
}

// Get Priorities from the labels
void pcss_GetBondInfo(int size,struct LabelPosition *labels,int nlabels, struct BondInfo *bond) {
	int i;
	t_real* priorities = calloc(size,sizeof(t_real));
	for (i=0;i<nlabels;i++) {
		priorities[labels[i].id] = labels[i].weight;
	}
	bond->rootpriorities = priorities;
	bond->branches_head = pcss_NumberChildrenParent(labels,nlabels,bond->head);
	bond->branches_tail = pcss_NumberChildrenParent(labels,nlabels,bond->tail);
	bond->children_head =  pcss_GetChildrenParent(labels,nlabels,bond->branches_head,bond->head);
	bond->children_tail =  pcss_GetChildrenParent(labels,nlabels,bond->branches_tail,bond->tail);

}

// NOTE: In principle, the priorities will be given in function of the depth, because we can have the same atom
// at different layers, and we want to know the priorities for a specific layer
int BranchPriorities(t_real* Priorities,int size,int nlabels,struct LabelPosition* AtomLabels,int depth,t_bool RESET) {
	int i;
	int natoms = 0;
	if (RESET) {
		for (i=0;i<size;i++) Priorities[i] = 0.0;
	}
	for(i=0;i<nlabels;i++){
		if (AtomLabels[i].depth == depth) {
			Priorities[AtomLabels[i].id] = AtomLabels[i].weight;
			natoms++;
		}
	}
	return(natoms);
}

// Parse the signature with atom weights
struct LabelPosition* pcss_ParseSignature(int nlabels,char* s,t_bool SIGSTRING) {

	t_bool WEIGHTS = FALSE;
	t_bool DBOND = FALSE;
	char* sigcursor;
	sigcursor = s;

	struct LabelPosition* AtomLabels = calloc(nlabels,sizeof(struct LabelPosition));
	int labelpos = 0;
	int depth = 0;
	while(sigcursor[0] != '\0') {
		switch(sigcursor[0])
		{
			case '(':
				depth++;
				break;
			case ')':
				depth--;
				break;
			case '>':
				WEIGHTS = TRUE;
				labelpos = 0;
				break;
			case '=':
				DBOND = TRUE;
				break;
			case '[':
				if (SIGSTRING ==  TRUE) {
					if (WEIGHTS == FALSE) {
						sigcursor = NextAtomLabel(sigcursor,(struct LabelPosition*) &AtomLabels[labelpos],depth);
					} else {
						sigcursor = NextAtomWeight(sigcursor,(struct LabelPosition*) &AtomLabels[labelpos]);
					}
				} else {
					NextAtomLabel(sigcursor,(struct LabelPosition*) &AtomLabels[labelpos],depth);
					sigcursor = NextAtomWeight(sigcursor,(struct LabelPosition*) &AtomLabels[labelpos]);
				}
				if (DBOND==TRUE) {
					AtomLabels[labelpos].dbond = 1;
					DBOND = FALSE;
				}
				labelpos++;
				break;
		}
		sigcursor++;
	}
	return(AtomLabels);

}

// ADDING STEREO LABELS TO THE SIGNATURE (LOCAL ASSIGNMENT)

// This function is used for Local stereo assignment
// Once the parity of the atoms in the atomic signature are determined, it recreates the signature with the parity symbols
char* pcss_AddSignatureAtomLabels(struct SigGraph *sg,int* parities,int signature)  {

	int i;
	struct LabelPosition* AtomLabels = sg->sig[signature].labels;
	char* sig = calloc(strlen(sg->sig[signature].s),2*sizeof(char)); // Length of signature with priorities is enough
	char* sigcursor = sg->sig[signature].s;

	// Add parities
	t_bool WEIGHTS = FALSE;
	i = 0; // first label
	while((sigcursor[0] != '\0')&&(!WEIGHTS)) {
		switch(sigcursor[0])
		{
			case '<':
				WEIGHTS = TRUE;
				break;
			case '[':
				strncat(sig,sigcursor,1);
				char* labelcursor;
				for(labelcursor=AtomLabels[i].start;labelcursor!=AtomLabels[i].end;labelcursor++) {
					if (labelcursor[0]==',') break; //Colored atom
					strncat(sig,labelcursor,1);
				}
				if (parities[AtomLabels[i].id]) {
					if (parities[AtomLabels[i].id]==1) {
						strncat(sig,"@",1);
					} else {
						strncat(sig,"@@",2);
					}
				}
				// Continue if it's a colored atom
				for(;labelcursor!=AtomLabels[i].end;labelcursor++) {
					strncat(sig,labelcursor,1);
				}
				i++;
				while ( (sigcursor[0]!='\0')&&(sigcursor[1]!=']') ) sigcursor++;
				break;
			default:
				strncat(sig,sigcursor,1);
				break;
		}
		if (WEIGHTS) break;
		sigcursor++;
	}
	strcat(sig,sigcursor); // Add again the weights
	free(sg->sig[signature].s);
	sg->sig[signature].s = sig;
	return(sig);
}

// This function is used for Local stereo bond assignment
// Once the parity of the bonds in the atomic signature are determined, it recreates the signature with the bond parity symbols
char* pcss_AddSignatureBondLabels(struct SigGraph *sg,struct GRAPH Graph, int* parities,int signature)  {

	int i,j;
	struct LabelPosition* AtomLabels = sg->sig[signature].labels;
	char* sig = calloc(strlen(sg->sig[signature].s),2*sizeof(char)); // Length of signature with priorities is enough
	char* sigcursor = sg->sig[signature].s;
	int *newatomlabel = calloc(sg->sig[signature].nlabels,sizeof(int)); // Actually there might be less atoms than labels occurrences in the signature


	for(j=0;j<Graph.NumbEdges;j++)
		if (parities[j]!=0) {
			if (parities[j]==1) {
				newatomlabel[Graph.head[j]] = 1;
				newatomlabel[Graph.end[j]] = -1;
		}
		if (parities[j]==-1) {
				newatomlabel[Graph.head[j]] = 1;
				newatomlabel[Graph.end[j]] = 1;
		}
	}
	// Add parities
	t_bool WEIGHTS = FALSE;
	i = 0; // first label
	while((sigcursor[0] != '\0')&&(!WEIGHTS)) {
		switch(sigcursor[0])
		{
			case '<':
				WEIGHTS = TRUE;
				break;
			case '[':
				strncat(sig,sigcursor,1);
				char* labelcursor;
				for(labelcursor=AtomLabels[i].start;labelcursor!=AtomLabels[i].end;labelcursor++) {
					if (labelcursor[0]==',') break; //Colored atom
					strncat(sig,labelcursor,1);
				}
				if (newatomlabel[AtomLabels[i].id]==1) {
					strcat(sig,"\\");
				}
				if (newatomlabel[AtomLabels[i].id]==-1) {
					strcat(sig,"/");
				}
				// Continue if it's a colored atom
				for(;labelcursor!=AtomLabels[i].end;labelcursor++) {
					strncat(sig,labelcursor,1);
				}
				i++;
				while ( (sigcursor[0]!='\0')&&(sigcursor[1]!=']') ) sigcursor++;
				break;
			default:
				strncat(sig,sigcursor,1);
				break;
		}
		if (WEIGHTS) break;
		sigcursor++;
	}

	strcat(sig,sigcursor); // Add again the weights
	free(sg->sig[signature].s);
	sg->sig[signature].s = sig;
	free(newatomlabel);
	return(sig);
}


// STEREOCENTERS AND STEREOBONDS SECTION


// Check for stereocenters
void pcss_GlobalStereoCenters(struct SigGraph *sg,struct GRAPH Graph,int *Parity) {
	int i;
	for(i=0;i<sg->size;i++)
		if(Graph.atom_identity[i].candidate!=0) {
			int state = pcss_SterocenterParity(Graph,sg,i,i);
			Parity[i] = state;
			if (PCINFO) {
				if (state!=0) pcss_PrintStereocenterInfo(sg,i,state,Graph);
			}
		}
}


// Find the bond with the priorities info
int pcss_FindBond(struct SigGraph *sg,int head, int tail) {
	int i;
	for (i=0;i<sg->nbonds;i++){
		if ((sg->bond[i].head==head)&&(sg->bond[i].tail==tail))
			return(i);
		if ((sg->bond[i].head==tail)&&(sg->bond[i].tail==head))
			return(i);
	}
	return(-1);
}

void pcss_ExchangeBond(struct SigGraph *sg, int nbond) {
	int t1,*t2;
	t1 = sg->bond[nbond].head;
	sg->bond[nbond].head = sg->bond[nbond].tail;
	sg->bond[nbond].tail = t1;
	t1 = sg->bond[nbond].branches_head;
	sg->bond[nbond].branches_head = sg->bond[nbond].branches_tail;
	sg->bond[nbond].branches_tail = t1;
	t2 = sg->bond[nbond].children_head;
	sg->bond[nbond].children_head = sg->bond[nbond].children_tail;
	sg->bond[nbond].children_tail = t2;
}

// This code is based on Davide's, but only the essential has been kept
int pcss_StereobondParity(struct GRAPH Graph,struct SigGraph *sg,int bond){
	int Parity = 1;

	int head = Graph.head[bond];
	int end = Graph.end[bond];
	int sgbond;
	if ((sgbond=pcss_FindBond(sg,head,end))<0) // Bond not found
	          return(Parity);
	if (end == sg->bond[sgbond].head) {
		pcss_ExchangeBond(sg,sgbond);
//		Parity *= -1;
//		head = Graph.end[bond];
//		end = Graph.head[bond];
	}

	t_real* Priorities_head = sg->bond[sgbond].rootpriorities;
	t_real* Priorities_end = sg->bond[sgbond].rootpriorities;

	int NChildHead = sg->bond[sgbond].branches_head;
	int NChildEnd = sg->bond[sgbond].branches_tail;
	int* ChildrenHead = sg->bond[sgbond].children_head;
	int* ChildrenEnd = sg->bond[sgbond].children_tail;

	// Set to FALSE if it is a fragment
	t_bool BRANCHESOK = TRUE;


			//We calculate first the parity of the site on the top, we distinguish btw there are two
			//other edges (beside the =C)
			//Consider the first not end neighbours of head


	if( ((NChildHead<2) || (NChildEnd<2)) ) { //||((NChildHead<=2)&&(NChildHead<=2)) ){
		BRANCHESOK = FALSE;
		// PC: Fragment, remove parity
	}else{

		int atom1, atom2;
		int nbranch;
		int TestDegeneration = 1;
		// Determine orientation at the head
		atom1 = 0; atom2 = 0;
		if (NChildHead == 3) {
			if(ChildrenHead[atom2]== end) atom2++;
			atom1 = atom2++;
			if(ChildrenHead[atom2] == end) atom2++;

			if (Priorities_head[ChildrenHead[atom1]] < Priorities_head[ChildrenHead[atom2]])
				Parity *= -1;
			TestDegeneration *= (Priorities_head[ChildrenHead[atom1]]
					- Priorities_head[ChildrenHead[atom2]]);


		} else {
			atom1 = 0;
			if(ChildrenHead[atom1]== end) atom1++;
		}
		nbranch = 0;
		if (Priorities_head[ChildrenHead[atom1]]!=0) nbranch++;
		if ((atom2>0)&&(Priorities_head[ChildrenHead[atom2]]!=0)) nbranch++;
		if (nbranch != (NChildHead-1)) BRANCHESOK = FALSE;


		atom1 = 0; atom2 = 0;
		// Determine orientation at the end
		if (NChildEnd == 3) {

			if(ChildrenEnd[atom2]== head) atom2++;
			atom1 = atom2++;
			if(ChildrenEnd[atom2] == head) atom2++;

			if (Priorities_end[ChildrenEnd[atom1]] < Priorities_end[ChildrenEnd[atom2]])
				Parity *= -1;

			TestDegeneration *= (Priorities_end[ChildrenEnd[atom1]]
					- Priorities_end[ChildrenEnd[atom2]]);
		} else {
			atom1 = 0;
			if(ChildrenEnd[atom1]== head) atom1++;
		}

		nbranch = 0;
		if (Priorities_end[ChildrenEnd[atom1]]!=0) nbranch++;
		if ((atom2>0)&&(Priorities_end[ChildrenEnd[atom2]]!=0)) nbranch++;
		if (nbranch != (NChildEnd-1)) BRANCHESOK = FALSE;

		if(Graph.bond_identity[bond].align != 0){//condition immutated if end has not two identical branches.
			Parity *= Graph.bond_identity[bond].align;
		}
		if (TestDegeneration == 0) {//branches are not different
			BRANCHESOK = FALSE;
		}


	}
	if ( BRANCHESOK == FALSE ) {
		// PC: Fragment, remove parity
		Graph.bond_identity[bond].align = 0;
		Graph.bond_identity[bond].candidate = 0;
		Parity = 0;
	}

	return(Parity);
}

// Check for stereobonds
void pcss_GlobalStereoBonds(struct SigGraph *sg,struct GRAPH Graph,int *StateLink) {
	int i;
	for (i = 0; i < sg->nbonds; i++)
		if (Graph.bond_identity[i].align != 0) {
			// Get priorities of head's children in signature rooted at head
			int state = pcss_StereobondParity(Graph, sg, i);
			StateLink[i] = state;
			if (PCINFO) {
				if (state != 0)
					pcss_PrintStereobondInfo(sg,Graph,i, state);
			}
		}
}


void pcss_LocalStereoCenters(struct SigGraph *sg,int* ListCandidates,struct GRAPH Graph) {

	// First we root each atom in order to store in RootPriorities the list of their children
	int i,j,z;
	int* Parity = calloc(sg->size,sizeof(int));

	// CHECK FOR STEREOCENTERS for each signature and each atom
	for(j=0;j<sg->size;j++) {// For each signature
		z = sg->sig[j].labels[0].id; // Find out who is the root
		for(i=0;i<sg->size;i++) {// For each atom
			Parity[i] = 0; // Reset
			if(Graph.atom_identity[i].candidate!=0) { // If the atom is a candidate
				//Now we compute the parity of atom i in signature j (rooted in atom z)
				Parity[i] = pcss_SterocenterParity(Graph,sg,z,i);
			}
		}
		// Now we update the signature with the parity labels
		pcss_AddSignatureAtomLabels(sg,Parity,j);
	}
	free(Parity);
}



void pcss_LocalStereoBonds(struct SigGraph *sg,int* ListCandidates,struct GRAPH Graph) {

	int i,j;
	int* Parities = calloc(Graph.NumbEdges,sizeof(int));

	// CHECK FOR STEREOBONDS for each signature and each bond
	for(j=0;j<sg->size;j++) {// For each signature
		for(i=0;i<Graph.NumbEdges;i++) { // For each bond
			Parities[i] = 0; // Reset
			if(ListCandidates[i] != 0){
				Parities[i] = pcss_StereobondParity(Graph,sg,i);
			}
		}
		// Now we update the signature with the parity labels
		pcss_AddSignatureBondLabels(sg,Graph,Parities,j);

	}
	free(Parities);
}

void pcss_UpdateSignature(struct SigGraph *sg, signature_p s) {
	int i;
	for(i=0;i<sg->size;i++) {
		s[i].as[0].s = sg->sig[i].s;
	}
}

void pcss_PrintStereobondInfo(struct SigGraph *sg,	struct GRAPH Graph,int i,int state) {
		int SWITCH = FALSE;
        int dir = 'E';
        if (state == 1) dir = 'Z';
        int x = pcss_FindBond(sg,Graph.head[i],Graph.end[i]);
        if (x<0) return;
        if (Graph.end[i] == sg->bond[x].head) SWITCH = TRUE;

        int na,nb;
        int *a,*b;
        if (!SWITCH) {
        		na = sg->bond[x].branches_head;
        		a = sg->bond[x].children_head;
        		nb = sg->bond[x].branches_tail;
        		b = sg->bond[x].children_tail;
        } else {
        		na = sg->bond[x].branches_tail;
    			a = sg->bond[x].children_tail;
        		nb = sg->bond[x].branches_head;
        		b = sg->bond[x].children_head;
        }

        printf("BOND %d HEAD %d END %d PARITY %c HEAD PRIORITIES: ", i,
                        Graph.head[i] + 1, Graph.end[i] + 1, dir);
        int j;
        for (j = 0; j < na; j++) {
				if (a[j] != Graph.end[i]) {
                        printf("%d:%.2f ", a[j] + 1, sg->bond[x].rootpriorities[a[j]]);
				}
        }
        printf("END PRIORITIES: ");
        for (j = 0; j < nb; j++) {
				if (b[j] != Graph.head[i]) {
                        printf("%d:%.2f ", b[j] + 1, sg->bond[x].rootpriorities[b[j]]);
				}
        }
        printf("\n");
}

void pcss_PrintStereocenterInfo(struct SigGraph *sg, int i,int state,struct GRAPH Graph) {
        int dir = 'S';
        if (state == 1)
                dir = 'R';
        printf("ATOM  %d PARITY %c ORDER %d CHILDREN PRIORITIES: ", i + 1, dir,Graph.atom_identity[i].order);
        int j;
        printf("%d BRANCHES ",sg->atom[i].branches);
        for (j = 0; j < sg->atom[i].branches; j++) {
                        printf("%d:%.2f ", sg->atom[i].children[j] + 1, sg->atom[i].rootpriorities[sg->atom[i].children[j]]);
        }
        printf("\n");
}


// Get AtomList
DEFINE 	t_err   	pcss_AssignAtom(a,MyAtom)
/***********************************************************
 Assign atoms to the array according to their IDs
 **********************************************************/
atom_p	a;
atom_p *MyAtom;
{
	MyAtom[a->ID] = a;
	return(OK);
	}

t_err   	pcss_AtomList(m,Graph)
/***********************************************************
 Populate array with atoms ordered by their IDs
 **********************************************************/
molecule_p	m;
struct GRAPH Graph;
{

	if (m == NIL) return(ERROR);
	if (m->size < 1) return(ERROR);
	dast_all_atom(m,NIL,pcss_AssignAtom,FALSE,(t_pointer) Graph.MyAtom);

	return(OK);
	}





