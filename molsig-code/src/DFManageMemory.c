/*
 *  ManageMemory.h
 *  Functions to manage the memory of Graphs, create bonds...
 *
 *  Created by Davide Fichera on 01/10.
 *
 */
#include "DFgeneral.h"

//Given the number of atoms of the molecule and the number of edges 
//allocate the space for the Graph describing the molecule. 
//Many pointers are given but many other have to be fixed once we know the 
//coordination of atoms.
int dfmm_CreateGRAPHMemory(int na, int nb, struct GRAPH *MyGraph){
	//ALLOC SPACE
	MyGraph[0].atom_identity = calloc(na , sizeof (struct ATOM_IDENTITY));
	MyGraph[0].bond_identity = calloc(nb , sizeof (struct BOND_IDENTITY));
	MyGraph[0].Neigh = calloc(na + na,sizeof(int*));
	MyGraph[0].Edge = MyGraph[0].Neigh + na;//na
	MyGraph[0].atom_identity[0].idNeigh = calloc( na + 8 * nb,sizeof(int));//2 * nb
	MyGraph[0].NNeigh = MyGraph[0].atom_identity[0].idNeigh + 2 * nb;//na
	MyGraph[0].end = MyGraph[0].NNeigh + na;//nb
	MyGraph[0].head = MyGraph[0].end + nb;//nb
	MyGraph[0].Neigh[0] = MyGraph[0].head + nb;//2 * nb
	MyGraph[0].Edge[0] = MyGraph[0].Neigh[0] + 2 * nb;//2 * nb
	MyGraph[0].NumbBond = 0;//Useful to construct the Graph
	{int i;
		for(i=0;i<na;i++)
			MyGraph[0].NNeigh[i] = 0;
	}
	
	MyGraph[0].size = na;
	MyGraph[0].NumbEdges = nb;
	MyGraph[0].NumbDirectedBonds = 2 * nb;
	
	MyGraph[0].MyAtom =  (atom_p *)calloc(na,sizeof(atom_p));  // PC Lookup table with the DAG

	return 0;
}

//Free what was allocated by CreateGRAPHMemory
int dfmm_DeleteGraphMemory(struct GRAPH Graph){
	if (Graph.size != 0) {
		free(Graph.atom_identity[0].idNeigh);
		free(Graph.Neigh);
		free(Graph.bond_identity);
		free(Graph.atom_identity);
		free(Graph.MyAtom);
	}
	return 0;
}

//Write the identity of the i-th atom as given by the mol file: its coordinates xyz, its name, its mass, 
//its charge, its stereo.
int dfmm_MakeIdentity(struct ATOM_IDENTITY *ai, int i,double x, double y, double z, char *name, int mass, int charge, int stereo){
	ai[i].candidate = 0;//1 if it is a candidate to be stereocenter
	ai[i].ordered = 0; //1 if its neighs have been ordered
	ai[i].order = -1; //the tag of the ordering 4vector
	ai[i].charName[0] = name[0];
	ai[i].charName[1] = name[1];
	ai[i].charName[2] = name[2];
	ai[i].Name = Name2Number_atom(name);
	ai[i].Coord.x =x ;
	ai[i].Coord.y = y;
	ai[i].Coord.z = z;
	ai[i].mass = mass;
	ai[i].charge = charge;
	ai[i].stereo = stereo;
	return 0;
}

//This function has to be called for each edge. Then we will have a description of the graph 
//as a directed one and we will get the right coordination for atoms (NNeigh). The identity of 
//the bond is stored into bond_identity.
int dfmm_AddBond(int i, int j, struct GRAPH *Graph,int order, int stereo){
	int Aux = Graph[0].NumbBond;
	Graph[0].head[Aux] = j - 1;
	Graph[0].end[Aux] = i - 1;
	Graph[0].NNeigh[i - 1]++;
	Graph[0].NNeigh[j - 1]++;
	Graph[0].NumbBond++;
	Graph[0].bond_identity[Aux].order = order; 
	Graph[0].bond_identity[Aux].order_mod = order; 
	if(order == 4)
		Graph[0].bond_identity[Aux].order_mod = 1; //Sometimes we want to access to the number of shared electrons.
	Graph[0].bond_identity[Aux].stereo = stereo; 
	Graph[0].bond_identity[Aux].candidate = 0;
	Graph[0].bond_identity[Aux].align = 0;
	return 0;	
}

//This function has to be called only when NNeigh rightly describe the coordination of atoms as it is given by 
//the info in Graph[].head .end. Pointers are setted and Neigh[][] are correctly filled.
int dfmm_MakeNeigh(struct GRAPH *Graph){
	int i;
	int size = Graph[0].size;
	int *AuxNNeigh = calloc(size, sizeof(int));//they have to be setted to zero
	for(i=1;i<size;i++){
		Graph[0].atom_identity[i].idNeigh = Graph[0].atom_identity[i-1].idNeigh + Graph[0].NNeigh[i-1];  
		Graph[0].Neigh[i] = Graph[0].Neigh[i-1] + Graph[0].NNeigh[i-1];
		Graph[0].Edge[i] = Graph[0].Edge[i-1]+Graph[0].NNeigh[i-1];
	}
	int NumbEdges = Graph[0].NumbEdges;
	int Auxhead, Auxend;
	for(i=0;i<NumbEdges;i++){
		Auxhead = Graph[0].head[i];
		Auxend = Graph[0].end[i];		
		Graph[0].Neigh[Auxhead][AuxNNeigh[Auxhead]] = Auxend;
		Graph[0].Neigh[Auxend][AuxNNeigh[Auxend]] = Auxhead;
		Graph[0].Edge[Auxhead][AuxNNeigh[Auxhead]] = i;
		Graph[0].Edge[Auxend][AuxNNeigh[Auxend]] = i;
		AuxNNeigh[Auxhead]++;
		AuxNNeigh[Auxend]++;
	}
	for(i=0;i<size;i++){
		Graph[0].atom_identity[i].tag = i;
		Graph[0].atom_identity[i].NumbNeigh = Graph[0].NNeigh[i];
	}
	free(AuxNNeigh);
	return 0;
}


