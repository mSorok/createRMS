/*
 *  ManageMemory.h
 *  
 *
 *  Created by Davide Fichera on 06/01/10.
 *
 */

#include "DFgeneral.h"
int dfmm_CreateGRAPHMemory(int na, int nb, struct GRAPH *MyGraph);
int dfmm_DeleteGraphMemory(struct GRAPH Graph);
int dfmm_MakeIdentity(struct ATOM_IDENTITY *ai, int i,double x, double y, double z, char *name, int mass, int charge, int stereo);
int dfmm_AddBond(int i, int j, struct GRAPH *Graph,int order, int stereo);
int dfmm_MakeNeigh(struct GRAPH *Graph);
enum atomic_number Name2Number_atom(char* atom);









