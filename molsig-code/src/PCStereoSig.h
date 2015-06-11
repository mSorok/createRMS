/*
 * PCStereoSig.h
 *
 *  Created on: March 2012
 *      Author: Carbonell
 */

#ifndef PCSTEREOSIG_H_
#define PCSTEREOSIG_H_

// Print debug info
EXTERN t_bool PCINFO;

void pcss_InitParities(int **Parity,int **BondParity, int **OldParity, int **OldBondParity,struct GRAPH Graph);
void pcss_FreeParities(int *Parity,int *BondParity, int *OldParity, int *OldBondParity);

void pcss_UpdateAtomParity(int* BondParity,int* OldBondParity,struct GRAPH Graph);
void pcss_UpdateBondParity(int* BondParity,int* OldBondParity,struct GRAPH Graph);
void pcss_StereocenterParities(int* Parity,int* OldParity,struct GRAPH Graph, set_molecule_p SM);
void pcss_StereobondParities(int* BondParity,int* OldBondParity,struct GRAPH Graph, set_molecule_p SM);
void pcss_ResetFlags(struct GRAPH Graph, set_molecule_p SM);

#endif /* PCSTEREOSIG_H_ */
