/*
 *  FindGeometry.h
 *  
 *
 *  Created by Davide Fichera on 06/01/10.
 *
 */
#ifndef FINDGEOMETRY_H
#define FINDGEOMETRY_H



int ThreePlusH(struct GRAPH Graph, int i, struct ORDER *Order);
int dffg_StereoBonds(struct GRAPH Graph);
int dffg_StereoCenters(struct GRAPH Graph);
int FourNeigh(struct GRAPH Graph, int i, struct ORDER *Order);
int TwoMore(struct GRAPH Graph, int i,int valence);
int SameDir(struct GRAPH Graph, int i, int head, int end);		
struct COORD DiffVectors(struct COORD vecA, struct COORD vecB);
struct COORD TestPlanarVectors(struct COORD a, struct COORD b, struct COORD c);
int AdjustCoordinates(struct GRAPH Graph, struct COORD *vec1, int i, int j);
int BitTestPlanarVectors(struct COORD a, struct COORD b, struct COORD c);
struct COORD PerToFirst_OppToSecond(struct COORD vec1, struct COORD vec2);	
int SwitchTestPlanarVectors(struct COORD a, struct COORD b, struct COORD c, struct COORD d);
struct COORD ExternProduct(struct COORD vecA, struct COORD vecB);
double InternalProduct(struct COORD vecA, struct COORD vecB);
int	StereoCenterIsa(enum atomic_number NameA, int valence,struct GRAPH Graph,int i, struct ORDER* Order);

#endif				

			




