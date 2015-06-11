/*
 *  FindGeometry.c
 *  
 *
 *  Created by Davide Fichera on 01/10.
 *  Updated by Pablo Carbonell on 03/12
 *
 */
/*
 Functions needed to determine the relative positions of atoms. 
 */

#include "DFgeneral.h"
#include "DFFindGeometry.h"

struct COORD ExternProduct(struct COORD vecA, struct COORD vecB){
	struct COORD veco;
	veco.x = vecA.y * vecB.z - vecB.y * vecA.z;
 	veco.y = vecA.z * vecB.x - vecB.z * vecA.x;
	veco.z = vecA.x * vecB.y - vecB.x * vecA.y;
	return veco;
}

int StereoCenterIsa(enum atomic_number NameA, int valence,struct GRAPH Graph,int i, struct ORDER* Order){
	int j, NumVal, NumbH;
	if((Graph.atom_identity[i].Name == NameA) && (Graph.NNeigh[i] > 2)){//An atom linked with at most one hydrogen
		for(j=0;j<Graph.NNeigh[i];j++){
			if (Graph.bond_identity[Graph.Edge[i][j]].order == 4) {
				Graph.atom_identity[i].candidate = 0;
				return 0; // PC: atom in an aromatic ring
			}
		}
		NumVal = 0;//Count the number of links.
		for(j=0;j<Graph.NNeigh[i];j++){
			NumVal+= Graph.bond_identity[Graph.Edge[i][j]].order_mod ;
		}
		NumbH = 0;//Count the number of links.
		for(j=0;j<Graph.NNeigh[i];j++){
			if(Graph.atom_identity[Graph.Neigh[i][j]].Name != H)
				NumbH += Graph.bond_identity[Graph.Edge[i][j]].order_mod ;
		}
		NumbH = valence - NumbH;
		//valence - NumVal = NumbOfHydrogen
		if (Graph.NNeigh[i] + valence - NumVal == 4){//Real coordination is 4 is the same as two lines below
			NumVal +=  4 - valence;//unchanged for Carbons, -1 for Phosphorus
			//valence 
			if(NumVal == Graph.NNeigh[i]){//Then the carbon is connected with four atoms and at most one is an hydrogen. 
				//As two lines above with NumVal = either 3 or 4
				Order[0] = EMPTY_O;
				Graph.atom_identity[i].candidate = 1;
				if(Verbose)
					printf("Candidate Carbon %d NV: %d \n",i,NumVal);
				if(NumVal == 3){//There is exactly one hydrogen
					ThreePlusH(Graph, i, Order);
				}
				if(NumVal == 4){//There is no hydrogen, we have coordinates of four atoms. 
					FourNeigh(Graph, i, Order);
				}
				if(Graph.atom_identity[i].ordered == 0)  // An atom order of 5 means unknown and not flagged in the file
					Graph.atom_identity[i].order = 5 +  Graph.atom_identity[i].stereo;
				if(VerboseDebug)
					printf("Atom %d Order: %d (>4 means unknown)\n",i+1,Graph.atom_identity[i].order);
					// printf("The order we found is : %d if it is larger then 4 we did not find it and we just report the stereoch + 5\n",Graph.atom_identity[i].order);
			}
		}
	}
	return 0;
}


//We look atom by atom for a carbon with four neighbours, since we are not interested in indistinguishable 
//branches we discard Carbons with more than one hydrogen as neighbours.
//Then, if the atom is a candidate for being a stereocenter we look at its geometry by calling 
//a test function: TreePlusH (if there is one hydrogen) or FourNeigh if there are no hydrogens. 
//If we did not find any information about geometry we copy the one written into the molfile as 
//atom stereo-parity and shift it by 5 so that 5=not stereo, 6 odd, 7 even, 8 either or unmarked 
int dffg_StereoCenters(struct GRAPH Graph){
	int size = Graph.size;
	int i;
	struct ORDER Order;
	
	for(i=0;i<size;i++){
		StereoCenterIsa(C,4,Graph,i,&Order);
		StereoCenterIsa(Si,4,Graph,i,&Order);
		StereoCenterIsa(Nplus,4,Graph,i,&Order);
		StereoCenterIsa(P,5,Graph,i,&Order);
		StereoCenterIsa(S,5,Graph,i,&Order);//Trivalent S, invents a phantom hydrogen
		StereoCenterIsa(S,6,Graph,i,&Order);
	
	}
	return 0;
}

//This function set the bond_identity[].align of all C=C edges such that the two C are connected to 
//two more atoms, each with at most one hydrog as neighbour. 
int dffg_StereoBonds(struct GRAPH Graph){
	int i;
	int size = Graph.NumbEdges;
	int valence1, valence2;

	for(i=0;i<size;i++)
		if(Graph.bond_identity[i].order == 2){
			int head , end;
			head = Graph.head[i];
			end = Graph.end[i];
			if((Graph.NNeigh[head] > 1) && (Graph.NNeigh[end] > 1) && (Graph.NNeigh[head] < 4) && (Graph.NNeigh[end] < 4)){
			if(((Graph.atom_identity[head].Name == C ) || (Graph.atom_identity[head].Name == N) || (Graph.atom_identity[head].Name == Nplus)) &&
				   ((Graph.atom_identity[end].Name == C) || (Graph.atom_identity[end].Name == N)  || (Graph.atom_identity[head].Name == Nplus) )){

					valence1 = 4;
					valence2 = 4;
					if(Graph.atom_identity[head].Name == N) 
						valence1 = 3;
					if(Graph.atom_identity[end].Name == N) 
						valence2 = 3;
					if(Graph.atom_identity[head].charge == 3)
						valence1 ++;
					if(Graph.atom_identity[end].charge == 3)
						valence2 ++;
				
					if(TwoMore(Graph, head, valence1) && TwoMore(Graph, end, valence2)){ //The double bond i is linked to 
						//two carbons or nitrogen 
						//(head and end) each connected 
						//to two(if C or 1 if N) more atoms and at most one of its neighb is a hydrog. 
						if(Verbose)
							printf("Candidate bond %d extremes %d %d ",i,head,end);
						Graph.bond_identity[i].align = SameDir(Graph,i,head,end);
						if(Verbose)
							printf("it is of kind %d \n",Graph.bond_identity[i].align);
					}
				}
			}
			

		}
	return 0;	
}

//Make vectors ingoing onto the stereocenter and then check if they are planar if they aren't the position of atoms is 
//corrected with informations about bonds (up or down), planarity is testes once more.
//if we still know nothing about the position of H we set.ordered = 0 and leave the order structure unchanged.
//Otherwise the order is given according to the sign of (a^b)*c 
int ThreePlusH(struct GRAPH Graph, int i, struct ORDER *Order){
	struct COORD *vec1;
	vec1 = calloc(3,sizeof(struct COORD));
	vec1[0] = DiffVectors(Graph.atom_identity[Graph.Neigh[i][0]].Coord,Graph.atom_identity[i].Coord);
	vec1[1] = DiffVectors(Graph.atom_identity[Graph.Neigh[i][1]].Coord,Graph.atom_identity[i].Coord);
	vec1[2] = DiffVectors(Graph.atom_identity[Graph.Neigh[i][2]].Coord,Graph.atom_identity[i].Coord);
	struct COORD vec4 = TestPlanarVectors(vec1[0],vec1[1],vec1[2]);
	if( fabs(vec4.x) + fabs(vec4.y) + fabs(vec4.z) == 0 ){
		//Check if one of them is single bound up. 
		int j;
		for(j=0;j<3;j++)		
			AdjustCoordinates(Graph, vec1 + j, i, j);
	}
	vec4 = TestPlanarVectors(vec1[0],vec1[1],vec1[2]);
	if(  fabs(vec4.x) + fabs(vec4.y) + fabs(vec4.z) == 0 ){
		if(Verbose)
			printf("We do not know where the H is\n");
		Graph.atom_identity[i].ordered = 0;
	}
	else{
		int sign;	
		sign =	BitTestPlanarVectors(vec1[0] ,vec1[1] ,vec1[2]);
		if(sign == 1){
			Order[0] = DEXTRO;
			Graph.atom_identity[i].order = 1;	
		}
		if(sign == -1){
			Order[0] = LEVO;
			Graph.atom_identity[i].order = -1;	
		}
		if(sign != 0)
			Graph.atom_identity[i].ordered = 1;
		if(Verbose)
			printf("Geometrical parity %d\n",sign);
	}
	free(vec1);
	return 0;
}


//Shift up or down the z coordinate of an amount between 0.5 and 0.7 the length of the vector.
int AdjustCoordinates(struct GRAPH Graph, struct COORD *vec1, int i, int j){
	double Factor = 1.0;
	if(Graph.end[Graph.Edge[i][j]] == i)
		Factor = 1.0;
	else
		Factor = -1.0;
	if(Graph.bond_identity[Graph.Edge[i][j]].order == 1) 
		switch(Graph.bond_identity[Graph.Edge[i][j]].stereo){
			case 1 :
				vec1[0].z += 0.5 * (fabs(vec1[0].x) + fabs(vec1[0].y)) * Factor;
				break;
			case 6:
				vec1[0].z -= 0.5 * (fabs(vec1[0].x) + fabs(vec1[0].y)) * Factor;
				break;
			default: ;
				break;
		}
	return Graph.bond_identity[Graph.Edge[i][j]].stereo;
}

//Given the coordinates-vectors from stereocenters looks for a non planar triple. If there aren't includes 
//informations about bond directions up/down the plane and look for a non planar triple.  
//The tag 1,...,4 of the link not included in the non planar triple is returned toghether with a sign stating 
//if the geometry il D or L.
//If there are no planar triple the carbon is listed as not ordered. 
int FourNeigh(struct GRAPH Graph, int i, struct ORDER *Order){
	struct COORD vec1, vec2, vec3, vec4;
	vec1 = DiffVectors(Graph.atom_identity[Graph.Neigh[i][0]].Coord, Graph.atom_identity[i].Coord);
	vec2 = DiffVectors(Graph.atom_identity[Graph.Neigh[i][1]].Coord, Graph.atom_identity[i].Coord);
	vec3 = DiffVectors(Graph.atom_identity[Graph.Neigh[i][2]].Coord, Graph.atom_identity[i].Coord);
	vec4 = DiffVectors(Graph.atom_identity[Graph.Neigh[i][3]].Coord, Graph.atom_identity[i].Coord);
	

	int which = SwitchTestPlanarVectors(vec1,vec2,vec3,vec4);	
	if(which == 0){
		AdjustCoordinates(Graph, &vec1, i, 0);	
		AdjustCoordinates(Graph, &vec2, i, 1);
		AdjustCoordinates(Graph, &vec3, i, 2);
		AdjustCoordinates(Graph, &vec4, i, 3);
		which = SwitchTestPlanarVectors(vec1,vec2,vec3,vec4);
	}
	if(which > 0){
		Order[0] = DEXTRO;
		Graph.atom_identity[i].order = 1;
	}
	if(which < 0){
		Order[0] = LEVO;
		Graph.atom_identity[i].order = -1;	
	}
	Graph.atom_identity[i].order = which; // PC March 2013
	if(which == 0)
		Graph.atom_identity[i].ordered = 0;
	else
		Graph.atom_identity[i].ordered = 1;
	
	if(Verbose)
		printf("Geometrical parity %d\n",which);
	
	return 0;
}

//If there is exactly one double edge and there is at most one hydrogen then returns 1.
int TwoMore(struct GRAPH Graph, int i, int valence){
	int j;
	int NBonds = 0;
	int NNeigh = Graph.NNeigh[i];
	for(j=0;j<NNeigh;j++)
		if(Graph.atom_identity[Graph.Neigh[i][j]].Name != H)
			NBonds += Graph.bond_identity[Graph.Edge[i][j]].order_mod;
		else
			NNeigh -=1;  // PC 05/12 Fixed when working with explicit H
	if((NBonds - NNeigh) == 1){
		NBonds += 4 - valence;//The same for Carbon
		if((NBonds == 4) ||(NBonds == 3))
			return 1;
	}
	return 0;
}

//We look at the double link and the first (different) edge of both ends. In this way we avoid 
//complications coming by the presence of hydrogens. 
//We then construct two auxiliar vectors and compare their directions. 
//The function returns their alignment.
int SameDir(struct GRAPH Graph, int i, int head, int end){
	struct COORD vec1, vec2, vec3, vec4, AuxV1, AuxV2, auxvec_h, auxvec_e; 
	int first = 0;
	if(Graph.Neigh[head][first] == end)
		first++;
	int second = 0;
	if(Graph.Neigh[end][second] == head)
		second++;
	
	if(Graph.NNeigh[head] == 3){
		int auxfirst = 0;
		while((Graph.Neigh[head][auxfirst] == end) || auxfirst == first)
			auxfirst++;
		auxvec_h = DiffVectors(Graph.atom_identity[Graph.Neigh[head][auxfirst]].Coord,Graph.atom_identity[head].Coord);
	}
	
	if(Graph.NNeigh[end] == 3){
		int auxsecond = 0;
		while((Graph.Neigh[end][auxsecond] == head) || auxsecond == second)
			auxsecond++;
		auxvec_e = DiffVectors(Graph.atom_identity[Graph.Neigh[end][auxsecond]].Coord,Graph.atom_identity[end].Coord);
	}
	
	vec1 = DiffVectors(Graph.atom_identity[head].Coord,Graph.atom_identity[end].Coord);
	vec2 = DiffVectors(Graph.atom_identity[Graph.Neigh[end][second]].Coord,Graph.atom_identity[end].Coord);
	vec3 = DiffVectors(Graph.atom_identity[end].Coord,Graph.atom_identity[head].Coord);
	vec4 = DiffVectors(Graph.atom_identity[Graph.Neigh[head][first]].Coord,Graph.atom_identity[head].Coord);

	if(Graph.NNeigh[head] == 3)
		vec4 = DiffVectors(vec4,auxvec_h);
	if(Graph.NNeigh[end] == 3)
		vec2 = DiffVectors(vec2,auxvec_e);

		
	AuxV1 = PerToFirst_OppToSecond(vec1, vec2) ;
	AuxV2 = PerToFirst_OppToSecond(vec3, vec4) ;
		
	double Prod = InternalProduct(AuxV1, AuxV2);
	if(fabs(Prod) < 0.0001){
		AdjustCoordinates(Graph,&vec1,end,1 - second);
		AdjustCoordinates(Graph,&vec2,end,second);
		AdjustCoordinates(Graph,&vec3,head,1 - first);
		AdjustCoordinates(Graph,&vec4,end,first);
		
		AuxV1 = PerToFirst_OppToSecond(vec1, vec2) ;
		AuxV2 = PerToFirst_OppToSecond(vec3, vec4) ;
		Prod = InternalProduct(AuxV1, AuxV2);
	}
	else{
		if(Prod > 0)
			return 1;
		if(Prod < 0)
			return -1;
	}
	return 0;
}

struct COORD DiffVectors(struct COORD vecA, struct COORD vecB){
	struct COORD vec1 = ZeroVector;
	vec1.x = vecA.x - vecB.x;
	vec1.y = vecA.y - vecB.y;
	vec1.z = vecA.z - vecB.z;
	return vec1;
}

//If the vectors are planar returns zerovector else returns -Average
struct COORD TestPlanarVectors(struct COORD a, struct COORD b, struct COORD c){
	struct COORD Aux;
	Aux = ExternProduct(a,b);
	if(InternalProduct(Aux,c)){
		Aux.x = -(a.x + b.x + c.x)/ 3.0;
		Aux.y = -(a.y + b.y + c.y)/ 3.0;
		Aux.z = -(a.z + b.z + c.z)/ 3.0;
	}
	else
		return ZeroVector;
	return Aux;
}

//Returns the sign of (a^b)*c eventually zero.
int BitTestPlanarVectors(struct COORD a, struct COORD b, struct COORD c){
	struct COORD Aux;
	double Prod;
	Aux = ExternProduct(a,b);
	Prod = InternalProduct(Aux,c);
	if(Prod > 0)
		return 1;
	if(Prod < 0)
		return -1;
	return 0;
}

struct COORD PerToFirst_OppToSecond(struct COORD vec1, struct COORD vec2){
	struct COORD vec3, vec4;
	vec3 = ExternProduct(vec1,vec2);
	vec4 = ExternProduct(vec1,vec3);
	return vec4;
}


//We make several hypothesys: all that works if we have not pathology on the geometry 
//i.e. if the four atoms do not belong to the same half-space. 
//The four possible triples are analyzed one by one and the first non planar one is returned with a sign 
//about their geometry.
//This sign depends on the triple so that we can infer the geometry just by looking at the sign we do not need 
//to look at the specific triple. 
int SwitchTestPlanarVectors(struct COORD a, struct COORD b, struct COORD c, struct COORD d){
	int Aux;
	int Aux1,Aux2,Aux3,Aux4;

	// Compute signed volume of the tetrahedron
	Aux1 = BitTestPlanarVectors(b,c,d);
	Aux2 = BitTestPlanarVectors(c,d,a);
	Aux3 = BitTestPlanarVectors(d,a,b);
	Aux4 = BitTestPlanarVectors(a,b,c);



		// Check if 3 of the ligands are coplanar
	if (Aux1*Aux2*Aux3*Aux4 !=0) { // PC test March 2013: Encode all combinations
		Aux = 0;
		// Reorient the z plane
		if (a.z>0) Aux1 *= -1;
		if (b.z>0) Aux2 *= -1;
		if (c.z>0) Aux3 *= -1;
		if (d.z>0) Aux4 *= -1;
		if (Aux1==1)	Aux += 2000;
		if (Aux2==1)	Aux += 200;
		if (Aux3==1) 	Aux += 20;
		if (Aux4==1) 	Aux += 2;
		if (Aux1==-1)	Aux += 1000;
		if (Aux2==-1)	Aux += 100;
		if (Aux3==-1) 	Aux += 10;
		if (Aux4==-1) 	Aux += 1;

		return(Aux);
	}

//	Aux = BitTestPlanarVectors(a,b,c);
	if(Aux4 == 1)
		return 4;
	if(Aux4 == -1)
		return -4;
	
//	Aux = BitTestPlanarVectors(a,b,d);
	if(Aux3 == 1)
		return 3;
	if(Aux3 == -1)
		return -3;
	
//	Aux = BitTestPlanarVectors(a,c,d);
	if(Aux2 == 1)
		return 2;
	if(Aux2 == -1)
		return -2;
	
//	Aux = BitTestPlanarVectors(b,c,d);
	if(Aux1 == 1)
		return 1;
	if(Aux1 == -1)
		return -1;
	
	return 0;
}	


double InternalProduct(struct COORD vecA, struct COORD vecB){
	return vecA.x * vecB.x + vecA.y * vecB.y + vecA.z * vecB.z;  
}

