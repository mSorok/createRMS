/*
 *  Const.c
 *  
 *
 *  Created by Davide Fichera on 06/01/10.
 *
 */
/*
 In this file all global variables declared in DFgeneral.h are given
 and some more routines too.
 */ 

#include "DFgeneral.h"

VerboseDebug = 0;
Verbose = 0;
atomnameP = 1;


const struct ORDER EMPTY_O = { 0,0,0,0};
const struct ORDER LEVO = { 1,2,3,4};
const struct ORDER DEXTRO = { 1,2,4,3};
const struct COORD ZeroVector = {0,0,0 }; 

int Len_Table = 0;
char** Table_name = NULL;
int* Table_int = NULL;

enum atomic_number (*Name2Number)(char* atom) = Name2Number_atom;

enum atomic_number Name2Number_table(char* atom){
	int counter = 0;
	while(counter < Len_Table){
		if(strcmp(atom,Table_name[counter]) == OK) return Table_int[counter];
		counter++;
	}
	return 0;
}

enum atomic_number Name2Number_atom(char* atom){
	char		element[MaxNom];
	enum atomic_number number;
	
	number = 0;
	
	if (atom == NIL) return(0.0);
	strcpy(element,atom);	
	if (strcmp(element,"H") == OK) number = H; 
	if (strcmp(element,"He") == OK) number = He;
	if (strcmp(element,"Li") == OK) number = Li;
	if (strcmp(element,"Be") == OK) number = Be;
	if (strcmp(element,"B") == OK) number = B;
	if (strcmp(element,"C") == OK) number = C;
	if (strcmp(element,"N") == OK) number = N; 
	if (strcmp(element,"O") == OK) number = O;
	if (strcmp(element,"F") == OK) number = F;
	if (strcmp(element,"Ne") == OK) number =Ne;
	if (strcmp(element,"Na") == OK) number = Na;
	if (strcmp(element,"Mg") == OK) number = Mg;
	if (strcmp(element,"Al") == OK) number = Al;
	if (strcmp(element,"Si") == OK) number = Si;
	if (strcmp(element,"P") == OK) number = P;
	if (strcmp(element,"S") == OK) number = S;
	if (strcmp(element,"Cl") == OK) number = Cl;
	if (strcmp(element,"Ar") == OK) number = Ar;
	if (strcmp(element,"K") == OK) number = K;
	if (strcmp(element,"Ca") == OK) number = Ca;
	if (strcmp(element,"Sc") == OK) number = Sc;
	if (strcmp(element,"Ti") == OK) number = Ti;
	if (strcmp(element,"V") == OK) number = V;
	if (strcmp(element,"Cr") == OK) number = Cr;
	if (strcmp(element,"Mn") == OK) number = Mn;
	if (strcmp(element,"Fe") == OK) number = Fe;
	if (strcmp(element,"Ni") == OK) number = Ni;
	if (strcmp(element,"Co") == OK) number = Co;
	if (strcmp(element,"Cu") == OK) number = Cu;
	if (strcmp(element,"Zn") == OK) number = Zn;
	if (strcmp(element,"Ga") == OK) number = Ga;
	if (strcmp(element,"Ge") == OK) number = Ge;
	if (strcmp(element,"As") == OK) number = As;
	if (strcmp(element,"Se") == OK) number = Se;
	if (strcmp(element,"Br") == OK) number = Br;
	if (strcmp(element,"Kr") == OK) number = Kr;
	if (strcmp(element,"Rb") == OK) number = Rb;
	if (strcmp(element,"Sr") == OK) number = Sr;
	if (strcmp(element,"Y") == OK) 	number = Y;
	if (strcmp(element,"Zr") == OK) number = Zr;
	if (strcmp(element,"Nb") == OK) number = Nb;
	if (strcmp(element,"Mo") == OK) number = Mo;
	if (strcmp(element,"Tc") == OK) number = Tc;
	if (strcmp(element,"Ru") == OK) number = Ru;
	if (strcmp(element,"Rh") == OK) number = Rh;
	if (strcmp(element,"Pd") == OK) number = Pd;
	if (strcmp(element,"Ag") == OK) number = Ag;
	if (strcmp(element,"Cd") == OK) number = Cd;
	if (strcmp(element,"In") == OK) number = In;
	if (strcmp(element,"Sn") == OK) number = Sn;
	if (strcmp(element,"Sb") == OK) number = Sb;
	if (strcmp(element,"Te") == OK) number = Te;
	if (strcmp(element,"I") == OK) number = I;
	if (strcmp(element,"Xe") == OK) number = Xe;
	if (strcmp(element,"Cs") == OK) number = Cs;
	if (strcmp(element,"Ba") == OK) number = Ba;
	if (strcmp(element,"La") == OK) number = La;
	if (strcmp(element,"Ce") == OK) number = Ce;
	if (strcmp(element,"Pr") == OK) number = Pr;
	if (strcmp(element,"Nd") == OK) number = Nd;
	if (strcmp(element,"Pm") == OK) number = Pm;
	if (strcmp(element,"Sm") == OK) number = Sm;
	if (strcmp(element,"Eu") == OK) number = Eu;
	if (strcmp(element,"Gd") == OK) number = Gd;
	if (strcmp(element,"Tb") == OK) number = Tb;
	if (strcmp(element,"Dy") == OK) number = Dy;
	if (strcmp(element,"Ho") == OK) number = Ho;
	if (strcmp(element,"Er") == OK) number = Er;
	if (strcmp(element,"Tm") == OK) number = Tm;
	if (strcmp(element,"Yb") == OK) number = Yb;
	if (strcmp(element,"Lu") == OK) number = Lu;
	if (strcmp(element,"Hf") == OK) number = Hf;
	if (strcmp(element,"Ta") == OK) number = Ta;
	if (strcmp(element,"W") == OK) number = W;
	if (strcmp(element,"Re") == OK) number = Re;
	if (strcmp(element,"Os") == OK) number = Os;
	if (strcmp(element,"Ir") == OK) number = Ir;
	if (strcmp(element,"Pt") == OK) number = Pt;
	if (strcmp(element,"Au") == OK) number = Au;
	if (strcmp(element,"Hg") == OK) number = Hg;
	if (strcmp(element,"Tl") == OK) number = Tl;
	if (strcmp(element,"Pb") == OK) number = Pb;
	if (strcmp(element,"Bi") == OK) number = Bi;
	if (strcmp(element,"Po") == OK) number = Po;
	if (strcmp(element,"At") == OK) number = At;
	if (strcmp(element,"Rn") == OK) number = Rn;
	if (strcmp(element,"Fr") == OK) number = Fr;
	if (strcmp(element,"Ra") == OK) number = Ra;
	if (strcmp(element,"Ac") == OK) number = Ac;
	if (strcmp(element,"Th") == OK) number = Th;
	if (strcmp(element,"Pa") == OK) number = Pa;
	if (strcmp(element,"U") == OK) number = U;
	if (strcmp(element,"Np") == OK) number = Np;
	if (strcmp(element,"Pu") == OK) number = Pu;
	if (strcmp(element,"Am") == OK) number = Am;
	if (strcmp(element,"Cm") == OK) number = Cm;
	if (strcmp(element,"Bk") == OK) number = Bk;
	if (strcmp(element,"Cf") == OK) number = Cf;
	if (strcmp(element,"Es") == OK) number = Es;
	if (strcmp(element,"Fm") == OK) number = Fm;
	if (strcmp(element,"Md") == OK) number = Md;
	if (strcmp(element,"No") == OK) number = No;
	if (strcmp(element,"Lr") == OK) number = Lr;
	if (strcmp(element,"Rf") == OK) number = Rf;
	if (strcmp(element,"Db") == OK) number = Db;
	if (strcmp(element,"Sg") == OK) number = Sg;
	if (strcmp(element,"Bh") == OK) number = Bh;
	if (strcmp(element,"Hs") == OK) number = Hs;
	if (strcmp(element,"Mt") == OK) number = Mt;
	if (strcmp(element,"Ds") == OK) number = Ds;
	if (strcmp(element,"Rg") == OK) number = Rg;
	if (strcmp(element,"Uub") == OK) number = Uub;
	if (strcmp(element,"Uut") == OK) number = Uut;
	if (strcmp(element,"Uuq") == OK) number = Uuq;
	if (strcmp(element,"Uup") == OK) number = Uup;
	if (strcmp(element,"Uuh") == OK) number = Uuh;
	if (strcmp(element,"Uuo") == OK) number = Uuo;
	if(strcmp(element,"X") == OK)number = X;
	if(strcmp(element,"R") == OK)number = R;
	if(strcmp(element,"N+") == OK)number = Nplus;
	if (strcmp(element,"C@") == OK) number = Cat;
	if (strcmp(element,"C@@") == OK) number = Catat;
	if (strcmp(element,"N@") == OK) number = Nat;
	if (strcmp(element,"N@@") == OK) number = Natat;
	if (strcmp(element,"P@") == OK) number = Pat;
	if (strcmp(element,"P@@") == OK) number = Patat;
	if (strcmp(element,"S@") == OK) number = Sat;
	if (strcmp(element,"S@@") == OK) number = Satat;
	if (strcmp(element,"N+@") == OK) number = Nplusat;
	if (strcmp(element,"N+@@") == OK) number = Nplusatat;

	return(number);	
}


