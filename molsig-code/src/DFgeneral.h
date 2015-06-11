/*
 *  
 *  
 *
 *  Created by Davide Fichera on 11/12/09.
 *
 */
#ifndef	DFGENERAL_H
#define DFGENERAL_H


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "c.h"
#include "general.h"

# define MaxNom                256
# define MAXCHAINE            256
# define MaxChaine            256
#define MAXCHARLINE     1024 /* MAX NUMBER OF CHAR PER LINE */
int VerboseDebug;// = 0;
int Verbose ;//= 0;
int atomnameP ;//= 1;
int Len_Table;
char** Table_name;
int* Table_int;

# define ERROR                -1L     
# define NIL                  0L
# define FALSE         0
# define OK                   NIL    /* can be also a pointer     */
# define TRUE                1
#define  D3      3
#define EPS         1.e-3  
# define TAS_WORD_ALIGN(T)                   ((T) &= 0xfffffffe)
/* alignement sur le mot de 16 bit */
#define	NORM(v)	(sqrt( (v)[0] * (v)[0] + (v)[1] * (v)[1] + (v)[2] * (v)[2] ))

struct Vertices{
	int processed;
	int NChild;
	int *Child;
	int *bond_order;
	int When;
	int Father; //it is the current father, we do not need to know more.
	int ifrepeat;
	int whenrepeat;
	int numbOfAt;
	int imported_invariant;
};
enum signfile {
	Bra
	,Ket
	,Atom
	,Bond
};
struct ReadVertexOut{
	enum signfile type;
	int atomtag;
	int bond_order;
	int invariant;
};
struct ReadVertexOut ReadChars();


struct ORDER {
	int num1, num2, num3, num4;	
};
extern const struct ORDER EMPTY_O;
extern const struct ORDER LEVO;
extern const struct ORDER DEXTRO;

//struct  COORD {
  //  double  x,y,z;
//};// t_coord, *p_coord;
const struct COORD ZeroVector; 

//enum atomic_number {
//	H = 1, He,
//	Li,  Be,  B,   C = 1060,   N = 1070,  O = 1080,  F,   Ne,
//	Na,  Mg,  Al,  Si,  P,  S,  Cl,  Ar,
//	K,  Ca,  Sc,  Ti,  V,  Cr,  Mn,  Fe,  Co, Ni, Cu,  Zn,  Ga,  Ge,  As,  Se,  Br,  Kr,
//	Rb,  Sr,  Y,  Zr,  Nb,  Mo,  Tc,  Ru,  Rh,  Pd,  Ag,  Cd,  In,  Sn,  Sb,  Te,  I,  Xe,
//	Cs,  Ba,
//	La,  Ce,  Pr,  Nd,  Pm,  Sm,  Eu,  Gd,  Tb,  Dy,  Ho,  Er,  Tm,  Yb,  Lu,
//	Hf,  Ta,  W, Re,  Os,  Ir,  Pt,  Au,  Hg,  Tl,  Pb,  Bi,  Po,  At,  Rn,
//	Fr,  Ra,
//	Ac,  Th,  Pa,  U, Np,  Pu,  Am,  Cm,  Bk,  Cf,  Es,  Fm,  Md,  No,  Lr,
//	Rf,  Db,  Sg,  Bh, Hs,  Mt,  Ds, Rg,  Uub,  Uut,  Uuq = 114,  Uup,  Uuh,  Uuo,
//	X, R, Nplus , Cat = 1059, Catat = 1061
//}
enum atomic_number {
	H = 1010, He = 1020,
	Li = 1030,  Be = 1040,  B=1050,   C = 1060,   N = 1070,  O = 1080,  F = 1090,   Ne = 1100,
	Na = 1110,  Mg = 1120,  Al = 1130 ,  Si = 1140 ,  P = 1150 ,  S = 1160,  Cl = 1170,  Ar = 1180,
	K = 1190, Ca = 1200,  Sc = 1210,  Ti = 1220,  V = 1230,  Cr = 1240,  Mn = 1250,  Fe = 1260,  Co = 1270, Ni = 1280, Cu = 1290,  Zn = 1300 ,  Ga = 1310,  Ge = 1320,  As = 1330,  Se = 1340,  Br = 1350,  Kr = 1360,
	Rb = 1370,  Sr = 1380,  Y = 1390,  Zr = 1400,  Nb = 1410,  Mo = 1420,  Tc = 1430,  Ru = 1440,  Rh = 1450,  Pd = 1460,  Ag = 1470,  Cd = 1480,  In = 1490,  Sn = 1500,  Sb = 1510,  Te = 1520,  I = 1530,  Xe = 1540,
	Cs = 1550,  Ba = 1560,
	La = 1570,  Ce = 1580,  Pr = 1590,  Nd = 1600,  Pm = 1610,  Sm = 1620,  Eu = 1630,  Gd = 1640,  Tb = 1650,  Dy = 1660,  Ho = 1670,  Er = 1680,  Tm = 1690,  Yb = 1700,  Lu = 1710,
	Hf = 1720,  Ta = 1730,  W = 1740, Re = 1750,  Os = 1760,  Ir = 1780,  Pt = 1790,  Au = 1800,  Hg = 1810,  Tl = 1820,  Pb = 1830,  Bi = 1840,  Po = 1850,  At = 1860,  Rn = 1870,
	Fr = 1880,  Ra = 1890,
	Ac = 1900,  Th = 1910,  Pa = 1920,  U = 1930, Np = 1940,  Pu = 1950,  Am = 1960,  Cm = 1970,  Bk = 1980,  Cf = 1990,  Es = 2000,  Fm = 2010,  Md = 2020,  No = 2030,  Lr = 2040,
	Rf = 2050,  Db = 2060,  Sg = 2070,  Bh = 2080, Hs = 2090,  Mt = 2100,  Ds = 2110, Rg = 2120,  Uub = 2130,  Uut = 2140,  Uuq = 2150,  Uup = 2160,  Uuh = 2170,  Uuo = 2180,
	X = 0, R = 0, Nplus = 2210 , Cat = 1060, Catat = 1060, Pat = 1149, Patat = 1148, Sat = 1159, Satat = 1158, Nplusat = 2209, Nplusatat = 2208, Nat = 1069, Natat = 1068
};
//1059 1058
//Begin: Graph Structure (as a directed graph)

//mass_diff, charge and stereoatom as defined in .mol files; 
enum mass_diff{
m3 = -3, m2  ,m1  , m0, p1, p2, p3 	
};
enum charge{
uncharged = 0, plus3, plus2, plus1, doubleradical, minus1, minus2, minus3 = 7	
};
enum stereo_atom{
notstereo = 0, odd, even, either = 3	
};

struct TESTPLANAR {
	struct COORD Coord;
	int Test;
};

struct ATOM_IDENTITY{
	int tag;//id_number in the graph
	//INTRINSIC PROPERTIES
	enum atomic_number Name;
	char charName[12];
	char completeName[12];
	int NumbNeigh;
	struct COORD Coord;
	enum mass_diff mass;
	enum charge charge;
	enum stereo_atom stereo; 	
	
	//AFTER CANONIZATION
	int identity;
	int *idNeigh;
	
	//GEOMETRY
	int candidate;//Tell if it is a candidate as sereocenter
	int ordered;  //Tell if its neighbours have been geometrically ordered
	int order;
};

struct GRAPHMemory{
	int* idNeigh;
	struct ATOM_IDENTITY *atom_identity;
	struct BOND_IDENTITY *bond_identity;
	int *NNeigh;
	int *end;
	int *head;
	int **Neigh;
};

struct BOND_IDENTITY{
	int tag;
	int order;
	int order_mod;
	int stereo;
	//AFTER STEREOCH
	int candidate;
	int align;
};

typedef struct GRAPH {
	int size;
	int NumbEdges;
	int NumbDirectedBonds;
	int *NNeigh;
	int *end;
	int *head;
	int **Neigh;
	int **Edge;
	//AS A MOLECULE
	struct ATOM_IDENTITY *atom_identity;
	struct BOND_IDENTITY *bond_identity;
	//TO BE READ JUST BY CONSTRUCTORS:
	int NumbBond;//Useful to construct the Graph
	atom_p* MyAtom;  // PC Lookup table with the DAG
} d_graph_t , *d_graph_p;

//End: Graph Structure

//Usage: 
//struct GRAPH MyGraph;
//struct GRAPHMemory MyGraphMemory = MakeGRAPHMemory(10, 15, &MyGraph);
//AddBond(2,8,&MyGraph);
//MakeNeigh(&MyGraph);
//DeleteGraphMemory(MyGraph);

extern int MakeIdentity(struct ATOM_IDENTITY *ai, int i,double x, double y, double z, char *name, int mass, int charge, int stereo);

enum atomic_number Name2Number_atom(char* atom);
enum atomic_number Name2Number_table(char* atom);
enum atomic_number (*Name2Number)(char* atom);


extern enum atomic_number Name2Number_atom(char* atom);
int OnlyNumbers;
int PrintOriginalTag;

#endif


