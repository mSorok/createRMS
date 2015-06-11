/***********************************************************
 FILE VECTOR Vectorial functions
       Jean-loup Faulon - May 1991 Albuquerque N.M.
***********************************************************/
#include <general.h>
#include <eps.h>

#define	NORM(v)	(sqrt( (v)[0] * (v)[0] + (v)[1] * (v)[1] + (v)[2] * (v)[2] ))

/*------------- TRIVIAL FUNCTIONS -----------------------------------*/

DEFINE t_bool geve_null_vector(v)  
/**********************************************************************
  
**********************************************************************/   
t_real 	v[D3];
{ 
	if ( (EPS_EQ(v[0],ZERO)) && (EPS_EQ(v[1],ZERO)) && (EPS_EQ(v[2],ZERO)) )
	   return(TRUE);
	return(FALSE);
	}           

DEFINE t_err geve_normalise_vector(vect1,resultat)  
/**********************************************************************
  
**********************************************************************/   
t_real 	vect1[D3];
t_real 	resultat[D3];
{ 
 	t_real temp;    
 	short i;
 	temp = NORM(vect1);
 	if (fabs(temp) < EPS ) 
   	{temp = 0;  for (i=0;i<3; i++) resultat[i]= vect1[i];}
  	else {
   	if (fabs((temp - 1)) < EPS) temp  = 1;
   	for (i=0;i<3; i++) resultat[i]= vect1[i]/temp;
  	}
  	return(OK);
	}           

DEFINE t_err geve_equal_vector(coord1,coord2)
/**********************************************************************
  
**********************************************************************/  
t_real  coord1[D3],coord2[D3];
{              
  	t_real temp;
  	temp = fabs(coord1[0] - coord2[0]);
  	if (temp > EPS) return(ERREUR);
  	temp = fabs(coord1[1] - coord2[1]);
  	if (temp > EPS) return(ERREUR);
  	temp = fabs(coord1[2] - coord2[2]);
  	if (temp > EPS) return(ERREUR); 
  	return(OK);
	}


DEFINE t_err geve_compute_vector(coord1,coord2,result) 
/**********************************************************************
  compute vertor coord1->coord2
**********************************************************************/  
t_real  	coord1[D3],coord2[D3];
t_real      	result[D3];
{       
       	result[0] = coord2[0] - coord1[0];
       	result[1] = coord2[1] - coord1[1];
       	result[2] = coord2[2] - coord1[2];
       	geve_normalise_vector(result,result);
       	return(OK);
}

/*------------- VECTOR PRODUCTS -------------------------------------*/

DEFINE t_err geve_product_scalar(vect1,vect2,result)
/**********************************************************************
  
**********************************************************************/  
t_real 	vect1[D3],vect2[D3];
t_real  *result;
{ 
	if (result == NIL) return(ERROR);    
 	*result = vect1[0]*vect2[0]+vect1[1]*vect2[1]+vect1[2]*vect2[2];
 	if (fabs(*result) < EPS ) *result = 0; 
 	if (fabs((*result - 1)) < EPS) *result = 1;
 	if (fabs((*result + 1)) < EPS) *result = -1;
  	return(OK);
	}

DEFINE t_err geve_multiply_scalar(x,vect,result)
/**********************************************************************
  
**********************************************************************/  
t_real 	vect[D3],result[D3];
t_real	x;
{     
 	result[0] = x * vect[0];
 	if (fabs(result[0]) < EPS ) result[0] = 0; 
 	result[1] = x * vect[1];
 	if (fabs(result[1]) < EPS ) result[1] = 0; 
 	result[2] = x * vect[2];
 	if (fabs(result[2]) < EPS ) result[2] = 0; 
 	return(OK);
	}


DEFINE t_err geve_product_vector(vect1,vect2,resultat) 
/**********************************************************************
  
**********************************************************************/   
t_real 	vect1[D3],vect2[D3];
t_real 	resultat[D3];
{     
	t_real	r[D3];
	t_integer	i;
 	r[0] =  vect2[2]*vect1[1] - vect2[1]*vect1[2];
 	if (fabs(r[0]) < EPS ) r[0] = 0; 
 	if (fabs((r[0] - 1)) < EPS) r[0] = 1;
 	r[1] =  vect2[0]*vect1[2] - vect2[2]*vect1[0];
 	if (fabs(r[1]) < EPS ) r[1] = 0; 
 	if (fabs((r[1] - 1)) < EPS) r[1] = 1;
 	r[2] =  vect2[1]*vect1[0] - vect2[0]*vect1[1];
 	if (fabs(r[2]) < EPS ) r[2] = 0; 
 	if (fabs((r[2] - 1)) < EPS) r[2] = 1;
	for (i = 0; i < D3; i++) resultat[i] = r[i];
  	return(OK);
	}

DEFINE t_err geve_scalar_vector(vect1,vect2,resultat) 
/**********************************************************************
  
**********************************************************************/   
t_real 	vect1[D3],vect2[D3];
t_real 	resultat[D3];
{     
	t_real	r[D3];
	t_integer	i;
 	r[0] =  vect2[0]*vect1[0];
 	if (fabs(r[0]) < EPS ) r[0] = 0; 
 	if (fabs((r[0] - 1)) < EPS) r[0] = 1;
 	r[1] =  vect2[1]*vect1[1];
 	if (fabs(r[1]) < EPS ) r[1] = 0; 
 	if (fabs((r[1] - 1)) < EPS) r[1] = 1;
 	r[2] =  vect2[2]*vect1[2];
 	if (fabs(r[2]) < EPS ) r[2] = 0; 
 	if (fabs((r[2] - 1)) < EPS) r[2] = 1;
	for (i = 0; i < D3; i++) resultat[i] = r[i];
  	return(OK);
	}

DEFINE t_err geve_oppose_vector(vect,resultat) 
/**********************************************************************
  
**********************************************************************/   
t_real 	vect[D3],resultat[D3];
{     
 	resultat[0] = -vect[0];
 	resultat[1] = -vect[1];
 	resultat[2] = -vect[2];
 	return(OK);
	}

DEFINE t_err geve_substract_vector(vect1,vect2,resultat) 
/**********************************************************************
  r = v2 - v1
**********************************************************************/   
t_real 	vect1[D3],vect2[D3];
t_real 	resultat[D3];
{     
 	resultat[0] = vect2[0]-vect1[0];
 	if (fabs(resultat[0]) < EPS ) resultat[0] = 0; 
 	if (fabs((resultat[0] - 1)) < EPS) resultat[0] = 1;
 	resultat[1] = vect2[1]-vect1[1];
 	if (fabs(resultat[1]) < EPS ) resultat[1] = 0; 
 	if (fabs((resultat[1] - 1)) < EPS) resultat[1] = 1;
 	resultat[2] = vect2[2]-vect1[2];
 	if (fabs(resultat[2]) < EPS ) resultat[2] = 0; 
 	if (fabs((resultat[2] - 1)) < EPS) resultat[2] = 1;
 	return(OK);
	}

DEFINE t_err geve_add_vector(vect1,vect2,resultat)
/**********************************************************************
  
**********************************************************************/   
t_real 	vect1[D3],vect2[D3];
t_real 	resultat[D3];
{    
 	resultat[0] = vect2[0]+vect1[0];
 	if (fabs(resultat[0]) < EPS ) resultat[0] = 0; 
 	if (fabs((resultat[0] - 1)) < EPS) resultat[0] = 1;
 	resultat[1] = vect2[1]+vect1[1];
 	if (fabs(resultat[1]) < EPS ) resultat[1] = 0; 
 	if (fabs((resultat[1] - 1)) < EPS) resultat[1] = 1;
 	resultat[2] = vect2[2]+vect1[2];
 	if (fabs(resultat[2]) < EPS ) resultat[2] = 0; 
 	if (fabs((resultat[2] - 1)) < EPS) resultat[2] = 1;
  	return(OK);
	}

/*------------- ANGLE AND DISTANCE ----------------------------------*/

DEFINE t_err geve_distance_vector(vect1,vect2,resultat)  
/**********************************************************************
 WARNING the distance is the square euclidian distance. 
**********************************************************************/   
t_real 	vect1[D3],vect2[D3];
t_real 	*resultat;
{     
	t_real	vect3[D3];
	if(resultat == NIL) return(ERROR);
	geve_substract_vector(vect1,vect2,vect3);
	*resultat = vect3[0]*vect3[0]+vect3[1]*vect3[1]+vect3[2]*vect3[2];
	return(OK);
	}

DEFINE t_err	geve_angle_vector(v1,v2,n,a)
/**********************************************************************
 Compute angle of two vectors from their sinus and cosinus
 n is the normal to the plane containing the angle  
**********************************************************************/ 
t_real	 v1[D3],v2[D3];
p_coord	 n;
t_real   *a;
{
	t_real   n1,n2,s,c,v[D3],sign;

	if (a == NIL) return(ERROR);
	n1 = NORM(v1); if (EPS_EQ(n1,ZERO)) return(ERROR);
        n2 = NORM(v2); if (EPS_EQ(n2,ZERO)) return(ERROR);
        geve_product_scalar(v1,v2,&c); c = c / (n1 * n2);
        geve_product_vector(v1,v2,v);  s = NORM(v) / (n1 * n2);

	/* sinus sign */
	if (n) {
	   geve_product_scalar(v,n,&sign); 
	   sign = ((sign > 0) ? 1: -1);
	   s = sign * s;
	   }
       	if (EPS_EQ(c,ZERO)) c = ZERO;
       	if (EPS_EQ(s,ZERO)) s = ZERO;
       	*a = sqrt((c) * (c) + (s) * (s));
       	if (EPS_EQ(*a,ZERO)) return(ERREUR);
       	c = c / (*a); s = s / (*a);

       	if (EPS_EQ(c,ZERO)) {
          if (EPS_LT(s,ZERO)) *a = (t_real)3/(t_real)2 * PI;
          else                *a = PI/(t_real)2;
          return(OK);
       	  }   
       	*a = (t_real)atan((t_real)( s / c));
       	if (EPS_LT(c,ZERO)) {
          if (EPS_EQ(s,ZERO)) *a = PI;
          if (EPS_LT(s,ZERO)) *a = *a + PI;
          if (EPS_GT(s,ZERO)) *a = *a - PI;
       	  }
       	while (EPS_LT(*a,ZERO)) *a += 2 * PI;       
       	return(OK); 
	}

/*------------- MATRICES --------------------------------------------*/

DEFINE	t_err	geve_matrix_vector(v,m,result)
/**********************************************************************
 Compute r = v x m.
**********************************************************************/   
t_real	v[D3],m[D3][D3],result[D3];
{
	t_real		r[D3];
	t_integer	i,j;

	for (j = 0; j < D3; j++) { 
            r[j] = 0;
            for (i = 0; i < D3; i++) r[j] += v[i] * m[i][j];
            }
        for (i = 0; i < D3; i++) result[i] = r[i];
	return(OK);
	}

DEFINE	t_err	geve_oppose_matrix(m,resultat)
/**********************************************************************
 Oppose matrix m in resultat.
**********************************************************************/   
t_real	m[D3][D3],resultat[D3][D3];
{
	t_real		r[D3][D3];
	t_integer	i,j;
	for (i = 0; i < D3; i++) for (j = 0; j < D3; j++)
	    r[i][j] = m[j][i];
	for (i = 0; i < D3; i++) for (j = 0; j < D3; j++)
            resultat[i][j] = r[i][j];
	return(OK);
	}

DEFINE	t_err	geve_product_matrix(r1,r2,resultat)
/**********************************************************************
 Compute r = r1 x r2.
**********************************************************************/   
t_real	r1[D3][D3],r2[D3][D3];
t_real	resultat[D3][D3];
{
	t_integer	i,j,k;
	t_real		r[D3][D3];
	for (i = 0; i < D3; i++) for (j = 0; j < D3; j++) {
	    r[i][j] = 0;
	    for (k = 0; k < D3; k++)
                r[i][j] = r[i][j] + r1[i][k] * r2[k][j];
	    }
	for (i = 0; i < D3; i++) for (j = 0; j < D3; j++)
            resultat[i][j] = r[i][j];
	return(OK);
	}

/*------------- AXES FUNCTIONS --------------------------------------*/

LOCAL	t_coord	axe_x = {1,0,0};
LOCAL	t_coord	axe_y = {0,1,0};
LOCAL	t_coord	axe_z = {0,0,1};
LOCAL	t_coord	axe_0 = {0,0,0};

DEFINE	t_err	geve_normal_axe(axe)
/**********************************************************************
 Return the normelized system axe
***********************************************************************/
t_coord		axe[D3];
{
	if (axe == NIL) return(ERROR);
	axe[0] = axe_x; axe[1] = axe_y; axe[2] = axe_z;
	return(OK);
	}

DEFINE	t_err	geve_compute_axe(d,xyz,axe)
/**********************************************************************
 Compute a system of orthogonal axes form one vector.
***********************************************************************/
t_coord		d;
t_name		xyz;
t_coord		axe[D3];
{
	t_coord X,Y,Z;

	/* compute axe */
	X = d;
	geve_product_vector(&axe_x,&X,&Y);
        if (geve_equal_vector(&Y,&axe_0) == OK)
	geve_product_vector(&axe_y,&X,&Y);
        if (geve_equal_vector(&Y,&axe_0) == OK) return(ERROR);
	geve_normalise_vector(&Y,&Y);
	geve_product_vector(&X,&Y,&Z);
	if (strcmp(xyz,"x") == OK) {
	   axe[0] = X; axe[1] = Y; axe[2] = Z; return(OK); }
	else if (strcmp(xyz,"y") == OK) {
	   axe[0] = Z; axe[1] = X; axe[2] = Y; return(OK); }
	else if (strcmp(xyz,"z") == OK) {
	   axe[0] = Y; axe[1] = Z; axe[2] = X; return(OK); }
	else return(ERROR);
	}

DEFINE	t_err	geve_compute_rotation(axe1,axe2,ALPHA,rotation)
/**********************************************************************
 Compute rotation = (axe1 x axe2)t * r(ALPHA)
 note: ALPHA = ZERO if not needed
***********************************************************************/
t_coord		axe1[D3],axe2[D3],rotation[D3];
t_real		ALPHA;
{
	t_coord		axe1t[D3],r[D3];

        r[0].x = 1; r[0].y = 0; r[0].z = 0;
        r[1].x = 0; r[1].y = cos(ALPHA); r[1].z = -sin(ALPHA);
        r[2].x = 0; r[2].y = sin(ALPHA); r[2].z = cos(ALPHA);
	geve_oppose_matrix(axe1,axe1t);
	geve_product_matrix(axe1t,r,rotation);
	geve_product_matrix(rotation,axe2,rotation);
	return(OK);
	}

DEFINE	t_err	geve_change_axis(O1,O2,axis,P)
/**********************************************************************
 Compute the coordinates of P in [O2,axe]
	P = m (P - O2) + O1
***********************************************************************/
t_coord		O1,O2,axis[D3],*P;
{
        t_coord p;
        t_real  m[D3][D3];

        if (P == NIL) return(ERROR);
        p = *P;
        geve_substract_vector(&O2,&p,&p);
	m[0][0] = axis[0].x; m[0][1] = axis[0].y; m[0][2] = axis[0].z;
	m[1][0] = axis[1].x; m[1][1] = axis[1].y; m[1][2] = axis[1].z;
	m[2][0] = axis[2].x; m[2][1] = axis[2].y; m[2][2] = axis[2].z;
        geve_matrix_vector(&p,m,&p);
        geve_add_vector(&O1,&p,&p);
        *P = p;
        return(OK);
	}
