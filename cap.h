/*
*	srct.h		Head file for srct.c
*/


#ifndef __SRC_HEAD__
  #define __SRC_HEAD__

#include "inversion.h"

/***********************Constants********************************/

#define RAD	1.745329252e-2	/*degree to rad*/

#define NCP  6		/* 6 segments, 3Pnl 3Sur */

#define STN	1000		/* up to STN stations */
#define NRC	3		/* 3 components of records */
#define NGR	8		/* 8 com. of greens function */


/*********************** Data Structure***************************/

/* focal mechanism (strike, dip, rake) data structure */
typedef struct {
	float	stk;	/* strkie */
	float	dip;	/* dip */
	float	rak;	/* rake */
} MECA;

/* a portion of observed waveform and corresponding 3 components
of Green's functions, cross-correlations, and L2 norms, etc */
typedef struct {
	float	*rec;
	float	*syn;
	int	npt;
	int	on_off;
	int	cc;
	float	b;
	float	w;
	float	rec2;
	float	syn2;
	float	*crl;
	float   ampt;
 	float   shft;
	float	err;
} COMP;

typedef struct {
	char	stn[10];
	float	*rec[NRC];
	float	*grn[NGR];
	float	az;
	float	dist;
	float	alpha;			/* take-off angle */
	float	stla;
	float	stlo;
	COMP	com[6];
} DATA;

typedef struct {
	MECA	meca;
	float	dev[6];			/* uncertainty ellipsis */
	float	err;			/* error for this solution */
	int	cfg[STN][NCP];		/* correlation for each comp. */
	int	shft[STN][NCP];		/* time shift for each comp. */
        float   ampt[STN][NCP];
	float	error[STN][NCP];	/* error for each component */
	int	ms;			/* number of local minimums < 10 */
	int	others[10];		/* top 10 best solutions */
	int	flag;			/* =1 if the best soln is at boundary */
} SOLN;

/* first-motion data */
typedef struct {
	float	az;	/* azimuth */
	float	alpha;	/* take-off angle */
	int	type;	/* 1=P; 2=SV; 3=SH; positive=up; negative=down */
} FM;

/* function declaration */
void	fdtensor(MECA meca,float tensor[3][3]);
float	radpmt(float mom[3][3], FM fm);
void    dc_radiat(double stk,double dip,double rak,double rad[3][3]);
/*SOLN	error(int,DATA *,int,FM *,const int *,float,GRID,float,float,float,float); */
SOLN	error(int,DATA *,int,FM *,const int *,float,GRID,float,float);
void    taper(float *aa, int n);
float	*trap(float, float, float, float *);

#endif
