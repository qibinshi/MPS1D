/****************************************************************

               Multiple Point Sources Inversion Code
	      using Markov Chain Monte Carlo Optimizer
               


  Max. No. of src: 6
  Max. Chain Len.: 250000
                   
  Modify history:
   2016/ 8/ 5 Shengji Wei    -prepare forward calculation code
                              using FK GF and CAP fashion
   2017/ 7/10 Qibin   Shi    -use ENZ(6) instead of rtz(5) 
                              (also modified Structure DATA)
   2017/ 7/15 Qibin   Shi    -read and save data in the memory
                             -untie all six timeshifts
   2017/11/30 Qibin   Shi    -use station-dependent filtering
   2017/12/ 5 Qibin   Shi    -allow different length of GF,
                              maximum npts is MAXNPT
   2018/ 4/23 Qibin   Shi    -read weight file only once
   2018/ 5/ 9 Qibin   Shi    -Joint inversion on regional and
                              teleseismic waveform (Pv&SH),
			      simply rename SH.t with sur.HHE,
			      assuming teleseismic data is NOT
			      sensitive to the source locations 
   2018/ 7/24 Qibin   Shi    -add calibrating time shift into
                              weight file for six components

  References:
  1) Zhao LS, Helmberger DV. 1994. BSSA.
  2) Zhu L, Helmberger DV. 1996. BSSA.
  3) Shi Q, Wei S, Chen M. 2018. GJI.

****************************************************************/

#include <omp.h>
#include <math.h> 
#include <stdio.h>
#include <float.h> 
#include <stdlib.h>
#include <string.h>
#include "nr.h"
#include "cap.h"
#include "sac.h"
#include "nrutil.h"
#include "Complex.h"

#define  NUM_THREADS 2
#define  NRANSI
#define  MAXFLEN 200
#define  CHLEN 250000
#define  NSRC 3150
#define  MAXNPT 10000
#define  MAXDUR 500

int      kd[6]={0,1,2,0,1,2};      // SE, SN, SZ, PE, PN, PZ
int      kk[6]={2,1,0,2,1,0};      // index of synthetic components
char     cm[3]={'E','N','Z'};
float    srcpara[NSRC][CHLEN][9];

int main (int argc, char **argv) {
  int    i,j,k,k1,l,m,nda,npt[STN],plot,useDisp,distkm[STN],durlen,n_shft,thrd;
  int    m_pnl,m_sw,wndw_size,ns,ns_pnl,ns_sw,n_data,io_taper;
  int    t0[STN],te[STN],t1[STN],t2[STN],t3[STN],t4[STN],n1[STN],n2[STN],mm[2],t[6],n[6],max_shft[6];
  int    shft_pnl[3][STN],shft_sur[3][STN],shft0[STN][6],total_n,kmon[STN];
  long int order=4, nsects[STN];
  char	 tmp[280],tmp1[280],path[280],*c_pt[STN],*c_pt1[STN];
  char   type[2]={'B','P'},proto[2]={'B','U'};
  float  x[STN],x1[STN],x2,y[STN],y1,amp,dt,rad[6],arad[STN][3][3],m0[STN],dm,w_cmp[6],rms_cut[5];
  float  pnl_sn[STN][30],pnl_sd[STN][30],sw_sn[STN][30],sw_sd[STN][30],wd1,wd2,sf1,sf2;
  float  *distance,dmin=100.,vp,vs1,vs2,depSqr=110.25,con_shft[STN];
  float  *crlss1[6],*crlss2[6],*f_pt,*cc_pt,on_pt;
  float  *src[STN],*src1[STN],*co[STN],*co1[STN],*syn,*data[STN][3];
  float  bs_body,bs_surf,bs[6],weight,nof_per_samp,tau0,riseTime;
  double f1_pnl,f2_pnl,f1_sw,f2_sw,pie=3.141592654;
  FILE   *f_out,*ff_out;
  SACHEAD *hd0,*hd1;
  DATA   *obs;
  COMP   *spt;
  /* MPS */
  int    npsrc,npara,n_shift,n_corr[STN],ii,z0,sol_shft[STN][6];
  char   tmp_char1[200],tmp_char2[200];
  float  Strike[NSRC],Dip[NSRC],Rake[NSRC],rise[NSRC],slip[NSRC],t_rupture[NSRC];
  float  lat[NSRC],lon[NSRC],dep_src[NSRC],dist_in_fault[NSRC],mu[NSRC],vr[NSRC];
  float  ddd,tp_ref[STN],ts_ref[STN],b_ref[STN],e_ref[STN],dist_ref[STN],tp[STN],ts[STN],*grn[STN],mw,cfg[STN][6],t_shift=200;
  float  az[STN],gcarc[STN],baz[STN],*pt,*pt2,z1,z2,y2,*ppv,*ppr,*ppt,*ppz,*ppe,*ppn,coef;
  float  *syn_mps[STN][3],*syn_tmp2[STN][3],*syn_tmp3[STN][3],total_error=0.,total_m0;
  float  lo1,lo2,la1,la2;
  /* MCMC */
  int    nchain,jj,qq,pp,qp,count=0;
  float  sigma1=4,sigma=1,sigma2,beta=0.998,inter,Mw[6]={0},jump[53]={0};
  float  scale[9]={0.003, 0.03, 0.5, 0.5, 0.05, 0.3, 0.2, 0.3, 0.02};
  float  boundx[4],boundy[4],slop[4],slop1[4];
  double ratio[CHLEN], L[CHLEN];

  /* subroutines */
  void   distaz(float,float,float,float,float *,float *,float *);
  double NormalRandom(double,double,double,double);
  double Normalnew(double,double,double);
  double AverageRandom(double,double);
  double Normal(double,double,double);
  float  *cosine(float, float, int *);



  /* input control parameters */
  scanf("%f%f%f%f",&bs_body,&bs_surf,&x2,&nof_per_samp);
//  scanf("%f%f%f%f",boundx+0,boundx+1,boundx+2,boundx+3);
//  scanf("%f%f%f%f",boundy+0,boundy+1,boundy+2,boundy+3);
  scanf("%f%f%f%f",&lo1,&lo2,&la1,&la2);
  scanf("%f%f%f%f",&vp,&vs1,&vs2,&dt);
  scanf("%f%f%f%f",&wd1,&wd2,&sf1,&sf2);
  scanf("%d%d%d%d",&nchain,&plot,&io_taper,&useDisp);
  scanf("%s",path);

  max_shft[3]=max_shft[4]=max_shft[5]=2*rint(sf1/dt); // max Pnl shift
  max_shft[0]=max_shft[1]=max_shft[2]=2*rint(sf2/dt); // max sur shift
  w_cmp[3]=w_cmp[4]=w_cmp[5]=x2; // relative weight of Pnl
  w_cmp[0]=w_cmp[1]=w_cmp[2]=1.; // relative weight of sur
  bs[3]=bs[4]=bs[5]=bs_body; // distance scale power Pnl
  bs[0]=bs[1]=bs[2]=bs_surf; // distance scale power sur
  mm[0]=rint(wd1/dt); // max Pnl window
  mm[1]=rint(wd2/dt); // max sur window
  n_shift = rint(t_shift/dt);

  /* initial source parameters */
  scanf("%d%f",&npsrc,&ddd);
  npara=npsrc*9-1;
  for(ii=0;ii<npsrc;ii++){
    scanf("%e",lat+ii);
    scanf("%e",lon+ii);
    scanf("%e",dep_src+ii);
    scanf("%e",dist_in_fault+ii);
    scanf("%e",vr+ii);
    scanf("%e",t_rupture+ii);
    scanf("%e",slip+ii);
    scanf("%e",Strike+ii);
    scanf("%e",Dip+ii);
    scanf("%e",Rake+ii);
    scanf("%e",rise+ii);
    scanf("%e",mu+ii);
    srcpara[ii][0][0] = lat[ii];
    srcpara[ii][0][1] = lon[ii];
    srcpara[ii][0][2] = dep_src[ii];
    srcpara[ii][0][3] = t_rupture[ii];
    srcpara[ii][0][4] = log10(slip[ii]);
    srcpara[ii][0][5] = Strike[ii];
    srcpara[ii][0][6] = Dip[ii];
    srcpara[ii][0][7] = Rake[ii];
    srcpara[ii][0][8] = rise[ii];
  }

  scanf("%d",&nda);
  if(nda    > STN  ){fprintf(stderr,"!Max. %d stations!  \n",STN  );nda   =STN  ;}
  if(nchain > CHLEN){fprintf(stderr,"!Max. %d iterations!\n",CHLEN);nchain=CHLEN;}

  beta=pow(0.05,(float)npara/(float)nchain); // cooling rate
  fprintf(stderr,"Cooling rate =  %f\n",beta);

  /* allocate memory */
  obs = (DATA *) malloc(nda*sizeof(DATA));
  if(obs == NULL){fprintf(stderr,"No memory for the whole thing\n");return -1;}

  hd0  = (SACHEAD *) malloc(nda*sizeof(SACHEAD));
  if(hd0  == NULL){fprintf(stderr,"No memory for data head info \n");return -1;}

  hd1  = (SACHEAD *) malloc(nda*sizeof(SACHEAD));
  if(hd1  == NULL){fprintf(stderr,"No memory for synthetics head\n");return -1;}

  distance = (float *) malloc(nda*sizeof(float));
  for(i=0;i<STN;i++) {src[i] = (float *) malloc(MAXDUR*sizeof(float));}
  for(i=0;i<STN;i++) {co[i]  = (float *) malloc((MAXDUR+MAXNPT)*sizeof(float));}
  for(i=0;i<STN;i++) {grn[i] = (float *) malloc(MAXNPT*sizeof(float));}
  for(i=0;i<STN;i++) {for(k=0; k<3; k++) syn_mps[i][k] =(float *) calloc(MAXNPT+2*n_shift,sizeof(float));}
  for(i=0;i<STN;i++) {for(k=0; k<3; k++) syn_tmp2[i][k]=(float *) calloc(MAXNPT,sizeof(float));}
  for(i=0;i<STN;i++) {for(k=0; k<3; k++) syn_tmp3[i][k]=(float *) calloc(MAXNPT,sizeof(float));}
  for(k=0; k<3; k++) crlss1[k] = (float *) malloc(4*mm[1]*sizeof(float));
  for(k=0; k<3; k++) crlss2[k] = (float *) malloc(4*mm[1]*sizeof(float));
  for(k=3;k<NCP;k++) crlss1[k] = (float *) malloc(4*mm[0]*sizeof(float));
  for(k=3;k<NCP;k++) crlss2[k] = (float *) malloc(4*mm[0]*sizeof(float));



  /* read sta-info and data */
  for(i=0;i<nda;i++) {

    scanf("%s",(&(obs[i]))->stn);                              // station name
    scanf("%d",kmon+i);                   // km switch for reading grn
    for(j=0;j<6;j++) scanf("%d",&(&(obs[i]))->com[5-j].on_off);      // weight
    scanf("%d",shft_pnl[0]+i);               // calibrating time shift
    scanf("%d",shft_pnl[1]+i);
    scanf("%d",shft_pnl[2]+i);
    scanf("%d",shft_sur[0]+i);
    scanf("%d",shft_sur[1]+i);
    scanf("%d",shft_sur[2]+i);
    scanf("%lf%lf%lf%lf",&f1_pnl,&f2_pnl,&f1_sw,&f2_sw);// corner freq
    design(order,type,proto,1.,1.,f1_pnl,f2_pnl,(double)dt,pnl_sn[i],pnl_sd[i],&nsects[i]);
    design(order,type,proto,1.,1.,f1_sw, f2_sw, (double)dt,sw_sn[i], sw_sd[i], &nsects[i]);

    strcat(strcat(strcat(strcpy(tmp,argv[1]),"/"),(&(obs[i]))->stn),".HHE");
    c_pt[i] = strrchr(tmp,(int) 'E');
    for(n_data=0;n_data<NRC;n_data++){
      *(c_pt[i]) = cm[n_data];
      if((data[i][n_data] = read_sac(tmp,(&(hd0[i])))) == NULL){
        fprintf(stderr,"Skip reading %s !\n",tmp);
        goto next_sta;
      }
    }
    (&(obs[i]))->dist = distance[i] = (&(hd0[i]))->dist;
    (&(obs[i]))->az   = (&(hd0[i]))->az;
    (&(obs[i]))->stla = (&(hd0[i]))->stla;
    (&(obs[i]))->stlo = (&(hd0[i]))->stlo;

    x[i]  = (&(hd0[i]))->b-(&(hd0[i]))->o;
    y[i]  = (&(hd0[i]))->e-(&(hd0[i]))->o;
    x1[i] = (&(hd0[i]))->a-(&(hd0[i]))->o;
    t1[i] = rint(((&(hd0[i]))->t1-(&(hd0[i]))->b)/dt);
    t2[i] = rint(((&(hd0[i]))->t2-(&(hd0[i]))->b)/dt);
    t3[i] = rint(((&(hd0[i]))->t3-(&(hd0[i]))->b)/dt);
    t4[i] = rint(((&(hd0[i]))->t4-(&(hd0[i]))->b)/dt); 

    goto normal;
    next_sta:
      nda--; i--;
    normal:
      fprintf(stderr,"Read station NO.%3d   %s\n",i,(&(obs[i]))->stn);
  }

  /* Markov Chain Monte Carlo */
  ratio[0]=1;
  ff_out=fopen("pdf.out","w");

  for(jj=0;jj<nchain;jj++) {
    for(ii=0;ii<npsrc;ii++){                // bounds of source parameters
      if     (srcpara[ii][jj][0]>la2)  srcpara[ii][jj][0] = la2;
      else if(srcpara[ii][jj][0]<la1)  srcpara[ii][jj][0] = la1;
      if     (srcpara[ii][jj][1]>119.9-srcpara[ii][jj][0]/8)  
	      srcpara[ii][jj][1] = 119.9-srcpara[ii][jj][0]/8; // for palu
      else if(srcpara[ii][jj][1]<119.7-srcpara[ii][jj][0]/8)  
	      srcpara[ii][jj][1] = 119.7-srcpara[ii][jj][0]/8;
//      if     (srcpara[ii][jj][1]>lo2)  srcpara[ii][jj][1] = lo2;
//      else if(srcpara[ii][jj][1]<lo1)  srcpara[ii][jj][1] = lo1;
      if     (srcpara[ii][jj][2]>20. )  srcpara[ii][jj][2] = 20. ;
      else if(srcpara[ii][jj][2]<2.  )  srcpara[ii][jj][2] = 2.  ;
      if     (srcpara[ii][jj][3]<0.  )  srcpara[ii][jj][3] = 0.  ;
      if     (srcpara[ii][jj][4]>2.  )  srcpara[ii][jj][4] = 2.  ;
      else if(srcpara[ii][jj][4]<-1.5)  srcpara[ii][jj][4] = -1.5;
      if     (srcpara[ii][jj][6]>90. ) {
        srcpara[ii][jj][6] = 180. - srcpara[ii][jj][6];
        srcpara[ii][jj][7] = 0.   - srcpara[ii][jj][7];
        srcpara[ii][jj][5] = srcpara[ii][jj][5] - 180.;
      }else if(srcpara[ii][jj][6]<0.) {
        srcpara[ii][jj][6] = 0.   - srcpara[ii][jj][6];
        srcpara[ii][jj][7] = srcpara[ii][jj][7] - 180.;
        srcpara[ii][jj][5] = srcpara[ii][jj][5] - 180.;
      }
      if     (srcpara[ii][jj][5]>360. )  srcpara[ii][jj][5] -= 360.;
      else if(srcpara[ii][jj][5]<0.   )  srcpara[ii][jj][5] += 360.;
      if     (srcpara[ii][jj][5]>25.  && srcpara[ii][jj][5]<90. )  
	      srcpara[ii][jj][5] = 25. ;
      else if(srcpara[ii][jj][5]>90.  && srcpara[ii][jj][5]<155.)  
	      srcpara[ii][jj][5] = 155.;
      else if(srcpara[ii][jj][5]>205. && srcpara[ii][jj][5]<270.)  
	      srcpara[ii][jj][5] = 205.;
      else if(srcpara[ii][jj][5]>270. && srcpara[ii][jj][5]<335.)  
	      srcpara[ii][jj][5] = 335.;
      if     (srcpara[ii][jj][7]>180. )  srcpara[ii][jj][7] -= 360.;
      else if(srcpara[ii][jj][7]<-180.)  srcpara[ii][jj][7] += 360.;
      if     (srcpara[ii][jj][8]<0.2  )  srcpara[ii][jj][8]  = 0.2 ; 
      else if(srcpara[ii][jj][8]>10.  )  srcpara[ii][jj][8]  = 10. ;

      lat[ii]       = srcpara[ii][jj][0];
      lon[ii]       = srcpara[ii][jj][1];
      dep_src[ii]   = srcpara[ii][jj][2];
      t_rupture[ii] = srcpara[ii][jj][3];
      slip[ii]      = pow(10.0, srcpara[ii][jj][4]);
      Strike[ii]    = srcpara[ii][jj][5];
      Dip[ii]       = srcpara[ii][jj][6];
      Rake[ii]      = srcpara[ii][jj][7];
      rise[ii]      = srcpara[ii][jj][8];
    }


// OMP
omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for private (ii,j,l,k,k1,m,y1,tmp_char1,tmp_char2,coef,pt,pt2,durlen,ppv,ppr,ppt,ppz,ppn,ppe,weight,wndw_size,x2,f_pt,cc_pt,on_pt) reduction (+:total_error,total_n,n_shft)

    /* parallel loop over stations */
    for(i=0;i<nda;i++) {
	fprintf(stderr,"%5.1f  %5.1f   %s\n",(&(obs[i]))->stla,(&(obs[i]))->stlo,(&(obs[i]))->stn);
//	fprintf(stderr,"i= %3d ---sta: %s\n",i,(&(obs[i]))->stn);
      /* compute MPS synthetics */
      for(ii=0;ii<npsrc;ii++){
//        fprintf(stderr,"i= %3d ---ii= %3d ---sta: %s\n",i,ii,(&(obs[i]))->stn);
        distaz((&(obs[i]))->stla,(&(obs[i]))->stlo,lat[ii],lon[ii],gcarc+i,baz+i,az+i);
//fprintf(stderr,"%5.1f  %5.1f   %5.1f  %5.1f   %4.1f  %s\n",(&(obs[i]))->stla,(&(obs[i]))->stlo,lat[ii],lon[ii],gcarc[i],(&(obs[i]))->stn);
        strcat(strcpy(tmp_char1,path),"_");

        if (dep_src[ii]+ddd/2>99.9){
          sprintf(tmp_char2,"%3.0f",floor((dep_src[ii]+ddd/2)/ddd)*ddd);  
        }else if (dep_src[ii]+ddd/2>=10.0){
          sprintf(tmp_char2,"%2.0f",floor((dep_src[ii]+ddd/2)/ddd)*ddd);
        }else if (dep_src[ii]>0.){
          sprintf(tmp_char2,"%1.0f",floor((dep_src[ii]+ddd/2)/ddd)*ddd);
        }else{
          fprintf(stderr,"Depth is negative! %f\n",dep_src[ii]);
        }

        strcat(strcat(tmp_char1,tmp_char2),"/");

	if(kmon[i]>0){
          distkm[i] = gcarc[i]*pie*6371/180;
          if(distkm[i]<1) distkm[i]=1;
          sprintf(tmp_char2,"%d",distkm[i]);
	}else{
          sprintf(tmp_char2,"%4.1f",gcarc[i]);
	}

        strcat(strcat(tmp_char1,tmp_char2),".grn.0");
        c_pt[i] = strrchr(tmp_char1,(int) '.') + 1;

        dc_radiat(az[i]-Strike[ii],Dip[ii],Rake[ii],arad[i]);
        m0[i] = mu[ii] * slip[ii] * 1.0e-13;

        for(l=0; l<3; l++){ // nine components
	  for(j=0; j<3; j++){

	    coef = m0[i] * arad[i][l][j];

	    if((grn[i]=read_sac(tmp_char1,(&(hd1[i])))) == NULL) continue;

            if(ii==0 && l==0 && j==0){ // 1st source's head

                     npt[i] = (&(hd1[i]))->npts;
	           b_ref[i] = (&(hd1[i]))->b;
                   e_ref[i] = (&(hd1[i]))->e;
                  tp_ref[i] = (&(hd1[i]))->t1;
                  ts_ref[i] = (&(hd1[i]))->t2;
	        dist_ref[i] = (&(hd1[i]))->dist;
   	         (&(obs[i]))->alpha = (&(hd1[i]))->user1;
                  n_corr[i] = rint(t_rupture[ii]/dt) + n_shift;

              for(k=0; k<3; k++) {
                for(k1=0,pt=syn_mps[i][k];k1<npt[i]+2*n_shift;k1++,pt++) *pt=0;
              }

            }else if(l==0 && j==0){ // following sources' head
              if(npt[i] < (&(hd1[i]))->npts) npt[i] = (&(hd1[i]))->npts; // npt is the largest
              tp[i] = (&(hd1[i]))->t1;
              ts[i] = (&(hd1[i]))->t2;
              n_corr[i] = rint((tp[i] - tp_ref[i] + t_rupture[ii])/dt) + n_shift;
            }

            for(pt=syn_tmp2[i][j],k1=0;k1<((&(hd1[i]))->npts);k1++,pt++){
	      (*pt) += coef*grn[i][k]; grn[i][k]=0;
	    }

            (*(c_pt[i]))++;

          }
        } // finish nine components

        durlen = rint(2*rise[ii]/dt);
	src[i] = (float *) realloc(src[i], (durlen+1)*sizeof(float));
        co[i]  = (float *) realloc(co[i],(durlen+npt[i])*sizeof(float));
        trap(2*rise[ii], rise[ii], dt, &durlen, src[i]);

        for(k=0;k<3;k++) conv(src[i],durlen,syn_tmp2[i][k],npt[i],co[i]); // convolution

        src1[i] = (float *) realloc(src[i],MAXDUR*sizeof(float));
        if(!src1[i]) {fprintf(stderr,"Lost memory for src");}else{src[i]=src1[i];}
        co1[i]  = (float *) realloc(co[i],(MAXDUR+MAXNPT)*sizeof(float));
        if(!co1[i])  {fprintf(stderr,"Lost memory for co");}else{co[i]=co1[i];}

	if(kmon[i]>0){ // vrt to ZNE for regional stations (each source)
          for(ppv=syn_tmp2[i][0],
	      ppr=syn_tmp2[i][1],
	      ppt=syn_tmp2[i][2],
	      ppz=syn_tmp3[i][0],
	      ppn=syn_tmp3[i][1],
	      ppe=syn_tmp3[i][2],
	      j=0;j<npt[i];j++,
	      ppv++,ppr++,ppt++,
	      ppz++,ppe++,ppn++){

	    (*ppz)=(*ppv);
	    (*ppn)=(*ppr)*cos((baz[i]+180)*pie/180)-(*ppt)*sin((baz[i]+180)*pie/180);
            (*ppe)=(*ppr)*sin((baz[i]+180)*pie/180)+(*ppt)*cos((baz[i]+180)*pie/180);
            (*ppv)=0; // clean syn_tmp2
            (*ppr)=0;
            (*ppt)=0;
	 
	  }
	  for(k=0;k<3;k++){
            for(pt=syn_tmp3[i][k],pt2=syn_mps[i][k]+n_corr[i],j=0;j<npt[i];j++,pt2++,pt++){
              (*pt2) += (*pt); (*pt) = 0; // clean syn_tmp3
            }
          }
        }else{ // for teleseismic use vrt
          for(k=0;k<3;k++){
            for(pt=syn_tmp2[i][k],pt2=syn_mps[i][k]+n_corr[i],j=0;j<npt[i];j++,pt2++,pt++){
              (*pt2) += (*pt); (*pt) = 0; // clean syn_tmp2
            }
          }
	}

      } // finish MPS waveform

      /* update header of synthetics */
      (&(hd1[i]))->b    = b_ref[i]  - t_shift; // add more points before the beginning
      (&(hd1[i]))->e    = e_ref[i]  + t_shift; // add more points after the end
      (&(hd1[i]))->t1   = tp_ref[i] + t_rupture[0];
      (&(hd1[i]))->t2   = ts_ref[i] + t_rupture[0];
      (&(hd1[i]))->npts = npt[i]    + 2*n_shift;
      (&(hd1[i]))->dist = dist_ref[i];

      /* use first-arrival (if picked) to correct synthetics */
      con_shft[i] = 0;
      if(x1[i]>0.) con_shft[i] += x1[i] - (&(hd1[i]))->t1;
      t0[i]=rint((x[i]-con_shft[i]-(&(hd1[i]))->b)/dt) + 1;
      te[i]=rint((y[i]-con_shft[i]-(&(hd1[i]))->b)/dt) - 1;

      /* Pnl window (relative to b of syn) */
      if(t1[i] > 0 && t2[i] > 0){ // use picked t1 t2
        t1[i]+=t0[i];
        t2[i]+=t0[i]; 
      }else{
        if(vp > 0.) // determine t1 by Vp
	  t1[i]=rint((sqrt(distance[i]*distance[i]+depSqr)/vp-(&(hd1[i]))->b)/dt-0.3*mm[0]);	
        else // determine t1 from syn
	  t1[i]=rint(((&(hd1[i]))->t1-(&(hd1[i]))->b)/dt-0.1*mm[0]);	
        t2[i]  =rint(((&(hd1[i]))->t2-(&(hd1[i]))->b)/dt+0.2*mm[0]); // determine t2 from syn (after ts)
      }

      if(t1[i]<t0[i]) t1[i]=t0[i]; // let t1 > b(data)
      if(t1[i]<shft_pnl[0][i]) t1[i]=shft_pnl[0][i]; // let t1 > b(syn)
      if(t1[i]<shft_pnl[1][i]) t1[i]=shft_pnl[1][i];
      if(t1[i]<shft_pnl[2][i]) t1[i]=shft_pnl[2][i];
      if(t2[i]>te[i]) t2[i]=te[i]; // let t2 < e(data)
      if(t2[i]>(&(hd1[i]))->npts+shft_pnl[0][i]) t2[i]=(&(hd1[i]))->npts+shft_pnl[0][i]; // let t2 < e(syn)
      if(t2[i]>(&(hd1[i]))->npts+shft_pnl[1][i]) t2[i]=(&(hd1[i]))->npts+shft_pnl[1][i];
      if(t2[i]>(&(hd1[i]))->npts+shft_pnl[2][i]) t2[i]=(&(hd1[i]))->npts+shft_pnl[2][i];

      /* Sur window (relative to b of syn) */
      if(t3[i] > 0 && t4[i] >0){ // use picked t3 t4
        t3[i]+=t0[i];
        t4[i]+=t0[i];
      }else{
        if(vs1 >0. && vs2 > 0.){ // determine t3 t4 by Vlove and Vrayleigh
	  t3[i]=rint((sqrt(distance[i]*distance[i]+depSqr)/vs1-(&(hd1[i]))->b)/dt - 0.3*mm[1]);
	  t4[i]=rint((sqrt(distance[i]*distance[i]+depSqr)/vs2-(&(hd1[i]))->b)/dt + 0.7*mm[1]);
        }else{ // determine t3 t4 from syn
          t3[i]=rint(((&(hd1[i]))->t2-(&(hd1[i]))->b)/dt-0.2*mm[1]);
          t4[i]=t3[i]+mm[1];
        }
      }

      if(t3[i]<t1[i]) t3[i]=t1[i]; // let t3 > t1
      if(t3[i]<t1[i]-shft_pnl[0][i]+shft_sur[0][i]) t3[i]=t1[i]-shft_pnl[0][i]+shft_sur[0][i];
      if(t3[i]<t1[i]-shft_pnl[1][i]+shft_sur[1][i]) t3[i]=t1[i]-shft_pnl[1][i]+shft_sur[1][i];
      if(t3[i]<t1[i]-shft_pnl[2][i]+shft_sur[2][i]) t3[i]=t1[i]-shft_pnl[2][i]+shft_sur[2][i];
      if(t4[i]>te[i]) t4[i]=te[i]; // let t4 < e(data)
      if(t4[i]>(&(hd1[i]))->npts+shft_sur[0][i]) t4[i]=(&(hd1[i]))->npts+shft_sur[0][i];
      if(t4[i]>(&(hd1[i]))->npts+shft_sur[1][i]) t4[i]=(&(hd1[i]))->npts+shft_sur[1][i];
      if(t4[i]>(&(hd1[i]))->npts+shft_sur[2][i]) t4[i]=(&(hd1[i]))->npts+shft_sur[2][i];

      /* windows length **/
      n1[i] = t2[i] - t1[i];
      n2[i] = t4[i] - t3[i];
      if(n1[i]>mm[0]) n1[i]=mm[0];
      if(n2[i]>mm[1]) n2[i]=mm[1];

      /* discard too short windows */
      if(n1[i]<max_shft[3]){
        fprintf(ff_out,"Shift less for P (Pnl) of %s !\n", (&(obs[i]))->stn);
        max_shft[3]=n1[i];
      }
      if(n2[i]<max_shft[0]){
        fprintf(ff_out,"Shift less for S (sur) of %s !\n", (&(obs[i]))->stn);
        max_shft[0]=n2[i];
      }

      /* cut */
      t[0]=t[1]=t[2]=t3[i]; // start of Sur components   
      t[3]=t[4]=t[5]=t1[i]; // start of Pnl components
      n[0]=n[1]=n[2]=n2[i]; // length of Sur components
      n[3]=n[4]=n[5]=n1[i]; // length of Pnl components

      shft0[i][0]=shft_sur[2][i]; // time calibration for Sur 
      shft0[i][1]=shft_sur[1][i];
      shft0[i][2]=shft_sur[0][i];
      shft0[i][3]=shft_pnl[2][i]; // time calibration for Pnl
      shft0[i][4]=shft_pnl[1][i];
      shft0[i][5]=shft_pnl[0][i];

      if((&(obs[i]))->com[0].on_off>0 || (&(obs[i]))->com[1].on_off>0 || (&(obs[i]))->com[2].on_off>0) n_shft++;
      if((&(obs[i]))->com[3].on_off>0 || (&(obs[i]))->com[4].on_off>0 || (&(obs[i]))->com[5].on_off>0) n_shft++;

      for(spt=(&(obs[i]))->com,j=0;j<NCP;j++,spt++){
        
	if(spt->on_off) total_n+=npt[i];
        spt->npt = npt[i] = n[j];
        spt->b   = t[j]*dt+con_shft[i]+(&(hd1[i]))->b; 
        weight = w_cmp[j]*pow(distance[i]/dmin,bs[j]);
        wndw_size = npt[i]*sizeof(float);
	/* data */
        if(jj==0){

          spt->rec = (float *) malloc(wndw_size);
          spt->syn = (float *) malloc(wndw_size);
          memcpy(spt->rec, data[i][kd[j]] + t[j] - t0[i], wndw_size);
          if(j<3){apply(spt->rec,(long int) npt[i],0, sw_sn[i], sw_sd[i],nsects[i]);}
          else   {apply(spt->rec,(long int) npt[i],0,pnl_sn[i],pnl_sd[i],nsects[i]);}
          if(useDisp) cumsum(spt->rec, npt[i], dt); // integrate velocity to displacement
          if(io_taper) taper(spt->rec, npt[i]);

          for(x2=0.,f_pt=spt->rec,l=0;l<npt[i];l++,f_pt++){
	    *f_pt *= weight;
	    x2+=(*f_pt)*(*f_pt); // square data for L2 norm
          }  
          spt->rec2 = x2;

        }
	/* syn */
        f_pt = spt->syn;
        memcpy(f_pt, syn_mps[i][kk[j]] + t[j] - shft0[i][j], wndw_size);
        if(j<3){apply(f_pt,(long int) npt[i],0, sw_sn[i], sw_sd[i],nsects[i]);}
	else   {apply(f_pt,(long int) npt[i],0,pnl_sn[i],pnl_sd[i],nsects[i]);}
        if(useDisp) cumsum(f_pt, npt[i], dt);
        if(io_taper) taper(f_pt, npt[i]);

        for(x2=0.,l=0;l<npt[i];l++){
	  f_pt[l] *= weight;
          x2+=f_pt[l]*f_pt[l]; // square syn for L2 norm
	}
	spt->syn2 = x2;

        /* CC */
        spt->crl = crscrl(npt[i],spt->rec,f_pt,max_shft[j],crlss1[j],crlss2[j]);
        if(j<3){
          crlss1[j] = (float *) realloc(crlss1[j],4*mm[1]*sizeof(float));
          crlss2[j] = (float *) realloc(crlss2[j],4*mm[1]*sizeof(float));
        }else{
          crlss1[j] = (float *) realloc(crlss1[j],4*mm[0]*sizeof(float));
          crlss2[j] = (float *) realloc(crlss2[j],4*mm[0]*sizeof(float));
        }
      }
fprintf(stderr,"i= %3d ---ii= %3d ---sta: %s\n",i,ii,(&(obs[i]))->stn);
      /* find time shifts */
      for(spt=(&(obs[i]))->com,j=0;j<NCP;j++,spt++){
        cc_pt=spt->crl;
	on_pt=spt->on_off>0?1:0;
	for(y1=-FLT_MAX,l=0;l<=max_shft[j];l++){
          x[i] = *cc_pt++;
          y[i] = on_pt*x[i];
          if(y[i]>y1){
            y1     = y[i];
	    cfg[i][j] = x[i];
	    m      = l;
	  }
        }
        sol_shft[i][j] = m - max_shft[j]/2;
      }
      spt--;

      /* error, cc, time shift, amplitude ratio */
      for(j=5;j>=0;j--,spt--){
        spt->err     = spt->on_off * (spt->rec2 + spt->syn2 - 2.* cfg[i][j]); // L2 norm
	spt->cc      = rint(100 * cfg[i][j] / sqrt(spt->rec2 * spt->syn2));
	spt->shft    = dt * (shft0[i][j] + sol_shft[i][j]);
        spt->ampt    = sqrt(spt->syn2 / spt->rec2);  
        total_error += spt->err; 
      }

      /* output waveforms */
      if(jj == nchain-1){
        strcat(strcat(strcat(strcat(strcat(strcpy(tmp,argv[1]),"/"),"model_mps"), "_"),(&(obs[i]))->stn),".11");
        c_pt[i] = strrchr(tmp,(int) '1');
        c_pt1[i]= c_pt[i]-1;
        for(spt=(&(obs[i]))->com,j=0;j<NCP;j++,spt++){
          npt[i]    = spt->npt;
          (*(&(hd1[i])))    = sachdr(dt, npt[i], spt->b); // add b, dt, npt into head
          (&(hd1[i]))->dist = (&(obs[i]))->dist;
	  (&(hd1[i]))->az   = (&(obs[i]))->az;
	  (&(hd1[i]))->user1= (&(obs[i]))->alpha;
          (&(hd1[i]))->a    = (&(hd1[i]))->b;
          write_sac(tmp,(*(&(hd1[i]))),spt->rec);
          (*(c_pt1[i]))++;
          (&(hd1[i]))->b   -= (shft0[i][j]*dt + con_shft[i]);
          (&(hd1[i]))->a    = (&(hd1[i]))->b - sol_shft[i][j]*dt;
          write_sac(tmp,(*(&(hd1[i]))),spt->syn);
          (*(c_pt[i]))++;
          (*(c_pt1[i]))--;
        }
      }

      if(jj==0){
        for(j=0;j<NRC;j++) {free(data[i][j]);}
      }

    } // end parallel loop over stations

    /* goodness */
    L[jj] = Normalnew(total_error,0.,sigma1);
    if(jj > 0 ) {
      if(jj<nchain/5)         {ratio[jj] = 1  * jj * (L[jj]/L[jj-1] - 1) + 1;}
      else if(jj<nchain*2/5)  {ratio[jj] = 1.5* jj * (L[jj]/L[jj-1] - 1) + 1;}
      else if(jj<nchain*3/5)  {ratio[jj] = 2  * jj * (L[jj]/L[jj-1] - 1) + 1;}
      else if(jj<nchain*4/5)  {ratio[jj] = 3  * jj * (L[jj]/L[jj-1] - 1) + 1;}
      else if(jj<nchain*9/10) {ratio[jj] = 4 * jj * (L[jj]/L[jj-1] - 1) + 1;}
      else                    {ratio[jj] = 8 * jj * (L[jj]/L[jj-1] - 1) + 1;}

      if(ratio[jj] < AverageRandom(0.,1.)){
        L[jj] = L[jj-1];
	for(ii=0;ii<npsrc;ii++) {
          for(qq=0;qq<9;qq++) {

            srcpara[ii][jj][qq] = srcpara[ii][jj-1][qq];

	  }
        }
      }else{
        count++;
      }        
    }

    /* next MCMC sample */
    if(jj<nchain-1){
      for(qq=0;qq<npara;qq++) jump[qq] = 0.;
      qp = jj%npara;
      sigma2   = sigma * pow(beta, jj/npara);
      inter    = 6 * cos(pie*jj/(nchain*2.5));
      jump[qp] = NormalRandom(0.,sigma2,(0-inter),inter);
      for(qq=0;qq<3;qq++)
	srcpara[0][jj+1][qq] = jump[qq]  *scale[qq]+srcpara[0][jj][qq];
//        srcpara[0][jj+1][qq] = srcpara[0][jj][qq];//fix location1
      for(qq=4;qq<9;qq++) 
	srcpara[0][jj+1][qq] = jump[qq-1]*scale[qq]+srcpara[0][jj][qq];
      for(ii=1;ii<npsrc;ii++){
        for(qq=0;qq<9;qq++)
	  srcpara[ii][jj+1][qq] = jump[qq+ii*9-1]*scale[qq]+srcpara[ii][jj][qq]; 
      }
    }

    for(ii=0;ii<npsrc;ii++){
      fprintf(ff_out,"%7.4f %8.4f %5.2f %5.2f %8.4f %6.2f %5.2f %7.2f %5.2f ",
                   lat[ii],
                   lon[ii],
                   dep_src[ii],
                   t_rupture[ii],
                   slip[ii],
                   Strike[ii],
                   Dip[ii],
                   Rake[ii],
                   rise[ii]);
    }
    fprintf(ff_out,"%.5f %.5f\n",L[jj],ratio[jj]);
  } // end MCMC
  fclose(ff_out);


  /* final output */
  total_m0 = 0;
  for(ii=0;ii<npsrc;ii++) {
    total_m0 += mu[ii]*slip[ii];
    Mw[ii] = (log10(mu[ii]*slip[ii]) - 9.1)*2./3.;
  }
  mw = (log10(total_m0) - 9.1)*2./3.;
  fprintf(stderr,"Mw %8.2f\n",mw);

  strcat(strcat(strcat(strcpy(tmp1,argv[1]),"/"),"model_mps"),".out");
  f_out=fopen(tmp1,"w");
  fprintf(f_out,"Event %s Model MPS ",argv[1]);
  for(ii=0;ii<npsrc;ii++)
    fprintf(f_out,"FM %3.0f %3.0f %3.0f Dp %2.0f Mw %5.2f Delay %4.1f dura %4.1f ",
		    Strike[ii],
		    Dip[ii],
		    Rake[ii],
		    dep_src[ii],
		    Mw[ii],
		    t_rupture[ii],
		    2*rise[ii]);
  fprintf(f_out,"rms %9.2e loc ",total_error);
  for(ii=0;ii<npsrc;ii++) fprintf(f_out,"%8.5f %8.4f ",lat[ii],lon[ii]);
  fprintf(f_out,"\n");
  for(i=0;i<nda;i++) {
    fprintf(f_out,"%-9s %5.2f ",(&(obs[i]))->stn,con_shft[i]);
    for(j=5,spt=(&(obs[i]))->com+5;j>=0;j--,spt--){
      fprintf(f_out,"%1d %8.2e %2d %5.2f %5.2f ",
		      spt->on_off,
		      spt->err,
		      spt->cc,
		      spt->ampt,
		      spt->shft);
    }
    fprintf(f_out,"\n");
  }
  fclose(f_out);

  if ( ! plot ) return 0;
  return 0;
} // end inversion program

/***************************************************************/


void    taper(float *aa, int n)
{
  int	i, m;
  float	tt, pi1;
  m = rint(0.3*n);
  pi1 = 3.1415926/m;
  for (i=0; i<m; i++) {
    tt = 0.5*(1.-cos(i*pi1));
    aa[i] *= tt;
    aa[n-i-1] *= tt;
  }
}

int discard_bad_data(int nda, DATA *obs, SOLN sol, float sig, float rms_cut[]) {
   int i, j, n;
   COMP	*spt;
   n = 0;
   for(i=0;i<nda;i++,obs++) {
     spt = obs->com;
     for(j=0; j<5; j++,spt++) {
        if (sol.error[i][j]/sig > rms_cut[j]) {
	   spt->on_off = 0;
	   n++;
	}
     }
   }
   return(n);
}

float *trap(float t1, float t2, float dt, int *n, float *s) {
    int i, n1, n2;
    float slope;
    t1=t1-t2;
    n1 = rint(t1/dt); if (n1<1) n1 = 1;
    n2 = rint(t2/dt); if (n2<1) n2 = 1;
    if (n1 > n2) {
	i = n1;
	n1 = n2;
	n2 = i;
    }
    slope = 1./(n1*n2);
    s[0] = 0;
    for(i=1;i<=n1;i++) s[*n-i] = s[i]=s[i-1] + slope;
    for(;i<n2;i++) s[i]=s[i-1];
    s[*n]=0;

    return s;
}

void buttbp(float *h, int *m, float *g, float fl, float fh) {
    int n;
    float eps, ah, al, ap, as, fs;
//    SACHEAD hd;
//    float gout[512],pout[512];
//    void recres_();
    void butpas_(float *,int *,float *,int *,float *,float *,float *,float *,float *);
    eps = 0.3;
    ah = 0.1;
    al = 0.1;
    ap=sqrt(ah/(1.-ah));
    as=sqrt((1.-al)/al);
    fs=fh*(1.+eps); if (fs>0.5) fs=0.5;
    butpas_(h,m,g,&n,&fl,&fh,&fs,&ap,&as);
//    fprintf(stderr,"butterworth filter order = %d\n",*m);
//    fl=0;
//    n=512;
//    fh=0.5/n;
//    recres_(h,m,g,&fl,&fh,gout,pout,&n);
//    fh=fh/0.2;
//    hd = sachdr(fh,n,fl);
//    write_sac("am.sac",hd,gout);
//    write_sac("ph.sac",hd,pout);
}

void principal_values(float a[]) {
   int i;
   float **b, **v;
   b = matrix(1,3,1,3);
   v = matrix(1,3,1,3);
   for(i=0;i<3;i++) b[i+1][i+1]=a[i];
   b[1][2]=b[2][1]=a[3];b[1][3]=b[3][1]=a[4];b[2][3]=b[3][2]=a[5];
   jacobi(b,3,a,v,&i);
   eigsrt(a,v,3);
   for(i=0;i<3;i++) {a[i] = a[i+1]; if (a[i]<0.0001*a[1]) a[i]=0.0001*a[1];}
   free_convert_matrix(b,1,3,1,3);
   free_matrix(v,1,3,1,3);
}

void distaz(float stla,float stlo, float evla, float evlo, float *gcarc, float *baz, float *az){
/*
 * c
 * c Subroutine to calculate the Great Circle Arc distance
 * c    between two sets of geographic coordinates
 * c
 * c Given:  stla => Latitude of first point (+N, -S) in degrees
 * c         stlo => Longitude of first point (+E, -W) in degrees
 * c         evla => Latitude of second point
 * c         evlo => Longitude of second point
 * c
 * c Returns:  gcarc => Great Circle Arc distance in degrees
 * c           az    => Azimuth from pt. 1 to pt. 2 in degrees
 * c           baz   => Back Azimuth from pt. 2 to pt. 1 in degrees
 * c
 * c If you are calculating station-epicenter pairs, pt. 1 is the station
 * c
 * c Equations take from Bullen, pages 154, 155
 * c
 * c T. Owens, September 19, 1991
 * c           Sept. 25 -- fixed az and baz calculations
 * c
 *   P. Crotwell, Setember 27, 1994
 *               Converted to c to fix annoying problem of fortran giving wrong
 *                              answers if the input doesn't contain a decimal point.
 *                              */
//  double gcarc, az, baz;
//  double scolat, slon, ecolat, elon;
//  double a,b,c,d,e,aa,bb,cc,dd,ee,g,gg,h,hh,k,kk;
//  double rhs1,rhs2,sph,rad,del,daz,dbaz,pi,piby2;
  float scolat, slon, ecolat, elon;
  float a,b,c,d,e,aa,bb,cc,dd,ee,g,gg,h,hh,k,kk;
  float rhs1,rhs2,sph,rad,del,daz,dbaz,pi,piby2;
  
  pi=3.141592654;
  piby2=pi/2.0;
  rad=2.*pi/360.0;

/*
 * c
 * c scolat and ecolat are the geocentric colatitudes
 * c as defined by Richter (pg. 318)
 * c
 * c Earth Flattening of 1/298.257 take from Bott (pg. 3)
 * c
 * */
   sph=1.0/298.257;
   scolat=piby2 - atan((1.-sph)*(1.-sph)*tan(stla*rad));
   ecolat=piby2 - atan((1.-sph)*(1.-sph)*tan(evla*rad));
   slon=stlo*rad;
   elon=evlo*rad;
/*
 * c
 * c  a - e are as defined by Bullen (pg. 154, Sec 10.2)
 * c     These are defined for the pt. 1
 * c
 * */
   a=sin(scolat)*cos(slon);
   b=sin(scolat)*sin(slon);
   c=cos(scolat);
   d=sin(slon);
   e=-cos(slon);
   g=-c*e;
   h=c*d;
   k=-sin(scolat);
/*
 *  c
 *  c  aa - ee are the same as a - e, except for pt. 2
 *  c
 *  */
   aa=sin(ecolat)*cos(elon);
   bb=sin(ecolat)*sin(elon);
   cc=cos(ecolat);
   dd=sin(elon);
   ee=-cos(elon);
   gg=-cc*ee;
   hh=cc*dd;
   kk=-sin(ecolat);
/*
 * c
 * c  Bullen, Sec 10.2, eqn. 4
 * c
 * */
   del=acos(a * aa + b * bb + c * cc);
   *gcarc=del/rad;
/*
 * c
 * c  Bullen, Sec 10.2, eqn 7 / eqn 8
 * c
 * c    pt. 1 is unprimed, so this is technically the baz
 * c
 * c  Calculate baz this way to avoid quadrant problems
 * c
 * */
   rhs1=(aa-d)*(aa-d)+(bb-e)*(bb-e)+cc*cc - 2.;
   rhs2=(aa-g)*(aa-g)+(bb-h)*(bb-h)+(cc-k)*(cc-k) - 2.;
   dbaz=atan2(rhs1,rhs2);
   if (dbaz<0.0) {
      dbaz=dbaz+2*pi;
   } 
   *baz=dbaz/rad;
/*
 * c
 * c  Bullen, Sec 10.2, eqn 7 / eqn 8
 * c
 * c    pt. 2 is unprimed, so this is technically the az
 * c
 * */
   rhs1=(a-dd)*(a-dd)+(b-ee)*(b-ee)+c*c - 2.;
   rhs2=(a-gg)*(a-gg)+(b-hh)*(b-hh)+(c-kk)*(c-kk) - 2.;
   daz=atan2(rhs1,rhs2);
   if(daz<0.0) {
      daz=daz+2*pi;
   }
   *az=daz/rad;
/*
 * c
 * c   Make sure 0.0 is always 0.0, not 360.
 * c
 * */
   if(abs(*baz-360.) < .00001) *baz=0.0;
   if(abs(*az-360.) < .00001) *az=0.0;
//   fprintf(stderr,"stla,stlo,evla,evlo: %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n",stla,stlo,evla,evlo,gcarc,az);
}

/** construct a cosine source time function(area of 1*dt) **/
float *cosine(float dura, float dt, int *ns) {
   int i;
   float amp, *src;
   * ns = rint(dura/dt);
   src = malloc((1+*ns)*sizeof(float));
   if (src == NULL) {*ns=0; return src;}
   for(i=0;i<*ns;i++) src[i] = (1.0 - cos(2*PI*i*dt/dura))/(*ns);
   return src;
}

/** 3 functions for sampling **/
double AverageRandom(double min,double max) {  
       int minInteger = (int)(min*1000);  
       int maxInteger = (int)(max*1000);  
       int randInteger = rand();  
       int diffInteger = maxInteger - minInteger;  
       int resultInteger = randInteger % diffInteger + minInteger; 
       return resultInteger/1000.0;  
}  

double Normal(double x,double miu,double sigma) {
       return 1.0/sqrt(2*PI)/sigma * exp(-1*(x-miu)*(x-miu)/(2*sigma*sigma));
}

double Normalnew(double x,double miu,double sigma) {
       return exp(-1*(x-miu)*(x-miu)/(2*sigma*sigma));
}

double NormalRandom(double miu,double sigma,double min,double max) {
       double x;
       double dScope;
       double y;
       do {
           x = AverageRandom(min,max);
           y = Normal(x, miu, sigma);
           dScope = AverageRandom(0,Normal(miu,miu,sigma));
          }while( dScope > y);
       return x;
}
