/***********************************************************************************************

                              Multiple Point Source Inversion Code
                            using Markov Chain Monte Carlo Optimizer

  Modify history:
   2016/ 8/ 5 Shengji Wei    -Cut-And-Paste (CAP) forward calculation code
   2017/ 7/10 Qibin   Shi    -change coordnate to ENZ(6) for regional waveform, but 
                              keep rtz(5) for teleseismic waveform using N->r & E->t  
                              (Structure DATA is modifed)
   2017/ 7/15 Qibin   Shi    -save data in the memory; untie shifts of 6 components
   2017/11/30 Qibin   Shi    -use station-dependent filtering
   2017/12/ 5 Qibin   Shi    -allow different length of GF, maximum npts is MAXNPT
   2018/ 4/23 Qibin   Shi    -read weight file only once
   2018/ 5/ 9 Qibin   Shi    -Joint inversion on regional and teleseismic waveform,
			                  assuming unsensitive sources' finiteness for the latter
   2018/ 7/24 Qibin   Shi    -calibrate time in weight file using smaller events
   2018/12/17 Qibin   Shi    -remove the option of using GF tp/s. Only use manual picks
   2018/12/19 Qibin   Shi    -OpenMP version: multiple threads from the 2nd Markov unit
   2019/ 2/18 Qibin   Shi    -add variance reduction
   2019/ 3/ 3 Qibin   Shi    -set boundary for parameter depth
   2019/ 4/22 Qibin   Shi    -if teleseismic, weight * 50, for balancing with reginal
                              waveform in joint inversion

   Potential memory corruption: several grn is allocated with memory in parallel in 
                             function read_sac().

  Reference:
   Q Shi, S Wei, M Chen. 2018. GJI.
   Q Shi, S Wei. 2020. GRL.

***********************************************************************************************/

#include <omp.h>
#include <math.h> 
#include <stdio.h>
#include <float.h> 
#include <stdlib.h>
#include <string.h>
#include "Complex.h"
#include "nrutil.h"
#include "cap.h"
#include "sac.h"
#include "nr.h"

#define  MAXDUR      2000
#define  MAXNPT      35000
#define  CHLEN       1000000
#define  NSRC        15
#define  NM          280

int      kd[6]={0,1,2,0,1,2};      // SE, SN, SZ, PE, PN, PZ
int      kk[6]={2,1,0,2,1,0};      // index of synthetic components
char     cm[3]={'E','N','Z'};
double   srcpara[NSRC][CHLEN][11];
double   ratio[CHLEN], L[CHLEN], VR[CHLEN];

int main (int argc, char **argv) {
  int    i,j,k,k1,l,m,nda,plot,useDisp,n_shft,total_n,m_pnl,m_sw,wndw_size,ns_pnl,ns_sw,n_data;
  int    distkm[STN],npt[STN],t0[STN],te[STN],t1[STN],t2[STN],t3[STN],t4[STN],durlen,io_taper;
  int    n1[STN],n2[STN],nft[6],mm[6],t[6][STN],n[6][STN],max_shft[6],shft0[6][STN],kmon[STN];
  long int order=4, nsects[STN];
  char	 tmp[NM],tmp1[NM],path[NM],*c_pt[STN],*c_pt1[STN],type[2]={'B','P'},proto[2]={'B','U'};
  float  pnl_sn[STN][30],pnl_sd[STN][30],sw_sn[STN][30],sw_sd[STN][30],wd1,wd2,sf1,sf2,riseTime;
  float  *distance,dmin=100.,vp,vs1,vs2,depSqr=110.25,con_shft[STN],w_cmp[6],rms_cut[5],weight;
  float  *crlss1[STN][6],*crlss2[STN][6],*crlss3[STN][6],*f_pt,*cc_pt,on_pt,nof_per_samp,tau0;
  float  x[STN],x1[STN],x2,xv,y[STN],xx[STN],yy[STN],y1,amp,dt,rad[6],dm;
  float  *src[NSRC],*src1[STN],*co[STN],*co1[STN],*syn,*data[STN][3],bs_body,bs_surf,bs[6];
  double arad[STN][3][3],m0[STN],coef,f1_pnl,f2_pnl,f1_sw,f2_sw,pie=3.141592654;
  FILE   *f_out,*ff_out;
  SACHEAD *hd0,*hd1,*hd2,*hd3;
  DATA   *obs,*obs0;
  COMP   *spt;
  /* MPS */
  int    npsrc,npara,n_shift,n_corr[STN],ii,z0,sol_shft[STN][6],no_proc=1;
  char   tmp_char1[200],tmp_char2[200];
  double Strike[NSRC],Dip[NSRC],Rake[NSRC],rise[NSRC],slip[NSRC],t_rupture[NSRC];
  double lat[NSRC],lon[NSRC],dep_src[NSRC],mu[NSRC],alp[NSRC],bet[NSRC];
  float  tp_ref[STN],ts_ref[STN],b_ref[STN],e_ref[STN],dist_ref[STN],tp[STN],ts[STN],*grn[STN];
  float  *syn_mps[STN][3],*syn_tmp2[STN][3],*syn_tmp3[STN][3],t_shift=200.,cfg[STN][6];
  float  *pt,*pt2,*ppv,*ppr,*ppt,*ppz,*ppe,*ppn,ddd;
  double total_m0,mw,total_error;
  double lo1,lo2,la1,la2,dp1,dp2,slop0,slop1,az[STN],gcarc[STN],baz[STN],stlat,stlon;
  /* MCMC */
  int    nchain,jj,qq,pp,qp,count=0,no_thrds;
  double sigma1=15.,sigma=1.,sigma2,beta=0.998,inter;
  double Mw[10]={0},jump[99]={0},total_obs=0.,temp,expo,scale[11];

  /* subroutines */
  void   distaz(double,double,double,double,double *,double *,double *);
  double NormalRandom(double,double,double,double);
  double Normalnew(double,double,double);
  double AverageRandom(double,double);
  double Normal(double,double,double);
  float  *cosine(float, float, int *);
  float  *yoffe(float,double,double,float,float *);
  float  *yoffecos(float,double,double,float,float *);

  /* input control parameters */
  for(i=0;i<11;i++) scanf("%lf",scale+i);
  scanf("%f%f%f%f%f",&bs_body,&bs_surf,&x2,&xv,&nof_per_samp);
  scanf("%lf%lf%lf%lf%lf%lf%lf%lf",&lo1,&la1,&lo2,&la2,&dp1,&dp2,&slop0,&slop1);
  scanf("%f%f%f%f",&vp,&vs1,&vs2,&dt);
  scanf("%f%f%f%f%lf%lf",&wd1,&wd2,&sf1,&sf2,&temp,&expo);
  scanf("%d%d%d%d%d",&nchain,&plot,&io_taper,&useDisp,&no_proc);
  scanf("%s",path);

  max_shft[3]=max_shft[4]=max_shft[5]=2*rint(sf1/dt); // max Pnl shift
  max_shft[0]=max_shft[1]=max_shft[2]=2*rint(sf2/dt); // max sur shift
  w_cmp[3]=w_cmp[4]=w_cmp[5]=x2; // relative weight of Pnl
  w_cmp[0]=w_cmp[1]=w_cmp[2]=1.; // relative weight of sur
  w_cmp[5]=w_cmp[5]*xv; //vertical weighting
  w_cmp[2]=w_cmp[2]*xv;
  bs[3]=bs[4]=bs[5]=bs_body; // distance scale power Pnl
  bs[0]=bs[1]=bs[2]=bs_surf; // distance scale power sur
  mm[3]=mm[4]=mm[5]=rint(wd1/dt); // max Pnl window
  mm[0]=mm[1]=mm[2]=rint(wd2/dt); // max sur window
  nft[3]=2; while(nft[3]<mm[3]) nft[3]*=2;
  nft[0]=2; while(nft[0]<mm[0]) nft[0]*=2;
  nft[2]=nft[1]=nft[0];
  nft[5]=nft[4]=nft[3];
  n_shift = rint(t_shift/dt);

  /* initial source parameters */
  scanf("%d%f",&npsrc,&ddd);
  npara=npsrc*11-1;
  for(ii=0;ii<npsrc;ii++){
    scanf("%le",lat+ii);
    scanf("%le",lon+ii);
    scanf("%le",dep_src+ii);
    scanf("%le",alp+ii);
    scanf("%le",bet+ii);
    scanf("%le",t_rupture+ii);
    scanf("%le",slip+ii);
    scanf("%le",Strike+ii);
    scanf("%le",Dip+ii);
    scanf("%le",Rake+ii);
    scanf("%le",rise+ii);
    scanf("%le",mu+ii);
    srcpara[ii][0][0] = lat[ii];
    srcpara[ii][0][1] = lon[ii];
    srcpara[ii][0][2] = dep_src[ii];
    srcpara[ii][0][3] = t_rupture[ii];
    srcpara[ii][0][4] = log10(slip[ii]);
    srcpara[ii][0][5] = Strike[ii];
    srcpara[ii][0][6] = Dip[ii];
    srcpara[ii][0][7] = Rake[ii];
    srcpara[ii][0][8] = rise[ii];
    srcpara[ii][0][9] = 1.-alp[ii];
    srcpara[ii][0][10] = 1.-bet[ii];
  }

  scanf("%d",&nda);
  if(nda    > STN  ){fprintf(stderr,"!Max. %d stations!  \n",STN  );nda   =STN  ;}
  if(nchain > CHLEN){fprintf(stderr,"!Max. %d iterations!\n",CHLEN);nchain=CHLEN;}

  //beta=pow(temp,(float)npara/(float)nchain); // exponential cooling rate
  //fprintf(stderr,"Cooling rate =  %f\n",beta);

  /* allocate memory */
  obs0 = (DATA *) malloc(nda*sizeof(DATA));
  if(obs0 == NULL){fprintf(stderr,"No memory for the whole thing\n");return -1;}

  hd0  = (SACHEAD *) malloc(nda*sizeof(SACHEAD));
  if(hd0  == NULL){fprintf(stderr,"No memory for data head info \n");return -1;}

  hd1  = (SACHEAD *) malloc(nda*sizeof(SACHEAD));
  if(hd1  == NULL){fprintf(stderr,"No memory for synthetics head\n");return -1;}

  distance = (float *) malloc(nda*sizeof(float));
  for(i=0;i<NSRC;i++) {src[i] = (float *) malloc(MAXDUR*sizeof(float));}
  for(i=0;i<STN;i++) {co[i]  = (float *) malloc((MAXDUR+MAXNPT)*sizeof(float));}
  for(i=0;i<STN;i++) {for(k=0; k<3; k++) syn_mps[i][k] =(float *) calloc(MAXNPT+2*n_shift,sizeof(float));}
  for(i=0;i<STN;i++) {for(k=0; k<3; k++) syn_tmp2[i][k]=(float *) calloc(MAXNPT,sizeof(float));}
  for(i=0;i<STN;i++) {for(k=0; k<3; k++) syn_tmp3[i][k]=(float *) calloc(MAXNPT,sizeof(float));}
  for(i=0;i<STN;i++) {for(k=0;k<NCP;k++) crlss1[i][k]=(float *) malloc(2*nft[k]*sizeof(float));}
  for(i=0;i<STN;i++) {for(k=0;k<NCP;k++) crlss2[i][k]=(float *) malloc(2*nft[k]*sizeof(float));}
  for(i=0;i<STN;i++) {for(k=0;k<NCP;k++) crlss3[i][k]=(float *) malloc((max_shft[k]+1)*sizeof(float));}

  /* read sta-info and data */
  for(i=0;i<nda;i++) {

    obs = obs0+i;
    hd2 = hd0 +i; // data head

    scanf("%s",obs->stn); // station name
    scanf("%d",kmon+i); // distance unit switch for grn
    for(j=0;j<6;j++) scanf("%d",&obs->com[5-j].on_off); // weight
    for(j=0;j<6;j++) scanf("%d",shft0[5-j]+i); // calibrating time shift
    scanf("%lf%lf%lf%lf",&f1_pnl,&f2_pnl,&f1_sw,&f2_sw);// corner freq
    design(order,type,proto,1.,1.,f1_pnl,f2_pnl,(double)dt,pnl_sn[i],pnl_sd[i],&nsects[i]);
    design(order,type,proto,1.,1.,f1_sw, f2_sw, (double)dt, sw_sn[i], sw_sd[i],&nsects[i]);

    strcat(strcat(strcat(strcpy(tmp,argv[1]),"/"),obs->stn),".HHE");
    c_pt[i] = strrchr(tmp,(int) 'E');
    for(n_data=0;n_data<NRC;n_data++){
      *(c_pt[i]) = cm[n_data];
      if((data[i][n_data] = read_sac(tmp,hd2)) == NULL){
        fprintf(stderr,"Skip reading %s !\n",tmp);
        goto next_sta;
      }
    }
    obs->dist = distance[i] = hd2->dist;
    obs->az   = hd2->az;
    obs->stla = hd2->stla;
    obs->stlo = hd2->stlo;

    goto normal;
    next_sta:
      nda--; i--;
    normal:
      fprintf(stderr,"Read Record No. %3d %s\n",i+1,obs->stn);
  }
  fprintf(stderr,"See you later !\n");

  /* Markov Chain Monte Carlo */
  ratio[0]=0;
  ff_out=fopen("pdf.out","w");

  for(jj=0;jj<nchain;jj++) {

    for(ii=0;ii<npsrc;ii++){                // set bounds of model space

      if     (srcpara[ii][jj][0]>la1+slop0*(srcpara[ii][jj][1]-lo1))  
	      srcpara[ii][jj][0]=la1+slop0*(srcpara[ii][jj][1]-lo1);
      else if(srcpara[ii][jj][0]<la2+slop0*(srcpara[ii][jj][1]-lo2))  
	      srcpara[ii][jj][0]=la2+slop0*(srcpara[ii][jj][1]-lo2);
      if     (srcpara[ii][jj][1]>lo2+slop1*(srcpara[ii][jj][0]-la2))  
	      srcpara[ii][jj][1]=lo2+slop1*(srcpara[ii][jj][0]-la2);
      else if(srcpara[ii][jj][1]<lo1+slop1*(srcpara[ii][jj][0]-la1))  
	      srcpara[ii][jj][1]=lo1+slop1*(srcpara[ii][jj][0]-la1);
      if     (srcpara[ii][jj][2]>dp2)  
	      srcpara[ii][jj][2]=dp2;
      else if(srcpara[ii][jj][2]<dp1)  
	      srcpara[ii][jj][2]=dp1;

      if     (srcpara[ii][jj][3]<0.)  
	      srcpara[ii][jj][3]=0.;
      if     (srcpara[ii][jj][4]>1.5)  
	      srcpara[ii][jj][4]=1.5;
      else if(srcpara[ii][jj][4]<-3.)  
	      srcpara[ii][jj][4]=-3.;

      if     (srcpara[ii][jj][6]>90.) {
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
      if     (srcpara[ii][jj][7]>180. )  srcpara[ii][jj][7] -= 360.;
      else if(srcpara[ii][jj][7]<-180.)  srcpara[ii][jj][7] += 360.;
      if     (srcpara[ii][jj][8]<0.5  )  srcpara[ii][jj][8]  = 0.5 ; 
      else if(srcpara[ii][jj][8]>100. )  srcpara[ii][jj][8]  = 100.;
      if     (srcpara[ii][jj][9]<0.01 )  srcpara[ii][jj][9]  = 0.01;
      else if(srcpara[ii][jj][9]>1.   )  srcpara[ii][jj][9]  = 1.  ;
      if     (srcpara[ii][jj][10]<0.05)  srcpara[ii][jj][10]  = 0.05;
      else if(srcpara[ii][jj][10]>1.  )  srcpara[ii][jj][10]  = 1. ;

      lat[ii]       = srcpara[ii][jj][0];
      lon[ii]       = srcpara[ii][jj][1];
      dep_src[ii]   = srcpara[ii][jj][2];
      t_rupture[ii] = srcpara[ii][jj][3];
      slip[ii]      = pow(10.0, srcpara[ii][jj][4]);
      Strike[ii]    = srcpara[ii][jj][5];
      Dip[ii]       = srcpara[ii][jj][6];
      Rake[ii]      = srcpara[ii][jj][7];
      rise[ii]      = srcpara[ii][jj][8];
      alp[ii]       = srcpara[ii][jj][9];
      bet[ii]       = srcpara[ii][jj][10];
      
      for(k1=0;k1<MAXDUR;k1++) src[ii][k1]=0.;
    // trap(2*rise[ii], rise[ii], dt, src[ii]);
    //  yoffe(2*rise[ii], alp[ii], bet[ii], dt, src[ii]);
      yoffecos(2*rise[ii], alp[ii], bet[ii], dt, src[ii]);
    }

    total_error=0.; total_n=0; n_shft=0; total_obs=0.;
    
// OMP
if(jj==0) {no_thrds=1;} else {no_thrds=no_proc;}
omp_set_num_threads(no_thrds);
#pragma omp parallel for private (spt,obs,hd2,hd3,ii,j,l,k,k1,m,y1,tmp,tmp_char1,tmp_char2,coef,pt,pt2,durlen,ppv,ppr,ppt,ppz,ppn,ppe,weight,wndw_size,x2,f_pt,cc_pt,on_pt) reduction (+:total_error,total_obs,total_n,n_shft)

    /* parallel loop over stations */
    for(i=0;i<nda;i++) {
      obs=obs0+i;
      hd2=hd0 +i; // data head
      hd3=hd1 +i; //  grn head

      x[i]  = hd2->b-hd2->o;
      y[i]  = hd2->e-hd2->o;
      x1[i] = hd2->a-hd2->o;
      t1[i] = rint((hd2->t1-hd2->b)/dt);
      t2[i] = rint((hd2->t2-hd2->b)/dt);
      t3[i] = rint((hd2->t3-hd2->b)/dt);
      t4[i] = rint((hd2->t4-hd2->b)/dt);

      /* compute MPS synthetics */
      for(ii=0;ii<npsrc;ii++){
        distaz((double) obs->stla,(double) obs->stlo,lat[ii],lon[ii],gcarc+i,baz+i,az+i);
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
       // if(kmon[i]==0){
          distkm[i] = gcarc[i]*pie*6371/180;
          if(distkm[i]<1) distkm[i]=1;
          sprintf(tmp_char2,"%d",distkm[i]);
	    }else{
          sprintf(tmp_char2,"%4.1f",gcarc[i]);
	    }
        
        strcat(strcat(tmp_char1,tmp_char2),".grn.0");
        c_pt[i] = strrchr(tmp_char1,(int) '.') + 1;

        dc_radiat(az[i]-Strike[ii],Dip[ii],Rake[ii],arad[i]);
        m0[i]   = mu[ii] * slip[ii] * 1.0e-13;

        for(l=0; l<3; l++){ // nine components
            for(j=0; j<3; j++){

	        coef = m0[i] * arad[i][l][j];

	        if((grn[i]=read_sac(tmp_char1,hd3)) == NULL) continue;

            if(ii==0 && l==0 && j==0){ // 1st source's head

                     npt[i] = hd3->npts;
	               b_ref[i] = hd3->b;
                   e_ref[i] = hd3->e;
                  tp_ref[i] = hd3->t1;
                  ts_ref[i] = hd3->t2;
	            dist_ref[i] = hd3->dist;
   	             obs->alpha = hd3->user1;
                  n_corr[i] = rint(t_rupture[ii]/dt) + n_shift;

              for(k=0; k<3; k++) {
                for(k1=0,pt=syn_mps[i][k];k1<npt[i]+2*n_shift;k1++,pt++) *pt=0;
              }

            }else if(l==0 && j==0){ // following sources' head
              if(npt[i] < hd3->npts) npt[i] = hd3->npts; // npt is the largest
              tp[i] = hd3->t1;
              ts[i] = hd3->t2;
              n_corr[i] = rint((tp[i] - tp_ref[i] + t_rupture[ii])/dt) + n_shift;
            }

            for(pt=syn_tmp2[i][j],k1=0;k1<(hd3->npts);k1++,pt++) (*pt) += coef*grn[i][k1];
	        free(grn[i]);

            (*(c_pt[i]))++;

          }
        } // finish nine components

        durlen = rint(2*rise[ii]/dt);
        for(k=0;k<3;k++) conv(src[ii],durlen,syn_tmp2[i][k],npt[i],co[i]); // convolution
	    for(k1=0;k1<MAXDUR+MAXNPT;k1++) co[i][k1]=0;

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
//             (*pt2) += (*pt)*NormalRandom(1.,0.2,0.2,1.8);  (*pt) = 0; /*Gaussian Noise*/
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
     
      /*Output MPS synthetics for testing, remember to set t_shift=0 and Gaussian noise. Comment this block for inversion
      fprintf(stderr,"simulate sta---- %s\n",obs->stn);
      strcat(obs->stn,".z");
      hd2->b = b_ref[i];
      c_pt[i] = strrchr(obs->stn,(int) 'z');
      write_sac(obs->stn,*hd2,syn_mps[i][0]);
      *(c_pt[i]) = 'n';
      write_sac(obs->stn,*hd2,syn_mps[i][1]);
      *(c_pt[i]) = 'e';
      write_sac(obs->stn,*hd2,syn_mps[i][2]); */

      /* update header of synthetics */
      hd3->b    = b_ref[i]  - t_shift; // add more points before the beginning
      hd3->e    = e_ref[i]  + t_shift; // add more points after the end
      hd3->t1   = tp_ref[i] + t_rupture[0];
      hd3->t2   = ts_ref[i] + t_rupture[0];
      hd3->npts = npt[i]    + 2*n_shift;
      hd3->dist = dist_ref[i];

      /* use first-arrival (if picked) to correct synthetics */
      con_shft[i] = 0;
//      if(x1[i]>0.) con_shft[i] += x1[i] - hd3->t1;
      t0[i]=rint((x[i]-con_shft[i]-hd3->b)/dt) + 1;
      te[i]=rint((y[i]-con_shft[i]-hd3->b)/dt) - 1;

      /* Pnl window (relative to b of syn) */
      if(t1[i]>0 && t2[i]>0) {t1[i]+=t0[i]; t2[i]+=t0[i];} else{fprintf(stderr,"! No P window for %s",obs->stn);}
      if(t1[i]<shft0[5][i]) t1[i]=shft0[5][i];
      if(t1[i]<shft0[4][i]) t1[i]=shft0[4][i];
      if(t1[i]<shft0[3][i]) t1[i]=shft0[3][i];
      if(t2[i]>te[i]) t2[i]=te[i];
      if(t2[i]>hd3->npts+shft0[5][i]) t2[i]=hd3->npts+shft0[5][i];
      if(t2[i]>hd3->npts+shft0[4][i]) t2[i]=hd3->npts+shft0[4][i];
      if(t2[i]>hd3->npts+shft0[3][i]) t2[i]=hd3->npts+shft0[3][i];

      /* Sur window (relative to b of syn) */
      if(t3[i]>0 && t4[i]>0) {t3[i]+=t0[i]; t4[i]+=t0[i];} else{fprintf(stderr,"! No S window for %s",obs->stn);}
      if(t3[i]<t1[i]) t3[i]=t1[i];
      if(t3[i]<t1[i]-shft0[5][i]+shft0[2][i]) t3[i]=t1[i]-shft0[5][i]+shft0[2][i];
      if(t3[i]<t1[i]-shft0[4][i]+shft0[1][i]) t3[i]=t1[i]-shft0[4][i]+shft0[1][i];
      if(t3[i]<t1[i]-shft0[3][i]+shft0[0][i]) t3[i]=t1[i]-shft0[3][i]+shft0[0][i];
      if(t4[i]>te[i]) t4[i]=te[i];
      if(t4[i]>hd3->npts+shft0[2][i]) t4[i]=hd3->npts+shft0[2][i];
      if(t4[i]>hd3->npts+shft0[1][i]) t4[i]=hd3->npts+shft0[1][i];
      if(t4[i]>hd3->npts+shft0[0][i]) t4[i]=hd3->npts+shft0[0][i];

      /* windows length **/
      n1[i] = t2[i] - t1[i];
      n2[i] = t4[i] - t3[i];
      if(n1[i]>mm[3]) n1[i]=mm[3];
      if(n2[i]>mm[0]) n2[i]=mm[0];

      /* discard too short windows */
      if(n1[i]<max_shft[3]){
        fprintf(ff_out,"Shift less for P (Pnl) of %s !\n", obs->stn);
        max_shft[3]=n1[i];
      }
      if(n2[i]<max_shft[0]){
        fprintf(ff_out,"Shift less for S (sur) of %s !\n", obs->stn);
        max_shft[0]=n2[i];
      }

      /* cut */
      t[0][i]=t[1][i]=t[2][i]=t3[i]; // start of Sur components   
      t[3][i]=t[4][i]=t[5][i]=t1[i]; // start of Pnl components
      n[0][i]=n[1][i]=n[2][i]=n2[i]; // length of Sur components
      n[3][i]=n[4][i]=n[5][i]=n1[i]; // length of Pnl components

      if(obs->com[0].on_off>0 || obs->com[1].on_off>0 || obs->com[2].on_off>0) n_shft++;
      if(obs->com[3].on_off>0 || obs->com[4].on_off>0 || obs->com[5].on_off>0) n_shft++;

      for(spt=obs->com,j=0;j<NCP;j++,spt++){
	    if(spt->on_off) total_n+=npt[i];
        spt->npt  = npt[i] = n[j][i];
        spt->b    = t[j][i]*dt+con_shft[i]+hd3->b; 
        weight    = w_cmp[j]*pow(distance[i]/dmin,bs[j]);
	    if(kmon[i]==0){
          weight *=4.; // weight for teleseismic in joint inversion
	      if(j<3) {weight *=1.;} // extra weight for tele S
	    }else{
          if(j>2) {weight *=10.;} // extra weight for local P
        }
        wndw_size = npt[i]*sizeof(float);

	/* data */
        if(jj==0){

          spt->rec = (float *) malloc(wndw_size);
          spt->syn = (float *) malloc(wndw_size);
          memcpy(spt->rec, data[i][kd[j]] + t[j][i] - t0[i], wndw_size);
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
        memcpy(f_pt, syn_mps[i][kk[j]] + t[j][i] - shft0[j][i], wndw_size);
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
        spt->crl=crscrl(npt[i],nft[j],max_shft[j],spt->rec,f_pt,crlss1[i][j],crlss2[i][j],crlss3[i][j]);

      }

      /* find time shifts */
      for(spt=obs->com,j=0;j<NCP;j++,spt++){
        cc_pt=spt->crl;
        on_pt=spt->on_off>0?1:0;
	    for(y1=-FLT_MAX,l=0;l<=max_shft[j];l++){
          xx[i] = *cc_pt++;
          yy[i] = on_pt*xx[i];
          if(yy[i]>y1){
            y1     = yy[i];
	        cfg[i][j] = xx[i];
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
//        spt->shft    = dt * (shft0[j][i] + sol_shft[i][j]);
	    spt->shft    = dt * sol_shft[i][j];
        spt->ampt    = sqrt(spt->syn2 / spt->rec2);  
        total_error += spt->err;
        total_obs   += spt->on_off * spt->rec2;	
      }

      /* output waveforms */
      if(jj == nchain-1){
        strcat(strcat(strcat(strcat(strcat(strcpy(tmp,argv[1]),"/"),"model_mps"), "_"),obs->stn),".11");
        c_pt[i] = strrchr(tmp,(int) '1');
        c_pt1[i]= c_pt[i]-1;
        for(spt=obs->com,j=0;j<NCP;j++,spt++){
          npt[i]    = spt->npt;
          (*hd3)    = sachdr(dt, npt[i], spt->b); // add b, dt, npt into head
          hd3->dist = obs->dist;
          hd3->az   = obs->az;
          hd3->user1= obs->alpha;
          hd3->a    = hd3->b;
          write_sac(tmp,(*hd3),spt->rec);
          (*(c_pt1[i]))++;
          hd3->b   -= (shft0[j][i]*dt + con_shft[i]);
          hd3->a    = hd3->b - sol_shft[i][j]*dt;
          write_sac(tmp,(*hd3),spt->syn);
          (*(c_pt[i]))++;
          (*(c_pt1[i]))--;
        }
      }

      if(jj==0){
        for(j=0;j<NRC;j++) {free(data[i][j]);}
      }

    } // end parallel loops over stations

    /* Evaluate the quality of waveform fit */
//    L[jj]  = Normalnew(total_error,0.,sigma1);
    VR[jj] = 1 - (total_error/total_obs);
    sigma2   = 0.5*(1-temp)*cos(pie*jj/nchain) + 0.5*(1+temp);
    if(jj > 0 ) {

//      ratio[jj] =  sqrt(jj) * log(jj) * (VR[jj]/VR[jj-1] - 1) + 1.;
//      ratio[jj] =  jj/2 * (VR[jj]/VR[jj-1] - 1) + 1.;
//      ratio[jj] = (double) jj / (double) sigma2 * (VR[jj]-VR[jj-1]) + 1.;
//      ratio[jj] = (double) jj * (double) sigma2 * (VR[jj]-VR[jj-1]) + 1.;
      ratio[jj] = pow( ((double)jj+3000.) , expo ) * (double)sigma2 * (VR[jj]-VR[jj-1]) + 1.;
//      ratio[jj] = (double) jj * (VR[jj]-VR[jj-1]) + 1.;
//        ratio[jj] = sqrt((double)jj) / (double) sigma2 / (double) sigma2 * (VR[jj]-VR[jj-1]) + 1.;
        if(ratio[jj] < AverageRandom(0.,1.) ){
//          L[jj] = L[jj-1];
	      VR[jj] = VR[jj-1];
	      for(ii=0;ii<npsrc;ii++) {
            for(qq=0;qq<11;qq++) srcpara[ii][jj][qq] = srcpara[ii][jj-1][qq];
          }
        }        
    }
//    if(jj == 0 && VR[jj] < 0.0001) VR[jj] = 0.0001; 

    /* next MCMC sample */
    if(jj<nchain-1){
      if(ratio[jj] <= 1.){
        for(qq=0;qq<npara;qq++) jump[qq] = 0.;
        qp = count%npara; count++;
        jump[qp] = NormalRandom(0.,sigma2,(0-4*sigma2),4*sigma2);
        if(jump[qp] < (0.05*sigma2) && jump[qp] >= 0.0) jump[qp] = 0.05*sigma2;
        else if(jump[qp] > (0.-0.05*sigma2) && jump[qp] <= 0.0) jump[qp] = 0.-0.05*sigma2;
      }
//      sigma2   = 0.5*(1-temp)*cos(pie*jj/nchain) + 0.5*(1+temp);
//      sigma2   = sigma * pow(beta, jj/npara);
//      inter    = 6 * cos(pie*jj/(nchain*2.5));
//      jump[qp] = NormalRandom(0.,sigma2,(0-inter),inter);

      for(qq=0;qq<3;qq++)
	    srcpara[0][jj+1][qq] = jump[qq]  *scale[qq]+srcpara[0][jj][qq];
      srcpara[0][jj+1][3] = srcpara[0][jj][3];
      for(qq=4;qq<11;qq++) 
	    srcpara[0][jj+1][qq] = jump[qq-1]*scale[qq]+srcpara[0][jj][qq];
      for(ii=1;ii<npsrc;ii++){
        for(qq=0;qq<11;qq++)
	      srcpara[ii][jj+1][qq] = jump[qq+ii*11-1]*scale[qq]+srcpara[ii][jj][qq]; 
      }
    }

    for(ii=0;ii<npsrc;ii++){
      fprintf(ff_out,"%10.6lf %10.6lf %10.2lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf  ",
           srcpara[ii][jj][0],
		   srcpara[ii][jj][1],
		   srcpara[ii][jj][2],
		   srcpara[ii][jj][3],
           pow(10.0, srcpara[ii][jj][4]),
		   srcpara[ii][jj][5],
		   srcpara[ii][jj][6],
		   srcpara[ii][jj][7],
		   srcpara[ii][jj][8],
           srcpara[ii][jj][9],
           srcpara[ii][jj][10]);
    }
    fprintf(ff_out,"%.5lf %d %.5lf\n",ratio[jj],count,VR[jj]);
  } // end of MCMC
  fclose(ff_out);


  /* final output */
  total_m0 = 0.;
  for(ii=0;ii<npsrc;ii++) {
    total_m0 += mu[ii]*slip[ii];
    Mw[ii] = (log10(mu[ii]*slip[ii]) - 9.1)*2./3.;
  }
  mw = (log10(total_m0) - 9.1)*2./3.;
  fprintf(stderr,"Mw %8.2lf\n",mw);

  strcat(strcat(strcat(strcpy(tmp1,argv[1]),"/"),"model_mps"),".out");
  f_out=fopen(tmp1,"w");
  fprintf(f_out,"Event %s Model MPS ",argv[1]);
  for(ii=0;ii<npsrc;ii++)
    fprintf(f_out,"FM %3.0lf %3.0lf %3.0lf Dp %2.0lf Mw %5.2lf Delay %4.1lf dura %4.1lf ",
		    Strike[ii],
		    Dip[ii],
		    Rake[ii],
		    dep_src[ii],
		    Mw[ii],
		    t_rupture[ii],
		    2*rise[ii]);
  fprintf(f_out,"rms %9.2le loc ",total_error);
  for(ii=0;ii<npsrc;ii++) fprintf(f_out,"%8.5lf %8.4lf ",lat[ii],lon[ii]);
  fprintf(f_out,"\n");
  for(i=0,obs=obs0;i<nda;i++,obs++) {
    fprintf(f_out,"%9s %5.2f ",obs->stn,con_shft[i]);
    for(j=5,spt=obs->com+5;j>=0;j--,spt--){
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
} // The end of the whole thing

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

float *trap(float t1, float t2, float dt, float *s) {
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
    for(i=1;i<=n1;i++) s[n1+n2-i] = s[i]=s[i-1] + slope;
    for(;i<n2;i++) s[i]=s[i-1];
    s[n1+n2]=0;

    return s;
}

float *yoffe(float tr, double alp, double bet, float dt, float *s) {
    int i, nr, n1, n2;
    float slope1, slope2;

    nr = rint(tr/dt); if (nr<5) nr = 5;
    n1 = rint(bet*tr/2./dt); if (n1<2) n1 = 2;
    n2 = rint((1.-alp/2.)*bet*tr/dt); if (n2<3) n2 = 3;

    slope1 = 4./(nr*nr*bet*(alp+bet-alp*bet));
    slope2 = 2.*alp/(nr*nr*(alp+bet-alp*bet)*(1.-bet+alp*bet/2.));

    s[0] = 0;
    for(i=1;i<=n1;i++) s[2*n1-i] = s[i] = s[i-1] + slope1;
    for(i=n2+1;i<nr;i++) s[i] = s[i-1] - slope2; 
    s[nr] = 0;

    return s;
}

float *yoffecos(float tr, double alp, double bet, float dt, float *s) {
    int i, nr, n1, n2;
    float am;
    FILE *stfout;
    
    alp = alp/(1.01-alp);
    nr = rint(tr/dt); if (nr<5) nr = 5;
    n1 = rint(bet*tr/2./dt); if (n1<2) n1 = 2;
    n2 = nr - n1;
    am = 1./nr/(alp + bet + alp*bet*(4.0-PI)/PI/2.0);

    for(i=0;i<n1;i++) s[i] = am*(1.-cos(PI*i/n1)) + 2.*alp*am*sin(PI*i/n1/2.);
    for(i=n1;i<2*n1;i++) s[i] = am*(1.-cos(PI*i/n1)) + alp*am*(1.+cos(PI*(i-n1)/n2));
    for(i=2*n1;i<nr;i++) s[i] = alp*am*(1.+cos(PI*(i-n1)/n2));
    s[nr] = 0;
    
    //stfout = fopen("stf.out","a");
    //fprintf(stfout,">\n");
    //for(i=0;i<=nr;i++) fprintf(stfout,"%3d %f\n",i,s[i]);
    //fclose(stfout);
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

void distaz(double stla,double stlo,double evla,double evlo,double *gcarc,double *baz,double *az){
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
  double scolat, slon, ecolat, elon;
  double a,b,c,d,e,aa,bb,cc,dd,ee,g,gg,h,hh,k,kk;
  double rhs1,rhs2,sph,rad,del,daz,dbaz,pi,piby2;
  
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
       return exp((miu-x)/sigma);
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
