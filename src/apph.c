//=============================================================================
// in cartesian coordinates t,x,y,z for x in [0,1], y in [0,1], z in [-1,1]
//
// NOTES: ==> uses the global variables of CFE4D
//        ==> hierarchy between Lmin and Lmax must be in-sync
//        ==> assumed time structure ... n,nm1,np1
//        ==> the current AH structure is specified by c_AH (0 indexed)
//
// The following routine searches for an apparent horizon over levels
// AH_Lmin to AH_Lmax.
//
// The hypersurface describing the AH is given by
//
// r = AH_R(chi,phi)
//
// An iterative flow-method is used, to absolute tolerance AH_tol, maximum
// iteration AH_max_iter. 
//
// If AH_tol<0 then the tolerance is dynamically adjusted (via change in R of min(dx)*-AH_tol)
// 
// If (use_R_ic), then existing value of R is used as an initial guess, else
// AH_rsteps values of R=AH_r0 to R=AH_r1 are used, with the one giving the
// smallest expansion used to start.
//
// AH_theta(AH_Nchi,AH_Nphi) is a real work array used to hold the null expansion.
// AH_own(AH_Nchi,AH_Nphi) and AH_lev(AH_Nchi,AH_Nphi) are integer work arrays describing the
// ownership of points.
// 
// returns (true) if found, and if so M & J are filled in to reflect the
// effect mass and angular momentum via the dynamical horizons approach i
// (not implemented yet)
//=============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <pamr.h>
#include <amrd.h>
#include <math.h>
#include <m_util_r8.h>
#include <bbhutil.h>
#include <mpi.h>
#include "CFE4D.h"
#include "apph.h"

#define UNOWNED -1e10
real AH_ct[MAX_BHS];

//=============================================================================
// some utility routines called by find_apph below
//=============================================================================

//-----------------------------------------------------------------------------
// initialize ownership array ... i.e., those points this node is
// responsible for
//-----------------------------------------------------------------------------
int fill_own(int Lmax, int ltrace, int *first)
{
   int i,np,L,valid;

   np=AH_Nchi[c_AH]*AH_Nphi[c_AH];

   for (i=0; i<np; i++) AH_own[c_AH][i]=-1;

   for (L=Lmax; L>=AH_Lmin[c_AH]; L--)
   {
      valid=PAMR_init_s_iter(L,PAMR_AMRH,1);  // loop over ALL grids!

      while(valid)
      {
         ldptr();

         ah_fill_own_(AH_R[c_AH],AH_xc[c_AH],AH_own[c_AH],AH_lev[c_AH],bbox,&dx,&dy,&g_rank,&L,&AH_Nchi[c_AH],&AH_Nphi[c_AH], &axisym);
         valid=PAMR_next_g();
      }
   }
   *first=0;

   for (i=0; i<np ; i++) if (AH_own[c_AH][i]==-1) 
   {
      if (my_rank==0 && ltrace)
      printf("find_apph: error ... point %i,%i is not owned; at AH_R=%f\n,chi=%f\n",i%AH_Nchi[c_AH]+1,i/AH_Nchi[c_AH]+1,AH_R[c_AH][i],(i%AH_Nchi[c_AH])*M_PI/(AH_Nchi[c_AH]-1));
      return 0;        
   }

   return 1;
}

//-----------------------------------------------------------------------------
// computes theta, assuming ownership calculated ... alters AH_w1
// is_ex set to 1 if any point couldn't be calculated due to closeness
// of excision zone ... returns theta average
//-----------------------------------------------------------------------------
#define USE_SMOOTH_A 1
#define MAX_TRACE 5000
real fill_theta(double *AH_theta0, real eps0, real *area, real *c_equat, real *c_polar, int *is_ex)
{
   int i,np,valid,dvtrace=0,i0,j0,is_int;
   static int num_trace=0;
   char name[256];
   int AH_shape[2],rank;
   real AH_bbox[4],resid,da[3],area_owned[3],area_global[3];  // elements 2 & 3 for circumferences.

   int k;
   real tmp=1.0;

   // outputs AH_*_iter gfns
   dvtrace=1;

   np=AH_Nchi[c_AH]*AH_Nphi[c_AH];

   area_owned[0]=area_owned[1]=area_owned[2]=0;
   *is_ex=0;

   for (i=0; i<np; i++) 
   {
      AH_theta0[i]=UNOWNED;
      if (AH_own[c_AH][i]==my_rank)
      {
         i0=i%AH_Nchi[c_AH]+1;
         j0=i/AH_Nchi[c_AH]+1;
         valid=PAMR_init_s_iter(AH_lev[c_AH][i],PAMR_AMRH,0); 
         while(valid)
         {
            ldptr_bbox();
            ah_is_int_(&is_int,AH_R[c_AH],AH_xc[c_AH],&i0,&j0,bbox,&dx,&dy,&AH_Nchi[c_AH],&AH_Nphi[c_AH],&axisym);
            if (is_int)
            {
               ldptr(); 

               // compute full theta 
               //(NOTE: in _post_tstep, have cycled time sequence np1,n,nm1 to time sequence n,nm1,np1,
               // so here, time level n is the most advanced time level)
               calc_exp0_(AH_R[c_AH],AH_xc[c_AH],AH_theta0,&i0,&j0,&AH_Nchi[c_AH],&AH_Nphi[c_AH],
                       theta,f,&da[0],&da[1],&da[2],
                       gb_tt_n,gb_tt_nm1,gb_tt_np1,
                       gb_tx_n,gb_tx_nm1,gb_tx_np1,
                       gb_ty_n,gb_ty_nm1,gb_ty_np1,
                       gb_xx_n,gb_xx_nm1,gb_xx_np1,
                       gb_xy_n,gb_xy_nm1,gb_xy_np1,
                       gb_yy_n,gb_yy_nm1,gb_yy_np1,
                       psi_n,psi_nm1,psi_np1,
                       &AdS_L,x,y,z,&dt,chr,&AMRD_ex,&AMRD_do_ex,&Nx,&Ny,&Nz,&axisym);

               area_owned[0]+=da[0];
               area_owned[1]+=da[1];
               area_owned[2]+=da[2];
               valid=0;
            }
            else valid=PAMR_next_g();
         }
         // to quit "cleanly"
         if (AH_theta0[i]==UNOWNED) area_owned[1]=1e10;

//         //TEST7//
//         area_owned[2]=0;
//         if (i0!=1 && i0!=AH_Nchi[c_AH]) 
//         {
//           tmp=0.0;
//           for (k=0; k<np; k++) tmp+=AH_theta0[k];
//           if (tmp==-7*np) {printf("\nBREAKING\n"); for (k=0; k<np; k++) AH_theta0[k]=0; area_owned[2]=1; break;}
//         } 
//         //TEST7//

      }

   }

   // for each i point on the AH surface, save max{AH_theta0[i]}_allprocessors into AH_w1[i],
   // and sum{area_owned[i]} into area_global[i]
   MPI_Allreduce(AH_theta0,AH_w1[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
   MPI_Allreduce(area_owned,area_global,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

   *area=area_global[0];
   *c_equat=area_global[1];
   if (*c_equat>1e9) AMRD_stop("fill_theta: error ... point 'unowned'",""); //quitting cleanly (cf. above)
   *c_polar=area_global[2];

   if (*area<0) { *is_ex=1; *area=0; }

//   for (i=0; i<np; i++) {if (AH_w1[c_AH][i]!=0) printf("AH_w1[c_AH][i]=%f, i=%i\n",AH_w1[c_AH][i],i);} //TMP//

   // copy AH_w1 into AH_theta0
   for (i=0; i<np; i++) {AH_theta0[i]=AH_w1[c_AH][i];}

   //--------------------------------------------------------------------------
   // regularity is essential on the axis, and we need to do it before
   // smoothing as calc_exp0 does not fill in the axis, and afterwards
   // again to make sure that theta and hence R is always exactly regular
   //--------------------------------------------------------------------------
   if (AH_xc[c_AH][1]==0) 
   { 
     reg_ah_r_(AH_theta0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
   }

   // AH smoothing
   if (eps0>0 && eps0<1)
   {
      if (USE_SMOOTH_A) 
      {
        smooth_ah_r_(AH_theta0,AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
      }
      else 
      {
        smooth_ah_r_b_(AH_theta0,AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
      }
   }
   else if (eps0>1 && eps0<2) 
   {
      eps0-=1;
      if (USE_SMOOTH_A) 
      {
        smooth_ah_r_(AH_theta0,AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
      }
      else 
      {
        smooth_ah_r_b_(AH_theta0,AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        smooth_ah_r_(AH_theta0,AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        if (AH_xc[c_AH][1]==0) reg_ah_r_(AH_theta0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]); 
        smooth_ah_r_(AH_theta0,AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        if (AH_xc[c_AH][1]==0) reg_ah_r_(AH_theta0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]); 
        eps0+=1;
      }
   }
   if (AH_xc[c_AH][1]==0) 
   { 
     reg_ah_r_(AH_theta0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
   }

   // compute resid as the root mean square of AH_theta values on the AH surface
   for (i=0, resid=0; i<np; i++) resid+=AH_theta0[i]*AH_theta0[i];
   resid=sqrt(resid/np);   

   if (!my_rank && dvtrace && num_trace<MAX_TRACE)
   {
      if (!((int)AH_ct[c_AH] % dvtrace))
      {
         num_trace++;
         AH_shape[0]=AH_Nchi[c_AH];
         AH_shape[1]=AH_Nphi[c_AH];
         AH_bbox[0]=0;
         AH_bbox[1]=M_PI;
         AH_bbox[2]=0;
         if (AH_xc[c_AH][1]<dy) {AH_bbox[3]=M_PI;} else {AH_bbox[3]=2*M_PI;}
         rank=2; 
  
         if (AH_xc[c_AH][0]<dx) {AH_bbox[0]=0; AH_bbox[1]=1; AH_bbox[2]=-1; AH_bbox[3]=1;} //planar BH 
 
         sprintf(name,"%sAH_%i_R_iter",AMRD_save_tag,c_AH+1);
         gft_out_bbox(name,AH_ct[c_AH],AH_shape,rank,AH_bbox,AH_R[c_AH]);
   
         sprintf(name,"%sAH_%i_theta_iter",AMRD_save_tag,c_AH+1);
         gft_out_bbox(name,AH_ct[c_AH],AH_shape,rank,AH_bbox,AH_theta0);
      }

      AH_ct[c_AH]++;
   }

//   //TEST7//
//   if (*c_polar>0) resid=0;  
//   //TEST7//

   return resid;
}

//-----------------------------------------------------------------------------
// driver routine
//-----------------------------------------------------------------------------
#define SMOOTH_R 0
#define SMOOTH_R_AFT 1
int find_apph(real *M, real *J, real *c_equat, real *c_polar, int use_R_ic, real *AH_min_resid)
{
   int iter,i,j,l,np,Lmax,Lmax_AH,is_ex;
   real resid,prev_resid,min_resid,min_R,c_R;
   int ltrace=2,numinc,skip_rest;
   real eps0,resid0,resid1,tol0,dx0[3],prev_AH_xc[2],area,eps1;
   static real min_tol,max_tol,first_c=1;
   int first=1;
   static int first_tol=1;

   *M=*J=0;
   *AH_min_resid=1e10;

   if (first_c) { for (l=0; l<MAX_BHS; l++) AH_ct[l]=0; first_c=0; }

   Lmax=PAMR_get_max_lev(PAMR_AMRH);
   if (Lmax>AH_Lmax[c_AH]) Lmax=AH_Lmax[c_AH];
   if (AH_Lmin[c_AH]>Lmax) return 0;

   np=AH_Nchi[c_AH]*AH_Nphi[c_AH];

   if (ltrace && my_rank==0)
   {
      printf("\n=========================================================================\n\n");
      printf("Searching for AH[%i] between L=%i and %i\n\n",c_AH+1,AH_Lmin[c_AH],Lmax);
   }

   eps0=eps1=AH_eps[c_AH];
   tol0=AH_tol[c_AH];
   if (eps1>1) eps1=eps1-1;

   // for found_AH==use_R_ic=0, figure out an initial guess for AH_R
   if (!use_R_ic)
   {
      if (AH_rsteps[c_AH]==1) for (i=0; i<np; i++) AH_R[c_AH][i]=AH_r0[c_AH]; // for one radius iteration
      else                                                                    // for several radii iterations
      {
         min_R=AH_r0[c_AH]; 
         min_resid=1e10;
         for (iter=1; iter<=(AH_rsteps[c_AH]+1); iter++)
         {
            c_R=AH_r0[c_AH]+(iter-1)*(AH_r1[c_AH]-AH_r0[c_AH])/(AH_rsteps[c_AH]-1);

            for (i=0; i<np; i++) AH_R[c_AH][i]=c_R;

            // compute initial theta values
            if (!fill_own(Lmax,ltrace,&first)) return 0;
            resid=fill_theta(AH_theta[c_AH],eps0,&area,c_equat,c_polar,&is_ex);

            if (is_ex)
            {
               printf("   is_ex on initial iteration\n");
               return 0;
            }
            if (resid<min_resid && !isnan(resid)) 
            {
               min_R=c_R;
               min_resid=resid;
            }
            if (ltrace>1 && my_rank==0)
               printf("   init iter %i: AH_r0=%lf ... |Theta|=%lf, Theta(mid)=%lf\n",iter,c_R,resid,AH_theta[c_AH][(np+1)/2]);
         }
         for (i=0; i<np; i++) AH_R[c_AH][i]=min_R;
      }
   }

   // save AH with minimum tol over all iterations in w2
   for (i=0; i<np; i++) AH_w2[c_AH][i]=AH_R[c_AH][i];
   for (i=0; i<2; i++) prev_AH_xc[i]=AH_xc[c_AH][i];

   int restarts=0;

   prev_resid=resid=tol0+1;
   iter=0;
   numinc=0;

   // (MAIN LOOP) while not yet AH_max_iter, or while residual is above tolerance
   while (iter < AH_max_iter[c_AH] && resid > tol0 && (numinc < AH_maxinc[c_AH] || AH_lambda[c_AH]>AH_lambda_min[c_AH]))   
   {
      iter++;
      skip_rest=0;

      if (numinc==AH_maxinc[c_AH])
      {
          AH_lambda[c_AH]*=(0.75);
          if (AH_lambda[c_AH]<AH_lambda_min[c_AH]) AH_lambda[c_AH]=AH_lambda_min[c_AH];
          if (my_rank==0 && ltrace) printf ("|Theta| increased for %i steps or out of bounds ... \n"
             "restoring AH shape to initial shape and decreasing lambda to %lf\n",numinc,AH_lambda[c_AH]);
          for (i=0; i<np; i++) AH_R[c_AH][i]=AH_w2[c_AH][i];
          for (i=0; i<2; i++) AH_xc[c_AH][i]=prev_AH_xc[i];
          numinc=0;
      }

      // do the adjustment within the iteration to keep changes small ... otherwise transformation is invalid
      adjust_ah_xc_(AH_R[c_AH],AH_xc[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH],&dx,&dy,&dz,&axisym);

      if (!fill_own(Lmax,ltrace,&first))
      {
         skip_rest=1;
         numinc=AH_maxinc[c_AH];
      }
      
      if (!skip_rest)
      {
      
         // compute theta values 
         resid=fill_theta(AH_theta[c_AH],eps0,&area,c_equat,c_polar,&is_ex);
      
         if (my_rank==0 && 0)
         {
            printf("iter=%i:\n\n",iter);
            for (i=0; i<AH_Nchi[c_AH]; i++)
               for (j=0; j<AH_Nphi[c_AH]; j++) 
                  printf("  i,j=%i,%i: theta=%lf, R=%lf\n",i,j,(AH_theta[c_AH])[i+j*AH_Nchi[c_AH]],(AH_R[c_AH])[i+j*AH_Nchi[c_AH]]);
         }

         if (is_ex)
         {
            if (ltrace>1 && my_rank==0) printf("is_ex iter %i: AH_R(mid)[%i]=%lf ... |Theta|=%lf, Theta(mid)=%lf \n",iter,(np+1)/2,AH_R[c_AH][(np+1)/2],resid,AH_theta[c_AH][(np+1)/2]);
            if (ltrace>0 && my_rank==0) printf("AH too close to excision boundary (iter=%i,min_resid=%lf)\n",iter,*AH_min_resid);

            resid=tol0+10; goto fin;
         }

         // update AH_R[i], for each i point on the AH surface, by a fraction AH_lambda of AH_theta
         //(NOTE: AH_R decreases when theta positive, and AH_R increases when theta negative)
         for (i=0; i<np; i++) AH_R[c_AH][i]+=(-AH_lambda[c_AH]*AH_theta[c_AH][i]);

         if (SMOOTH_R && eps1>0)
         {
            if (AH_xc[c_AH][1]==0) reg_ah_r_(AH_R[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_R[c_AH],AH_w1[c_AH],&eps1,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
         }

//LOOKING//
         /*if (iter>2 && resid>prev_resid) numinc++; else*/ numinc=0;

         prev_resid=resid;

         // print out AH_R(middle), |Theta|, Theta(middle), (xc,yc,zc)
         if (ltrace>1 && my_rank==0) printf("   iter %i: AH_R(mid)[%i]=%lf ... |Theta|=%lf, Theta(mid)=%lf, (xc,yc)=(%lf,%lf) \n",iter,(np+1)/2,AH_R[c_AH][(np+1)/2],resid,AH_theta[c_AH][(np+1)/2],AH_xc[c_AH][0],AH_xc[c_AH][1]);

         // checks for not-allowed negative AH_R 
         for (i=0; i<np; i++)
         {
            if (AH_R[c_AH][i]<0)
            {
              if (ltrace>0 && my_rank==0) printf("AH flowed to negative AH_R (iter=%i,min_resid=%lf)\n",iter,*AH_min_resid); resid=tol0+10; goto fin;
            }
         }

         if (resid < (*AH_min_resid) && !isnan(resid)) *AH_min_resid=resid;
      }
   }

   *M=3*M_PI/8*(pow(area/2/M_PI/M_PI,2.0/3.0))*(1+(pow(area/2/M_PI/M_PI,2.0/3.0))/AdS_L/AdS_L);

   if (ltrace && my_rank==0)
   {
      if (resid>tol0) 
      {
         printf("\n ... failed to find an AH in %i iterations\n",iter);
         if (iter<AH_max_iter[c_AH]) printf("     (stopped because of increase in residual ... min_resid=%lf)\n",*AH_min_resid);
      }
      else 
      {
         printf("\n ... found an AH (to within %lf) in %i iterations ... \n",
                tol0,iter);
         printf("     horizon mass: %5.3lf, from horizon area: %5.3lf, areal radius: %5.3lf \n",
              *M,area,cbrt(area/2/M_PI/M_PI));
         printf("     equat circum: %5.3lf, and  polar circum: %5.3lf \n",
              *c_equat,*c_polar);
      }

      printf("\n=========================================================================\n");
   }

   if (SMOOTH_R_AFT && eps1>0 && resid<=tol0)
   {
      if (AH_xc[c_AH][1]==0) reg_ah_r_(AH_R[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
      smooth_ah_r_(AH_R[c_AH],AH_w1[c_AH],&eps1,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
   }

   AH_ct[c_AH]+=100;

fin:

   if (resid>tol0)
   {
      for (i=0; i<np; i++) AH_R[c_AH][i]=AH_w2[c_AH][i];
      for (i=0; i<2; i++) AH_xc[c_AH][i]=prev_AH_xc[i];
      return 0; 
   }
   else return 1;
}






