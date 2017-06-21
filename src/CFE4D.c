//============================================================================
// in cartesian coordinates t,x,y for x in [-1,1], y in [0,1]
// using r=2*rho/(1-rho^2) compactification for rho=sqrt(x^2+y^2)
//
// application interface functions for CFE4D
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

//=============================================================================
// if axisym=1, then 2+1 simulation in (x,y) plane 
//=============================================================================
int axisym=1;

//=============================================================================
// set in fparam for now
//=============================================================================
real AdS_L;

//=============================================================================
// Carsten's "constraint-damping" parameters
//=============================================================================
real kappa_cd,rho_cd;

//=============================================================================
// id (and other) parameters
//=============================================================================

// gaussians
real phi1_amp_1,phi1_B_1,phi1_r0_1,phi1_delta_1,phi1_x0_1[2],phi1_ecc_1[2];
real phi1_amp_2,phi1_B_2,phi1_r0_2,phi1_delta_2,phi1_x0_2[2],phi1_ecc_2[2];

// if > 0, initialize with exact BH
real ief_bh_r0;

//gauge parameters:
int gauge_t;
int gauge_i;
real c1_t,c2_t,c3_t;
real c1_i,c2_i,c3_i;
real rho1_t,rho2_t,rho3_t,rho4_t,xi1_t,xi2_t,cbulk_t;
real rho1_i,rho2_i,rho3_i,rho4_i,xi1_i,xi2_i,cbulk_i;
real rhoa,rhob;

int cp_version; 

// excision parameters x-ex_xc  ex_r
real ex_rbuf[MAX_BHS];
int ex_reset_rbuf;
int ex_max_repop,ex_repop_buf,ex_repop_io;

// "internal" excision parameters, set by AH finder (eventually)
real ex_r[MAX_BHS][2],ex_xc[MAX_BHS][2];

int background,skip_constraints;
int output_ires,output_quasiset;

// new parameters in rtfile
int interptype,i_shift,regtype,stype;

int harmonize;

//extra dissipation
real diss_eps_k,diss_eps_y_cutoff;
int diss_kmax,diss_eps_k_cutoff_n,diss_bdy_k,diss_all_past_k,diss_all;

// excision parameters (more)
real ex_rbuf_a[MAX_BHS];

// AH parameters
int AH_Nchi[MAX_BHS],AH_Nphi[MAX_BHS],AH_Lmin[MAX_BHS],AH_Lmax[MAX_BHS],AH_find_best_fit[MAX_BHS];
int AH_max_iter[MAX_BHS],AH_freq[MAX_BHS],AH_freq_aft[MAX_BHS],AH_rsteps[MAX_BHS],AH_maxinc[MAX_BHS];
real AH_tol[MAX_BHS],AH_tol_aft[MAX_BHS],AH_r0[MAX_BHS],AH_lambda[MAX_BHS],AH_lambda_min[MAX_BHS];
real AH_eps[MAX_BHS],AH_r1[MAX_BHS],AH_tol_scale[MAX_BHS],AH_reset_scale[MAX_BHS];
real AH_xc[MAX_BHS][2],AH_max_tol_inc[MAX_BHS],AH_tmin[MAX_BHS],AH_omt_scale[MAX_BHS];
int use_AH_new_smooth,use_AH_new;
int c_AH;

//=============================================================================
// some convenient, "local" global variables
//=============================================================================

real *cl_res;

real *phi1,*phi1_n,*phi1_np1,*phi1_nm1; // MGH, AMRH n/np1/nm1

real *eb_xx,*eb_xx_n,*eb_xx_np1,*eb_xx_nm1; 
real *eb_xy,*eb_xy_n,*eb_xy_np1,*eb_xy_nm1; 
real *eb_xz,*eb_xz_n,*eb_xz_np1,*eb_xz_nm1; 
real *eb_yy,*eb_yy_n,*eb_yy_np1,*eb_yy_nm1; 
real *eb_yz,*eb_yz_n,*eb_yz_np1,*eb_yz_nm1; 
real *eb_zz,*eb_zz_n,*eb_zz_np1,*eb_zz_nm1; 

real *gb_tt,*gb_tt_n,*gb_tt_np1,*gb_tt_nm1; 
real *gb_tx,*gb_tx_n,*gb_tx_np1,*gb_tx_nm1; 
real *gb_ty,*gb_ty_n,*gb_ty_np1,*gb_ty_nm1; 
real *gb_xx,*gb_xx_n,*gb_xx_np1,*gb_xx_nm1; 
real *gb_xy,*gb_xy_n,*gb_xy_np1,*gb_xy_nm1; 
real *gb_yy,*gb_yy_n,*gb_yy_np1,*gb_yy_nm1; 
real *psi,*psi_n,*psi_np1,*psi_nm1; 

real *Hb_t,*Hb_t_n,*Hb_t_np1,*Hb_t_nm1;
real *Hb_x,*Hb_x_n,*Hb_x_np1,*Hb_x_nm1;
real *Hb_y,*Hb_y_n,*Hb_y_np1,*Hb_y_nm1;

real *phi1_t,*phi1_t_n;

real *eb_xx_t,*eb_xx_t_n;
real *eb_xy_t,*eb_xy_t_n;
real *eb_xz_t,*eb_xz_t_n;
real *eb_yy_t,*eb_yy_t_n;
real *eb_yz_t,*eb_yz_t_n;
real *eb_zz_t,*eb_zz_t_n;

real *gb_tt_t,*gb_tt_t_n;
real *gb_tx_t,*gb_tx_t_n;
real *gb_ty_t,*gb_ty_t_n;
real *gb_xx_t,*gb_xx_t_n;
real *gb_xy_t,*gb_xy_t_n;
real *gb_yy_t,*gb_yy_t_n;
real *psi_t,*psi_t_n;
real *Hb_t_t,*Hb_t_t_n;
real *Hb_x_t,*Hb_x_t_n;
real *Hb_y_t,*Hb_y_t_n;

real *w1,*mg_w1;
real *w2,*mg_w2;
real *w3,*mg_w3;
real *w4,*mg_w4;

real *mask,*mask_mg,*chr,*chr_mg;
real *kg_ires,*alpha,*ricci,*theta,*f,*K;

real *phi1_res;
real *eb_res;
real *gb_res;

real *efe_all_ires;
real *efe_tt_ires,*efe_tx_ires,*efe_ty_ires;
real *efe_xx_ires,*efe_xy_ires,*efe_yy_ires,*efe_psi_ires;
real *quasiset_tt,*quasiset_tx,*quasiset_ty;
real *quasiset_xx,*quasiset_xy,*quasiset_yy;
real *quasiset_psi;
real *quasiset_mass;
real *kretsch;

real *tfunction,*test1,*test2,*test3,*test4;
real *iresall,*irestt,*irestx,*iresty,*iresxx,*iresxy,*iresyy,*irespsi;
real *qstt,*qstx,*qsty,*qsxx,*qsxy,*qsyy,*qspsi;
real *qsmass;
real *qsone;

real *zeta,*zeta_res,*zeta_lop,*zeta_rhs;

real *hb_t_res,*hb_i_res;
real *Hb_t_0,*Hb_x_0,*Hb_y_0;

real *g_norms;

real *x,*y,*z;
int shape[2],ghost_width[4],Nx,Ny,Nz,phys_bdy[4],size,g_rank;
real base_bbox[4],bbox[4],dx,dy,dz,dt,dx_Lc;
int g_L;

int cl_res_gfn;

int phi1_gfn,phi1_n_gfn,phi1_np1_gfn,phi1_nm1_gfn; 

int eb_xx_gfn,eb_xx_n_gfn,eb_xx_np1_gfn,eb_xx_nm1_gfn; 
int eb_xy_gfn,eb_xy_n_gfn,eb_xy_np1_gfn,eb_xy_nm1_gfn; 
int eb_xz_gfn,eb_xz_n_gfn,eb_xz_np1_gfn,eb_xz_nm1_gfn; 
int eb_yy_gfn,eb_yy_n_gfn,eb_yy_np1_gfn,eb_yy_nm1_gfn; 
int eb_yz_gfn,eb_yz_n_gfn,eb_yz_np1_gfn,eb_yz_nm1_gfn; 
int eb_zz_gfn,eb_zz_n_gfn,eb_zz_np1_gfn,eb_zz_nm1_gfn; 

int gb_tt_gfn,gb_tt_n_gfn,gb_tt_np1_gfn,gb_tt_nm1_gfn; 
int gb_tx_gfn,gb_tx_n_gfn,gb_tx_np1_gfn,gb_tx_nm1_gfn; 
int gb_ty_gfn,gb_ty_n_gfn,gb_ty_np1_gfn,gb_ty_nm1_gfn; 
int gb_xx_gfn,gb_xx_n_gfn,gb_xx_np1_gfn,gb_xx_nm1_gfn; 
int gb_xy_gfn,gb_xy_n_gfn,gb_xy_np1_gfn,gb_xy_nm1_gfn; 
int gb_yy_gfn,gb_yy_n_gfn,gb_yy_np1_gfn,gb_yy_nm1_gfn; 
int psi_gfn,psi_n_gfn,psi_np1_gfn,psi_nm1_gfn; 

int Hb_t_gfn,Hb_t_n_gfn,Hb_t_np1_gfn,Hb_t_nm1_gfn;
int Hb_x_gfn,Hb_x_n_gfn,Hb_x_np1_gfn,Hb_x_nm1_gfn;
int Hb_y_gfn,Hb_y_n_gfn,Hb_y_np1_gfn,Hb_y_nm1_gfn;

int phi1_t_gfn,phi1_t_n_gfn;

int eb_xx_t_gfn,eb_xx_t_n_gfn;
int eb_xy_t_gfn,eb_xy_t_n_gfn;
int eb_xz_t_gfn,eb_xz_t_n_gfn;
int eb_yy_t_gfn,eb_yy_t_n_gfn;
int eb_yz_t_gfn,eb_yz_t_n_gfn;
int eb_zz_t_gfn,eb_zz_t_n_gfn;

int gb_tt_t_gfn,gb_tt_t_n_gfn;
int gb_tx_t_gfn,gb_tx_t_n_gfn;
int gb_ty_t_gfn,gb_ty_t_n_gfn;
int gb_xx_t_gfn,gb_xx_t_n_gfn;
int gb_xy_t_gfn,gb_xy_t_n_gfn;
int gb_yy_t_gfn,gb_yy_t_n_gfn;
int psi_t_gfn,psi_t_n_gfn;
int Hb_t_t_gfn,Hb_t_t_n_gfn;
int Hb_x_t_gfn,Hb_x_t_n_gfn;
int Hb_y_t_gfn,Hb_y_t_n_gfn;

int mask_gfn,mask_mg_gfn,chr_gfn,chr_mg_gfn;
real *gu_tt,*gu_tx,*gu_ty,*gu_xx,*gu_xy,*gu_yy,*m_g_det;
int kg_ires_gfn,alpha_gfn,theta_gfn,f_gfn;

int phi1_res_gfn;
int eb_res_gfn;
int gb_res_gfn;
int efe_all_ires_gfn;
int efe_tt_ires_gfn,efe_tx_ires_gfn,efe_ty_ires_gfn;
int efe_xx_ires_gfn,efe_xy_ires_gfn,efe_yy_ires_gfn,efe_psi_ires_gfn;
int quasiset_tt_gfn,quasiset_tx_gfn,quasiset_ty_gfn;
int quasiset_xx_gfn,quasiset_xy_gfn,quasiset_yy_gfn;
int quasiset_psi_gfn;
int quasiset_mass_gfn;
int kretsch_gfn;

int tfunction_gfn,test1_gfn,test2_gfn,test3_gfn,test4_gfn;
int iresall_gfn,irestt_gfn,irestx_gfn,iresty_gfn,iresxx_gfn,iresxy_gfn,iresyy_gfn,irespsi_gfn;
int qstt_gfn,qstx_gfn,qsty_gfn,qsxx_gfn,qsxy_gfn,qsyy_gfn,qspsi_gfn;
int qsmass_gfn;
int qsone_gfn;

int zeta_gfn,zeta_res_gfn,zeta_lop_gfn,zeta_rhs_gfn;

int hb_t_res_gfn,hb_i_res_gfn;
int Hb_t_0_gfn,Hb_x_0_gfn,Hb_y_0_gfn;

int w1_gfn,mg_w1_gfn;
int w2_gfn,mg_w2_gfn;
int w3_gfn,mg_w3_gfn;
int w4_gfn,mg_w4_gfn;

int skip_ires=0;
int skip_exp=0;

int gu_tt_gfn,gu_tx_gfn,gu_ty_gfn,gu_xx_gfn;
int gu_xy_gfn,gu_yy_gfn,m_g_det_gfn;

//=============================================================================
// arrays holding AH shape and other hypersurface info
//=============================================================================
real *AH_theta[MAX_BHS],*AH_R[MAX_BHS],*AH_w1[MAX_BHS],*AH_w2[MAX_BHS],*AH_w3[MAX_BHS];
real *AH_theta_ads[MAX_BHS];
int *AH_lev[MAX_BHS],*AH_own[MAX_BHS];

//=============================================================================
// call after variables have been defined
//=============================================================================
void set_gfns(void)
{
    if ((cl_res_gfn   = PAMR_get_gfn("cl_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((phi1_gfn     = PAMR_get_gfn("phi1",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((phi1_nm1_gfn = PAMR_get_gfn("phi1",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((phi1_n_gfn   = PAMR_get_gfn("phi1",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((phi1_np1_gfn = PAMR_get_gfn("phi1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((eb_xx_gfn     = PAMR_get_gfn("eb_xx",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_xx_nm1_gfn = PAMR_get_gfn("eb_xx",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_xx_n_gfn   = PAMR_get_gfn("eb_xx",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_xx_np1_gfn = PAMR_get_gfn("eb_xx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((eb_xy_gfn     = PAMR_get_gfn("eb_xy",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_xy_nm1_gfn = PAMR_get_gfn("eb_xy",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_xy_n_gfn   = PAMR_get_gfn("eb_xy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_xy_np1_gfn = PAMR_get_gfn("eb_xy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((eb_xz_gfn     = PAMR_get_gfn("eb_xz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_xz_nm1_gfn = PAMR_get_gfn("eb_xz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_xz_n_gfn   = PAMR_get_gfn("eb_xz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_xz_np1_gfn = PAMR_get_gfn("eb_xz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((eb_yy_gfn     = PAMR_get_gfn("eb_yy",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_yy_nm1_gfn = PAMR_get_gfn("eb_yy",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_yy_n_gfn   = PAMR_get_gfn("eb_yy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_yy_np1_gfn = PAMR_get_gfn("eb_yy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((eb_yz_gfn     = PAMR_get_gfn("eb_yz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_yz_nm1_gfn = PAMR_get_gfn("eb_yz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_yz_n_gfn   = PAMR_get_gfn("eb_yz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_yz_np1_gfn = PAMR_get_gfn("eb_yz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((eb_zz_gfn     = PAMR_get_gfn("eb_zz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_zz_nm1_gfn = PAMR_get_gfn("eb_zz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_zz_n_gfn   = PAMR_get_gfn("eb_zz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_zz_np1_gfn = PAMR_get_gfn("eb_zz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_tt_gfn     = PAMR_get_gfn("gb_tt",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_nm1_gfn = PAMR_get_gfn("gb_tt",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_n_gfn   = PAMR_get_gfn("gb_tt",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_np1_gfn = PAMR_get_gfn("gb_tt",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_tx_gfn     = PAMR_get_gfn("gb_tx",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_nm1_gfn = PAMR_get_gfn("gb_tx",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_n_gfn   = PAMR_get_gfn("gb_tx",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_np1_gfn = PAMR_get_gfn("gb_tx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_ty_gfn     = PAMR_get_gfn("gb_ty",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_ty_nm1_gfn = PAMR_get_gfn("gb_ty",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_ty_n_gfn   = PAMR_get_gfn("gb_ty",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_ty_np1_gfn = PAMR_get_gfn("gb_ty",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_xx_gfn     = PAMR_get_gfn("gb_xx",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_nm1_gfn = PAMR_get_gfn("gb_xx",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_n_gfn   = PAMR_get_gfn("gb_xx",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_np1_gfn = PAMR_get_gfn("gb_xx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_xy_gfn     = PAMR_get_gfn("gb_xy",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xy_nm1_gfn = PAMR_get_gfn("gb_xy",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xy_n_gfn   = PAMR_get_gfn("gb_xy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xy_np1_gfn = PAMR_get_gfn("gb_xy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_yy_gfn     = PAMR_get_gfn("gb_yy",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yy_nm1_gfn = PAMR_get_gfn("gb_yy",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yy_n_gfn   = PAMR_get_gfn("gb_yy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yy_np1_gfn = PAMR_get_gfn("gb_yy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((psi_gfn     = PAMR_get_gfn("psi",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_nm1_gfn = PAMR_get_gfn("psi",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_n_gfn   = PAMR_get_gfn("psi",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_np1_gfn = PAMR_get_gfn("psi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((Hb_t_gfn      = PAMR_get_gfn("Hb_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_nm1_gfn  = PAMR_get_gfn("Hb_t",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_n_gfn    = PAMR_get_gfn("Hb_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_np1_gfn  = PAMR_get_gfn("Hb_t",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((phi1_t_gfn   = PAMR_get_gfn("phi1_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((phi1_t_n_gfn = PAMR_get_gfn("phi1_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((eb_xx_t_gfn   = PAMR_get_gfn("eb_xx_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_xx_t_n_gfn = PAMR_get_gfn("eb_xx_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_xy_t_gfn   = PAMR_get_gfn("eb_xy_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_xy_t_n_gfn = PAMR_get_gfn("eb_xy_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_xz_t_gfn   = PAMR_get_gfn("eb_xz_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_xz_t_n_gfn = PAMR_get_gfn("eb_xz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_yy_t_gfn   = PAMR_get_gfn("eb_yy_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_yy_t_n_gfn = PAMR_get_gfn("eb_yy_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_yz_t_gfn   = PAMR_get_gfn("eb_yz_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_yz_t_n_gfn = PAMR_get_gfn("eb_yz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_zz_t_gfn   = PAMR_get_gfn("eb_zz_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((eb_zz_t_n_gfn = PAMR_get_gfn("eb_zz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_tt_t_gfn   = PAMR_get_gfn("gb_tt_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_t_n_gfn = PAMR_get_gfn("gb_tt_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_t_gfn   = PAMR_get_gfn("gb_tx_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_t_n_gfn = PAMR_get_gfn("gb_tx_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_ty_t_gfn   = PAMR_get_gfn("gb_ty_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_ty_t_n_gfn = PAMR_get_gfn("gb_ty_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_t_gfn   = PAMR_get_gfn("gb_xx_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_t_n_gfn = PAMR_get_gfn("gb_xx_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xy_t_gfn   = PAMR_get_gfn("gb_xy_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xy_t_n_gfn = PAMR_get_gfn("gb_xy_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yy_t_gfn   = PAMR_get_gfn("gb_yy_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yy_t_n_gfn = PAMR_get_gfn("gb_yy_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_t_gfn     = PAMR_get_gfn("psi_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_t_n_gfn   = PAMR_get_gfn("psi_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_t_gfn    = PAMR_get_gfn("Hb_t_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_t_n_gfn  = PAMR_get_gfn("Hb_t_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_t_gfn    = PAMR_get_gfn("Hb_x_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_t_n_gfn  = PAMR_get_gfn("Hb_x_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_t_gfn    = PAMR_get_gfn("Hb_y_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_t_n_gfn  = PAMR_get_gfn("Hb_y_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((Hb_x_gfn      = PAMR_get_gfn("Hb_x",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_nm1_gfn  = PAMR_get_gfn("Hb_x",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_n_gfn    = PAMR_get_gfn("Hb_x",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_np1_gfn  = PAMR_get_gfn("Hb_x",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((Hb_y_gfn      = PAMR_get_gfn("Hb_y",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_nm1_gfn  = PAMR_get_gfn("Hb_y",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_n_gfn    = PAMR_get_gfn("Hb_y",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_np1_gfn  = PAMR_get_gfn("Hb_y",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((zeta_gfn     = PAMR_get_gfn("zeta",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((zeta_res_gfn = PAMR_get_gfn("zeta_res",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((zeta_lop_gfn = PAMR_get_gfn("zeta_lop",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((zeta_rhs_gfn = PAMR_get_gfn("zeta_rhs",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    if ((mask_mg_gfn = PAMR_get_gfn("cmask",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_gfn    = PAMR_get_gfn("cmask",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((chr_gfn     = PAMR_get_gfn("chr",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((chr_mg_gfn  = PAMR_get_gfn("chr",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    if ((phi1_res_gfn  = PAMR_get_gfn("phi1_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((eb_res_gfn    = PAMR_get_gfn("eb_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_res_gfn    = PAMR_get_gfn("gb_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((efe_all_ires_gfn    = PAMR_get_gfn("efe_all_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_tt_ires_gfn    = PAMR_get_gfn("efe_tt_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_tx_ires_gfn    = PAMR_get_gfn("efe_tx_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_ty_ires_gfn    = PAMR_get_gfn("efe_ty_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_xx_ires_gfn    = PAMR_get_gfn("efe_xx_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_xy_ires_gfn    = PAMR_get_gfn("efe_xy_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_yy_ires_gfn    = PAMR_get_gfn("efe_yy_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_psi_ires_gfn    = PAMR_get_gfn("efe_psi_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_tt_gfn    = PAMR_get_gfn("quasiset_tt",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_tx_gfn    = PAMR_get_gfn("quasiset_tx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_ty_gfn    = PAMR_get_gfn("quasiset_ty",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_xx_gfn    = PAMR_get_gfn("quasiset_xx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_xy_gfn    = PAMR_get_gfn("quasiset_xy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_yy_gfn    = PAMR_get_gfn("quasiset_yy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_psi_gfn    = PAMR_get_gfn("quasiset_psi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_mass_gfn    = PAMR_get_gfn("quasiset_mass",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((kretsch_gfn    = PAMR_get_gfn("kretsch",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((hb_t_res_gfn  = PAMR_get_gfn("hb_t_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((hb_i_res_gfn  = PAMR_get_gfn("hb_i_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_0_gfn  = PAMR_get_gfn("Hb_t_0",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_0_gfn  = PAMR_get_gfn("Hb_x_0",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_0_gfn  = PAMR_get_gfn("Hb_y_0",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((w1_gfn   = PAMR_get_gfn("w1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((w2_gfn   = PAMR_get_gfn("w2",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((w3_gfn   = PAMR_get_gfn("w3",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((w4_gfn   = PAMR_get_gfn("w4",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((kg_ires_gfn= PAMR_get_gfn("kg_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((alpha_gfn= PAMR_get_gfn("alpha",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((theta_gfn= PAMR_get_gfn("theta",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((f_gfn= PAMR_get_gfn("f",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((tfunction_gfn  = PAMR_get_gfn("tfunction",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((test1_gfn  = PAMR_get_gfn("test1",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((test2_gfn  = PAMR_get_gfn("test2",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((test3_gfn  = PAMR_get_gfn("test3",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((test4_gfn  = PAMR_get_gfn("test4",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((iresall_gfn  = PAMR_get_gfn("iresall",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((irestt_gfn  = PAMR_get_gfn("irestt",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((irestx_gfn  = PAMR_get_gfn("irestx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((iresty_gfn  = PAMR_get_gfn("iresty",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((iresxx_gfn  = PAMR_get_gfn("iresxx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((iresxy_gfn  = PAMR_get_gfn("iresxy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((iresyy_gfn  = PAMR_get_gfn("iresyy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((irespsi_gfn  = PAMR_get_gfn("irespsi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qstt_gfn  = PAMR_get_gfn("qstt",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qstx_gfn  = PAMR_get_gfn("qstx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qsty_gfn  = PAMR_get_gfn("qsty",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qsxx_gfn  = PAMR_get_gfn("qsxx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qsxy_gfn  = PAMR_get_gfn("qsxy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qsyy_gfn  = PAMR_get_gfn("qsyy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qspsi_gfn  = PAMR_get_gfn("qspsi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qsmass_gfn  = PAMR_get_gfn("qsmass",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qsone_gfn  = PAMR_get_gfn("qsone",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((mg_w1_gfn   = PAMR_get_gfn("mg_w1",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mg_w2_gfn   = PAMR_get_gfn("mg_w2",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mg_w3_gfn   = PAMR_get_gfn("mg_w3",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mg_w4_gfn   = PAMR_get_gfn("mg_w4",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    if ((gu_tt_gfn   = PAMR_get_gfn("gu_tt",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gu_tx_gfn   = PAMR_get_gfn("gu_tx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gu_ty_gfn   = PAMR_get_gfn("gu_ty",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gu_xx_gfn   = PAMR_get_gfn("gu_xx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gu_xy_gfn   = PAMR_get_gfn("gu_xy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gu_yy_gfn   = PAMR_get_gfn("gu_yy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((m_g_det_gfn = PAMR_get_gfn("m_g_det",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    g_norms=AMRD_get_global_norms();
}

//=============================================================================
// call with valid iter to set up globals:
//=============================================================================
void ldptr_bbox(void)
{
   real dx0[2];
   static int first=1;

   if (first) 
   {
      first=0; 
      set_gfns();
      PAMR_get_global_bbox(base_bbox);
      if (PAMR_get_max_lev(PAMR_AMRH)>1) PAMR_get_dxdt(2,dx0,&dt); else PAMR_get_dxdt(1,dx0,&dt);
      dx_Lc=dx0[0];
   }

   PAMR_get_g_rank(&g_rank);
   PAMR_get_g_shape(shape);
   PAMR_get_g_bbox(bbox);
   PAMR_get_g_ghost_width(ghost_width);
   PAMR_get_g_level(&g_L);
   PAMR_get_dxdt(g_L,dx0,&dt);
   dx=dx0[0];
   dy=dx0[1];
   dz=dx;

   if ((bbox[0]-base_bbox[0])<dx/2) phys_bdy[0]=1; else phys_bdy[0]=0;
   if ((base_bbox[1]-bbox[1])<dx/2) phys_bdy[1]=1; else phys_bdy[1]=0;
   if ((bbox[2]-base_bbox[2])<dy/2) phys_bdy[2]=1; else phys_bdy[2]=0;
   if ((base_bbox[3]-bbox[3])<dy/2) phys_bdy[3]=1; else phys_bdy[3]=0;

   Nx=shape[0];
   Ny=shape[1];
   Nz=1;

   size=Nx*Ny;
}

void ldptr(void)
{
   real *x0[2],*gfs[PAMR_MAX_GFNS];

   ldptr_bbox();

   PAMR_get_g_x(x0);

   x=x0[0];
   y=x0[1]; 

   PAMR_get_g_gfs(gfs);

   cl_res   = gfs[cl_res_gfn-1];

   phi1     = gfs[phi1_gfn-1];
   phi1_n   = gfs[phi1_n_gfn-1];
   phi1_np1 = gfs[phi1_np1_gfn-1];
   phi1_nm1 = gfs[phi1_nm1_gfn-1];

   eb_xx     = gfs[eb_xx_gfn-1];
   eb_xx_n   = gfs[eb_xx_n_gfn-1];
   eb_xx_np1 = gfs[eb_xx_np1_gfn-1];
   eb_xx_nm1 = gfs[eb_xx_nm1_gfn-1];

   eb_xy     = gfs[eb_xy_gfn-1];
   eb_xy_n   = gfs[eb_xy_n_gfn-1];
   eb_xy_np1 = gfs[eb_xy_np1_gfn-1];
   eb_xy_nm1 = gfs[eb_xy_nm1_gfn-1];

   eb_xz     = gfs[eb_xz_gfn-1];
   eb_xz_n   = gfs[eb_xz_n_gfn-1];
   eb_xz_np1 = gfs[eb_xz_np1_gfn-1];
   eb_xz_nm1 = gfs[eb_xz_nm1_gfn-1];

   eb_yy     = gfs[eb_yy_gfn-1];
   eb_yy_n   = gfs[eb_yy_n_gfn-1];
   eb_yy_np1 = gfs[eb_yy_np1_gfn-1];
   eb_yy_nm1 = gfs[eb_yy_nm1_gfn-1];

   eb_yz     = gfs[eb_yz_gfn-1];
   eb_yz_n   = gfs[eb_yz_n_gfn-1];
   eb_yz_np1 = gfs[eb_yz_np1_gfn-1];
   eb_yz_nm1 = gfs[eb_yz_nm1_gfn-1];

   eb_zz     = gfs[eb_zz_gfn-1];
   eb_zz_n   = gfs[eb_zz_n_gfn-1];
   eb_zz_np1 = gfs[eb_zz_np1_gfn-1];
   eb_zz_nm1 = gfs[eb_zz_nm1_gfn-1];

   gb_tt     = gfs[gb_tt_gfn-1];
   gb_tt_n   = gfs[gb_tt_n_gfn-1];
   gb_tt_np1 = gfs[gb_tt_np1_gfn-1];
   gb_tt_nm1 = gfs[gb_tt_nm1_gfn-1];

   gb_tx     = gfs[gb_tx_gfn-1];
   gb_tx_n   = gfs[gb_tx_n_gfn-1];
   gb_tx_np1 = gfs[gb_tx_np1_gfn-1];
   gb_tx_nm1 = gfs[gb_tx_nm1_gfn-1];

   gb_ty     = gfs[gb_ty_gfn-1];
   gb_ty_n   = gfs[gb_ty_n_gfn-1];
   gb_ty_np1 = gfs[gb_ty_np1_gfn-1];
   gb_ty_nm1 = gfs[gb_ty_nm1_gfn-1];

   gb_xx     = gfs[gb_xx_gfn-1];
   gb_xx_n   = gfs[gb_xx_n_gfn-1];
   gb_xx_np1 = gfs[gb_xx_np1_gfn-1];
   gb_xx_nm1 = gfs[gb_xx_nm1_gfn-1];

   gb_xy     = gfs[gb_xy_gfn-1];
   gb_xy_n   = gfs[gb_xy_n_gfn-1];
   gb_xy_np1 = gfs[gb_xy_np1_gfn-1];
   gb_xy_nm1 = gfs[gb_xy_nm1_gfn-1];

   gb_yy     = gfs[gb_yy_gfn-1];
   gb_yy_n   = gfs[gb_yy_n_gfn-1];
   gb_yy_np1 = gfs[gb_yy_np1_gfn-1];
   gb_yy_nm1 = gfs[gb_yy_nm1_gfn-1];

   psi     = gfs[psi_gfn-1];
   psi_n   = gfs[psi_n_gfn-1];
   psi_np1 = gfs[psi_np1_gfn-1];
   psi_nm1 = gfs[psi_nm1_gfn-1];

   Hb_t      = gfs[Hb_t_gfn-1];
   Hb_t_n    = gfs[Hb_t_n_gfn-1];
   Hb_t_nm1  = gfs[Hb_t_nm1_gfn-1];
   Hb_t_np1  = gfs[Hb_t_np1_gfn-1];

   Hb_x      = gfs[Hb_x_gfn-1];
   Hb_x_n    = gfs[Hb_x_n_gfn-1];
   Hb_x_nm1  = gfs[Hb_x_nm1_gfn-1];
   Hb_x_np1  = gfs[Hb_x_np1_gfn-1];
   Hb_y      = gfs[Hb_y_gfn-1];
   Hb_y_n    = gfs[Hb_y_n_gfn-1];
   Hb_y_nm1  = gfs[Hb_y_nm1_gfn-1];
   Hb_y_np1  = gfs[Hb_y_np1_gfn-1];

   phi1_t   = gfs[phi1_t_gfn-1];
   phi1_t_n = gfs[phi1_t_n_gfn-1];

   eb_xx_t   = gfs[eb_xx_t_gfn-1];
   eb_xx_t_n = gfs[eb_xx_t_n_gfn-1];
   eb_xy_t   = gfs[eb_xy_t_gfn-1];
   eb_xy_t_n = gfs[eb_xy_t_n_gfn-1];
   eb_xz_t   = gfs[eb_xz_t_gfn-1];
   eb_xz_t_n = gfs[eb_xz_t_n_gfn-1];
   eb_yy_t   = gfs[eb_yy_t_gfn-1];
   eb_yy_t_n = gfs[eb_yy_t_n_gfn-1];
   eb_yz_t   = gfs[eb_yz_t_gfn-1];
   eb_yz_t_n = gfs[eb_yz_t_n_gfn-1];
   eb_zz_t   = gfs[eb_zz_t_gfn-1];
   eb_zz_t_n = gfs[eb_zz_t_n_gfn-1];

   gb_tt_t   = gfs[gb_tt_t_gfn-1];
   gb_tt_t_n = gfs[gb_tt_t_n_gfn-1];
   gb_tx_t   = gfs[gb_tx_t_gfn-1];
   gb_tx_t_n = gfs[gb_tx_t_n_gfn-1];
   gb_ty_t   = gfs[gb_ty_t_gfn-1];
   gb_ty_t_n = gfs[gb_ty_t_n_gfn-1];
   gb_xx_t   = gfs[gb_xx_t_gfn-1];
   gb_xx_t_n = gfs[gb_xx_t_n_gfn-1];
   gb_xy_t   = gfs[gb_xy_t_gfn-1];
   gb_xy_t_n = gfs[gb_xy_t_n_gfn-1];
   gb_yy_t   = gfs[gb_yy_t_gfn-1];
   gb_yy_t_n = gfs[gb_yy_t_n_gfn-1];
   psi_t     = gfs[psi_t_gfn-1];
   psi_t_n   = gfs[psi_t_n_gfn-1];
   Hb_t_t    = gfs[Hb_t_t_gfn-1];
   Hb_t_t_n  = gfs[Hb_t_t_n_gfn-1];
   Hb_x_t    = gfs[Hb_x_t_gfn-1];
   Hb_x_t_n  = gfs[Hb_x_t_n_gfn-1];
   Hb_y_t    = gfs[Hb_y_t_gfn-1];
   Hb_y_t_n  = gfs[Hb_y_t_n_gfn-1];

   zeta     = gfs[zeta_gfn-1];
   zeta_lop = gfs[zeta_lop_gfn-1];
   zeta_res = gfs[zeta_res_gfn-1];
   zeta_rhs = gfs[zeta_rhs_gfn-1];

   mask    = gfs[mask_gfn-1];
   mask_mg = gfs[mask_mg_gfn-1];
   chr = gfs[chr_gfn-1]; 
   chr_mg = gfs[chr_mg_gfn-1]; 

   phi1_res  = gfs[phi1_res_gfn-1];

   eb_res    = gfs[eb_res_gfn-1];

   gb_res    = gfs[gb_res_gfn-1];

   efe_all_ires  = gfs[efe_all_ires_gfn-1];
   efe_tt_ires  = gfs[efe_tt_ires_gfn-1];
   efe_tx_ires  = gfs[efe_tx_ires_gfn-1];
   efe_ty_ires  = gfs[efe_ty_ires_gfn-1];
   efe_xx_ires  = gfs[efe_xx_ires_gfn-1];
   efe_xy_ires  = gfs[efe_xy_ires_gfn-1];
   efe_yy_ires  = gfs[efe_yy_ires_gfn-1];
   efe_psi_ires  = gfs[efe_psi_ires_gfn-1];
   quasiset_tt  = gfs[quasiset_tt_gfn-1];
   quasiset_tx  = gfs[quasiset_tx_gfn-1];
   quasiset_ty  = gfs[quasiset_ty_gfn-1];
   quasiset_xx  = gfs[quasiset_xx_gfn-1];
   quasiset_xy  = gfs[quasiset_xy_gfn-1];
   quasiset_yy  = gfs[quasiset_yy_gfn-1];
   quasiset_psi  = gfs[quasiset_psi_gfn-1];
   quasiset_mass  = gfs[quasiset_mass_gfn-1];
   kretsch  = gfs[kretsch_gfn-1];
   hb_t_res  = gfs[hb_t_res_gfn-1];
   hb_i_res  = gfs[hb_i_res_gfn-1];
   Hb_t_0  = gfs[Hb_t_0_gfn-1];
   Hb_x_0  = gfs[Hb_x_0_gfn-1];
   Hb_y_0  = gfs[Hb_y_0_gfn-1];
   w1 = gfs[w1_gfn-1];
   w2 = gfs[w2_gfn-1];
   w3 = gfs[w3_gfn-1];
   w4 = gfs[w4_gfn-1];
   kg_ires = gfs[kg_ires_gfn-1];
   alpha = gfs[alpha_gfn-1];
   theta = gfs[theta_gfn-1];
   f = gfs[f_gfn-1];

   tfunction = gfs[tfunction_gfn-1];
   test1 = gfs[test1_gfn-1];
   test2 = gfs[test2_gfn-1];
   test3 = gfs[test3_gfn-1];
   test4 = gfs[test4_gfn-1];
   iresall = gfs[iresall_gfn-1];
   irestt = gfs[irestt_gfn-1];
   irestx = gfs[irestx_gfn-1];
   iresty = gfs[iresty_gfn-1];
   iresxx = gfs[iresxx_gfn-1];
   iresxy = gfs[iresxy_gfn-1];
   iresyy = gfs[iresyy_gfn-1];
   irespsi = gfs[irespsi_gfn-1];
   qstt = gfs[qstt_gfn-1];
   qstx = gfs[qstx_gfn-1];
   qsty = gfs[qsty_gfn-1];
   qsxx = gfs[qsxx_gfn-1];
   qsxy = gfs[qsxy_gfn-1];
   qsyy = gfs[qsyy_gfn-1];
   qspsi = gfs[qspsi_gfn-1];
   qsmass = gfs[qsmass_gfn-1];
   qsone = gfs[qsone_gfn-1];

   mg_w1    =gfs[mg_w1_gfn-1]; 
   mg_w2    =gfs[mg_w2_gfn-1]; 
   mg_w3    =gfs[mg_w3_gfn-1]; 
   mg_w4    =gfs[mg_w4_gfn-1]; 

   gu_tt = gfs[gu_tt_gfn-1];
   gu_tx = gfs[gu_tx_gfn-1];
   gu_ty = gfs[gu_ty_gfn-1];
   gu_xx = gfs[gu_xx_gfn-1];
   gu_xy = gfs[gu_xy_gfn-1];
   gu_yy = gfs[gu_yy_gfn-1];
   m_g_det = gfs[m_g_det_gfn-1];

}

//=============================================================================
// PAMR_get_dxdt() only works with AMR hierarchy levels ... here we use
// lambda for dt, but this only works if rhosp=rhotm
//=============================================================================
void ldptr_mg(void)
{
   real lambda;

   ldptr();

   dx=x[1]-x[0]; dy=y[1]-y[0];
   PAMR_get_lambda(&lambda);
   dt=lambda*dx;
}

//=============================================================================
// utility routines
//=============================================================================
void const_f(real *f, real c)
{
   int i;

   for (i=0; i<Nx*Ny; i++) f[i]=c;
}

void zero_f(real *f)
{
   const_f(f,0);
}

void zero_f_ex(real *f, real *chr)
{
   int i;

   for (i=0; i<Nx*Ny; i++) if (chr[i]==AMRD_ex) f[i]=0;
}

real norm_l2(real *f, real *cmask, real *chr)
{
   int i;
   real norm=0;
   int sum=0;

   for (i=0; i<Nx*Ny; i++) 
      if (cmask[i]==AMRD_CMASK_ON && (chr[i]!=AMRD_ex)) { sum++; norm+=f[i]*f[i]; }

   if (!sum) sum=1;
   return (sqrt(norm/sum));
}

//=============================================================================
// the following zeros the AH_max_iter and ex_r (semiaxes of best-fit ellipses) 
// of engulfed AHs
//=============================================================================
void remove_redundant_AH()
{
   int i1,i2,is_1in2,is_2in1,intr;
   real d,min_bhr,bhr;
   int rem[MAX_BHS];

   for (i1=0; i1<MAX_BHS; i1++) rem[i1]=0;

   for (i1=0; i1<MAX_BHS; i1++)
   {
      for (i2=i1+1; i2<MAX_BHS; i2++)
      {
         if (ex_r[i1][0]>0 && ex_r[i2][0]>0)
         {
            is_inside_(&is_1in2,&intr,ex_r[i1],ex_xc[i1],ex_r[i2],ex_xc[i2],&AMRD_dim);
            is_inside_(&is_2in1,&intr,ex_r[i2],ex_xc[i2],ex_r[i1],ex_xc[i1],&AMRD_dim);

//(unless this is disabled, AH[4] is removed at first time step)
            if (is_1in2 || is_2in1)
            {
               if (my_rank==0) printf("\n\nremove_redundant_AH: the following two excised regions "
                  "overlap ... removing the smallest one (is_1in2=%i, is_2in1=%i).\n"
                  "1: ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf,%lf],  "
                  "2: ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf,%lf]\n\n",is_1in2,is_2in1,
                  ex_r[i1][0],ex_r[i1][1],ex_r[i1][2],ex_xc[i1][0],ex_xc[i1][1],ex_xc[i1][2],
                  ex_r[i2][0],ex_r[i2][1],ex_r[i2][2],ex_xc[i2][0],ex_xc[i2][1],ex_xc[i2][2]);
               if (is_2in1) rem[i2]=1; else rem[i1]=1;
            }

            // special case
            if (i2==MAX_BHS-1 && rem[i2]==1)
            { 
               if (my_rank==0) printf("\nremove_redundant_AH: 'unremoving' the tentative encapsulating BH\n");
               rem[i2]=0; 
            }
            if (i2==MAX_BHS-1 && i1!=MAX_BHS-2 && intr) 
            {
               if (my_rank==0) printf("\nremove_redundant_AH: encapsulating BH found, removing intersecting BHs.\n"
                  "1: ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf,%lf],  "
                  "2: ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf,%lf]\n\n",
                  ex_r[i1][0],ex_r[i1][1],ex_r[i1][2],ex_xc[i1][0],ex_xc[i1][1],ex_xc[i1][2],
                  ex_r[i2][0],ex_r[i2][1],ex_r[i2][2],ex_xc[i2][0],ex_xc[i2][1],ex_xc[i2][2]);
               rem[i1]=1;
            }

            if (is_1in2 && is_2in1 && my_rank==0) printf("WARNING ... both is_1in2 && is_2in1\n");
         }
      }
   }

   for (i1=0; i1<MAX_BHS; i1++) if (rem[i1])
   {
       ex_r[i1][0]=0;
       ex_r[i1][1]=0;
       ex_r[i1][2]=0;
       AH_max_iter[i1]=0;
   }
}

//=============================================================================
// the following returns 1 if the test ellipsoid does not intersect any
// existing horizons, *and* is entirely contained within the finest overlapping
// level
//=============================================================================
#define MAX_GRIDS 20
#define EX_MIN_PTS_IN 4
int no_AH_intersect(real ex_r0[3],real ex_xc0[3],int inew)
{
   int i,i1;
   real d;
   real sgh[6*MAX_GRIDS],dx0[3],dt0,bound[6],max_xy;
   int num,L,Lf,contained,ex0_in_exi,ex0_int_exi,ltrace=1;
   int exi_in_ex0,exi_int_ex0;

   for (i=0; i<MAX_BHS; i++)
   {
      if (ex_r[i][0]>0 && i!=inew)
      {
         is_inside_(&ex0_in_exi,&ex0_int_exi,ex_r0,ex_xc0,ex_r[i],ex_xc[i],&AMRD_dim);
         is_inside_(&exi_in_ex0,&exi_int_ex0,ex_r[i],ex_xc[i],ex_r0,ex_xc0,&AMRD_dim);

         if (ex0_int_exi)
         {
            // special case
            if (inew==MAX_BHS-1)
            {
              if (my_rank==0) printf("\nno_AH_intersect: encapsulating BH found, ignoring intersections\n");
              return 1;
            }
            
            if (my_rank==0) printf("\nno_AH_intersect: new BH intersects existing BH\n"
                   "new: ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf.%lf],"
                   "current(%i): ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf,%lf]\n",
                   ex_r0[0],ex_r0[1],ex_r0[2],ex_xc0[0],ex_xc0[1],ex_xc0[2],i,
                   ex_r[i][0],ex_r[i][1],ex_r[i][2],ex_xc[i][0],ex_xc[i][1],ex_xc[i][2]);
            if (my_rank==0) printf("ex0_in_exi=%i,ex0_int_exi=%i,exi_in_ex0=%i,exi_int_ex0=%i\n\n",ex0_in_exi,ex0_int_exi,exi_in_ex0,exi_int_ex0);
         }

         if (ex0_in_exi)
         {
            if (my_rank==0) printf("no_AH_intersect: WARNING ... new BH is *entirely* contained in existing BH\n"
                   "new: ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf.%lf],"
                   "current(%i): ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf,%lf]\n",
                   ex_r0[0],ex_r0[1],ex_r0[2],ex_xc0[0],ex_xc0[1],ex_xc0[2],i,
                   ex_r[i][0],ex_r[i][1],ex_r[i][2],ex_xc[i][0],ex_xc[i][1],ex_xc[i][2]);
            if (my_rank==0) printf("ex0_in_exi=%i,ex0_int_exi=%i,exi_in_ex0=%i,exi_int_ex0=%i\n\n",ex0_in_exi,ex0_int_exi,exi_in_ex0,exi_int_ex0);
            return 0;
         }

         if (ex0_int_exi && !exi_in_ex0) return 0;
      }
   }

   return 1;
}

//=============================================================================
// Routines required by amrd:
//=============================================================================

//=============================================================================
// Returns 0 to use default mechanism, or is expected to calculate
// the correct initial hierarchy and return 1:
//=============================================================================
int CFE4D_id(void)
{
   return 0;
}

//=============================================================================
// Sets custom parameters, variables, etc. Split up into two segments,
// one called before the pamr context is initialized and standard
// parameters are read, and the other afterwards
//=============================================================================
void CFE4D_var_pre_init(char *pfile)
{
   AMRD_echo_params=1;
   AMRD_int_param(pfile,"echo_params",&AMRD_echo_params,1);

   cp_version=ADS5D_CP_VERSION;
   AMRD_int_param(pfile,"cp_version",&cp_version,1);
   
   AdS_L=1.0; AMRD_real_param(pfile,"AdS_L",&AdS_L,1);

   regtype=1; AMRD_int_param(pfile,"regtype",&regtype,1);
   stype=1; AMRD_int_param(pfile,"stype",&stype,1);
   interptype=2; AMRD_int_param(pfile,"interptype",&interptype,1);
   i_shift=0; AMRD_int_param(pfile,"i_shift",&i_shift,1);

   return;
}

void CFE4D_var_post_init(char *pfile)
{
   int i,j,ind;
   char buf[64];
   real rmin,deltar;

   if (my_rank==0)
   {
      printf("===================================================================\n");
      printf("Reading CFE4D parameters:\n\n");
   }

   phi1_amp_1=phi1_B_1=phi1_r0_1=phi1_x0_1[0]=phi1_x0_1[1]=phi1_ecc_1[0]=phi1_ecc_1[1]=0;
   phi1_amp_2=phi1_B_2=phi1_r0_1=phi1_x0_2[0]=phi1_x0_2[1]=phi1_ecc_2[0]=phi1_ecc_2[1]=0;

   AMRD_real_param(pfile,"phi1_amp_1",&phi1_amp_1,1);
   AMRD_real_param(pfile,"phi1_B_1",&phi1_B_1,1);
   AMRD_real_param(pfile,"phi1_r0_1",&phi1_r0_1,1);
   AMRD_real_param(pfile,"phi1_delta_1",&phi1_delta_1,1);
   AMRD_real_param(pfile,"phi1_x0_1",phi1_x0_1,AMRD_dim);
   AMRD_real_param(pfile,"phi1_ecc_1",phi1_ecc_1,AMRD_dim);

   AMRD_real_param(pfile,"phi1_amp_2",&phi1_amp_2,1);
   AMRD_real_param(pfile,"phi1_B_2",&phi1_B_2,1);
   AMRD_real_param(pfile,"phi1_r0_2",&phi1_r0_2,1);
   AMRD_real_param(pfile,"phi1_delta_2",&phi1_delta_2,1);
   AMRD_real_param(pfile,"phi1_x0_2",phi1_x0_2,AMRD_dim);
   AMRD_real_param(pfile,"phi1_ecc_2",phi1_ecc_2,AMRD_dim);

   kappa_cd=0; AMRD_real_param(pfile,"kappa_cd",&kappa_cd,1);
   rho_cd=0; AMRD_real_param(pfile,"rho_cd",&rho_cd,1);

   diss_eps_k=0; AMRD_real_param(pfile,"diss_eps_k",&diss_eps_k,1);
   diss_eps_y_cutoff=1; AMRD_real_param(pfile,"diss_eps_y_cutoff",&diss_eps_y_cutoff,1);
   diss_kmax=0; AMRD_int_param(pfile,"diss_kmax",&diss_kmax,1);
   diss_bdy_k=0; AMRD_int_param(pfile,"diss_bdy_k",&diss_bdy_k,1);
   diss_all_past_k=0; AMRD_int_param(pfile,"diss_all_past_k",&diss_all_past_k,1);
   diss_eps_k_cutoff_n=0; AMRD_int_param(pfile,"diss_eps_k_cutoff_n",&diss_eps_k_cutoff_n,1);
   diss_all=1; AMRD_int_param(pfile,"diss_all",&diss_all,1);

   background=0; AMRD_int_param(pfile,"background",&background,1);
   skip_constraints=0; AMRD_int_param(pfile,"skip_constraints",&skip_constraints,1);
   output_ires=0; AMRD_int_param(pfile,"output_ires",&output_ires,1);
   output_quasiset=0; AMRD_int_param(pfile,"output_quasiset",&output_quasiset,1);

   harmonize=0; AMRD_int_param(pfile,"harmonize",&harmonize,1);

   gauge_t=0; AMRD_int_param(pfile,"gauge_t",&gauge_t,1);
   c1_t=0; AMRD_real_param(pfile,"c1_t",&c1_t,1);
   c2_t=0; AMRD_real_param(pfile,"c2_t",&c2_t,1);
   c3_t=0; AMRD_real_param(pfile,"c3_t",&c3_t,1);
   rho1_t=0; AMRD_real_param(pfile,"rho1_t",&rho1_t,1);
   rho2_t=0; AMRD_real_param(pfile,"rho2_t",&rho2_t,1);
   rho3_t=1; AMRD_real_param(pfile,"rho3_t",&rho3_t,1);
   rho4_t=1; AMRD_real_param(pfile,"rho4_t",&rho4_t,1);
   xi1_t=0; AMRD_real_param(pfile,"xi1_t",&xi1_t,1);
   xi2_t=0; AMRD_real_param(pfile,"xi2_t",&xi2_t,1);
   cbulk_t=2; AMRD_real_param(pfile,"cbulk_t",&cbulk_t,1);

   gauge_i=0; AMRD_int_param(pfile,"gauge_i",&gauge_i,1);
   c1_i=0; AMRD_real_param(pfile,"c1_i",&c1_i,1);
   c2_i=0; AMRD_real_param(pfile,"c2_i",&c2_i,1);
   c3_i=0; AMRD_real_param(pfile,"c3_i",&c3_i,1);
   rho1_i=0; AMRD_real_param(pfile,"rho1_i",&rho1_i,1);
   rho2_i=0; AMRD_real_param(pfile,"rho2_i",&rho2_i,1);
   rho3_i=1; AMRD_real_param(pfile,"rho3_i",&rho3_i,1);
   rho4_i=1; AMRD_real_param(pfile,"rho4_i",&rho4_i,1);
   xi1_i=0; AMRD_real_param(pfile,"xi1_i",&xi1_i,1);
   xi2_i=0; AMRD_real_param(pfile,"xi2_i",&xi2_i,1);
   cbulk_i=2; AMRD_real_param(pfile,"cbulk_i",&cbulk_i,1);

   rhoa=1; AMRD_real_param(pfile,"rhoa",&rhoa,1);
   rhob=1; AMRD_real_param(pfile,"rhob",&rhob,1);

   ex_reset_rbuf=0; AMRD_int_param(pfile,"ex_reset_rbuf",&ex_reset_rbuf,1);

   // set fraction, 1-ex_rbuf, of AH radius to be excised
   for (j=0; j<MAX_BHS; j++)
   {
      if (!AMRD_cp_restart || ex_reset_rbuf || ex_r[j][0]<=0)
      {
         if (j==0) { if (!AMRD_cp_restart) ex_rbuf[j]=0; sprintf(buf,"ex_rbuf"); }
         else { if (!AMRD_cp_restart) ex_rbuf[j]=ex_rbuf[0]; sprintf(buf,"ex_rbuf_%i",j+1); }
         AMRD_real_param(pfile,buf,&ex_rbuf[j],1);
         if (ex_rbuf[j]<0 || ex_rbuf[j]>1 ) printf("WARNING ... ex_rbuf[%i]=%lf is outside of standard range\n",j,ex_rbuf[j]);
      }

      if (!AMRD_cp_restart)
      {
         ex_xc[j][0]=ex_xc[j][1]=0;
         ex_r[j][0]=ex_r[j][1]=0;
      }
   }

   // set AH parameters
   int AH_count[MAX_BHS],found_AH[MAX_BHS],freq0[MAX_BHS];

   use_AH_new_smooth=0; AMRD_int_param(pfile,"use_AH_new_smooth",&use_AH_new_smooth,1);
   use_AH_new=0; AMRD_int_param(pfile,"use_AH_new",&use_AH_new,1);

   for (j=0; j<MAX_BHS; j++)
   {
      // because the AH shape is saved, we can't currently change that upon
      // a restart (unless we haven't yet found one)
      if (!AMRD_cp_restart || !found_AH[j])
      {
         if (AMRD_cp_restart) free(AH_R[j]);
         if (j==0) { AH_Nchi[j]=17; sprintf(buf,"AH_Nchi"); } 
         else { AH_Nchi[j]=AH_Nchi[0]; sprintf(buf,"AH_Nchi_%i",j+1); }
         AMRD_int_param(pfile,buf,&AH_Nchi[j],1);

         AH_Nphi[j]=1;

         if (!( ((AH_Nchi[j]-1)/2)==(AH_Nchi[j]/2) ) ||
             !( ((AH_Nphi[j]-1)/2)==(AH_Nphi[j]/2) ))
         { 
            printf("\n\n\n\n\n WARNING: SMOOTHING AND REGULARIZATION IN AH ROUTINES ASSUME\n"
                   " AN ODD NUMBER OF POINTS IN AH_Nchi and AH_Nphi\n\n\n\n\n");
         }
      }

      if (j==0) { AH_find_best_fit[j]=0; sprintf(buf,"AH_find_best_fit"); }
      else { AH_find_best_fit[j]=AH_find_best_fit[0]; sprintf(buf,"AH_find_best_fit_%i",j+1); }
      AMRD_int_param(pfile,buf,&AH_find_best_fit[j],1);

      if (j==0) { AH_Lmin[j]=2; sprintf(buf,"AH_Lmin"); }
      else { AH_Lmin[j]=AH_Lmin[0]; sprintf(buf,"AH_Lmin_%i",j+1); }
      AMRD_int_param(pfile,buf,&AH_Lmin[j],1);

      if (j==0) { AH_Lmax[j]=100; sprintf(buf,"AH_Lmax"); }
      else { AH_Lmax[j]=AH_Lmax[0]; sprintf(buf,"AH_Lmax_%i",j+1); }
      AMRD_int_param(pfile,buf,&AH_Lmax[j],1);

      if (j==0) { AH_max_iter[j]=0; sprintf(buf,"AH_max_iter"); }
      else { AH_max_iter[j]=AH_max_iter[0]; sprintf(buf,"AH_max_iter_%i",j+1); }
      AMRD_int_param(pfile,buf,&AH_max_iter[j],1);

      if (j==0) { AH_freq[j]=1; sprintf(buf,"AH_freq"); }
      else { AH_freq[j]=AH_freq[0]; sprintf(buf,"AH_freq_%i",j+1); }
      AMRD_int_param(pfile,buf,&AH_freq[j],1);

      if (j==0) { AH_freq_aft[j]=1; sprintf(buf,"AH_freq_aft"); }
      else { AH_freq_aft[j]=AH_freq_aft[0]; sprintf(buf,"AH_freq_aft_%i",j+1); }
      AMRD_int_param(pfile,buf,&AH_freq_aft[j],1);

      if (AMRD_cp_restart)
      {
         if (found_AH[j]) freq0[j]=AH_freq_aft[j];
         else freq0[j]=AH_freq[j];
      }

      if (j==0) { AH_rsteps[j]=1; sprintf(buf,"AH_rsteps"); }
      else { AH_rsteps[j]=AH_rsteps[0]; sprintf(buf,"AH_rsteps_%i",j+1); }
      AMRD_int_param(pfile,buf,&AH_rsteps[j],1);

      if (j==0) { AH_maxinc[j]=10; sprintf(buf,"AH_maxinc"); }
      else { AH_maxinc[j]=AH_maxinc[0]; sprintf(buf,"AH_maxinc_%i",j+1); }
      AMRD_int_param(pfile,buf,&AH_maxinc[j],1);

      if (j==0) { AH_tol[j]=1e-2; sprintf(buf,"AH_tol"); }
      else { AH_tol[j]=AH_tol[0]; sprintf(buf,"AH_tol_%i",j+1); }
      AMRD_real_param(pfile,buf,&AH_tol[j],1);

      AH_tol_aft[j]=AH_tol[j];
      if (j==0) sprintf(buf,"AH_tol_aft"); else sprintf(buf,"AH_tol_aft_%i",j+1);
      AMRD_real_param(pfile,buf,&AH_tol_aft[j],1);

      if (j==0) { AH_max_tol_inc[j]=1; sprintf(buf,"AH_max_tol_inc"); }
      else { AH_max_tol_inc[j]=AH_max_tol_inc[0]; sprintf(buf,"AH_max_tol_inc_%i",j+1); }
      AMRD_real_param(pfile,buf,&AH_max_tol_inc[j],1);

      if (j==0) { AH_tol_scale[j]=1; sprintf(buf,"AH_tol_scale"); }
      else { AH_tol_scale[j]=AH_tol_scale[0]; sprintf(buf,"AH_tol_scale_%i",j+1); }
      AMRD_real_param(pfile,buf,&AH_tol_scale[j],1);

      if (j==0) { AH_omt_scale[j]=1.1; sprintf(buf,"AH_omt_scale"); }
      else { AH_omt_scale[j]=AH_omt_scale[0]; sprintf(buf,"AH_omt_scale_%i",j+1); }
      AMRD_real_param(pfile,buf,&AH_omt_scale[j],1);

      if (j==0) { AH_reset_scale[j]=0; sprintf(buf,"AH_reset_scale"); }
      else { AH_reset_scale[j]=AH_reset_scale[0]; sprintf(buf,"AH_reset_scale_%i",j+1); }
      AMRD_real_param(pfile,buf,&AH_reset_scale[j],1);

      if (j==0) { AH_lambda_min[j]=0.1; sprintf(buf,"AH_lambda_min"); }
      else { AH_lambda_min[j]=AH_lambda_min[0]; sprintf(buf,"AH_lambda_min_%i",j+1); }
      AMRD_real_param(pfile,buf,&AH_lambda_min[j],1);

      if (j==0) { AH_lambda[j]=0.1; sprintf(buf,"AH_lambda"); }
      else { AH_lambda[j]=AH_lambda[0]; sprintf(buf,"AH_lambda_%i",j+1); }
      AMRD_real_param(pfile,buf,&AH_lambda[j],1);

      if (j==0) { AH_r0[j]=0.1; sprintf(buf,"AH_r0"); }
      else { AH_r0[j]=AH_r0[0]; sprintf(buf,"AH_r0_%i",j+1); }
      AMRD_real_param(pfile,buf,&AH_r0[j],1);

      if (j==0) { AH_r1[j]=0.2; sprintf(buf,"AH_r1"); }
      else { AH_r1[j]=AH_r1[0]; sprintf(buf,"AH_r1_%i",j+1); }
      AMRD_real_param(pfile,buf,&AH_r1[j],1);

      if (j==0) { AH_tmin[j]=0.0; sprintf(buf,"AH_tmin"); }
      else { AH_tmin[j]=AH_tmin[0]; sprintf(buf,"AH_tmin_%i",j+1); }
      AMRD_real_param(pfile,buf,&AH_tmin[j],1);

      if (!AMRD_cp_restart || !found_AH[j])
      {
         AH_xc[j][0]=AH_xc[j][1]=0;
         if (j==0) sprintf(buf,"AH_xc");
         else sprintf(buf,"AH_xc_%i",j+1);
         AMRD_real_param(pfile,buf,AH_xc[j],AMRD_dim);
      }

      if (j==0) { AH_eps[j]=0.0; sprintf(buf,"AH_eps"); }
      else { AH_eps[j]=AH_eps[0]; sprintf(buf,"AH_eps_%i",j+1); }
      AMRD_real_param(pfile,buf,&AH_eps[j],1);

      if (AH_Nchi[j]<5 || (AMRD_dim==3 && AH_Nphi[j]<5)) AMRD_stop("error ... AH_Nchi<5 || AH_Nphi<5 \n","");
      if (AH_rsteps[j]<1) AMRD_stop("error ... AH_rsteps<1 \n","");

      AH_theta[j]=(real *)malloc(sizeof(real)*AH_Nchi[j]*AH_Nphi[j]);
      if (!AMRD_cp_restart) AH_R[j]=(real *)malloc(sizeof(real)*AH_Nchi[j]*AH_Nphi[j]);
      AH_w1[j]=(real *)malloc(sizeof(real)*AH_Nchi[j]*AH_Nphi[j]);
      AH_w2[j]=(real *)malloc(sizeof(real)*AH_Nchi[j]*AH_Nphi[j]);
      AH_w3[j]=(real *)malloc(sizeof(real)*AH_Nchi[j]*AH_Nphi[j]);
      AH_theta_ads[j]=(real *)malloc(sizeof(real)*AH_Nchi[j]*AH_Nphi[j]);
      AH_own[j]=(int *)malloc(sizeof(int)*AH_Nchi[j]*AH_Nphi[j]);
      AH_lev[j]=(int *)malloc(sizeof(int)*AH_Nchi[j]*AH_Nphi[j]);

   }

   ex_max_repop=0; AMRD_int_param(pfile,"ex_max_repop",&ex_max_repop,1);
   ex_repop_buf=1; AMRD_int_param(pfile,"ex_repop_buf",&ex_repop_buf,1);
   ex_repop_io=1; AMRD_int_param(pfile,"ex_repop_io",&ex_repop_io,1);
   if (abs(ex_repop_io)<1 || abs(ex_repop_io)>4) AMRD_stop("invalid |ex_repop_io| ... must be 1,2,3 or 4","");

   ief_bh_r0=0; AMRD_real_param(pfile,"ief_bh_r0",&ief_bh_r0,1);


   // later have a proper routine to convert from AH to excision parameters.
   // rh is radius in uncompcatified coordicates, rhoh in compactified,
   // ief_bh_r0 is BH radius parameter, ex_r is excision radius
   int l,ah_finder_is_off=1; 
   for (l=0; l<MAX_BHS; l++) {if (AH_max_iter[l]!=0) ah_finder_is_off=0;}
   if (ah_finder_is_off) 
   {
     real rh,rhoh,mh;
     rh=AdS_L*sqrt(sqrt(1+4*ief_bh_r0*ief_bh_r0/AdS_L/AdS_L)/2-0.5);
     mh=3*M_PI/8*ief_bh_r0*ief_bh_r0;
     rhoh=(sqrt(1+rh*rh)-1)/rh;
     ex_r[0][0]=ex_r[0][1]=rhoh*(1-ex_rbuf[0]);
     if (my_rank==0) printf("\nBH initial data\n"
                            "r0/L=%lf, rh/L=%lf, mass M = 3*PI/8*r0^2 = 3*PI/8*rh^2*(1+rh^2) = %lf\n"
                            "Initial BH radius=%lf, (%lf in compactified (code) coords)\n"
                            "Initial excision radius=%lf\n\n",ief_bh_r0/AdS_L,rh/AdS_L,mh,rh,rhoh,ex_r[0][0]);
   }

   if (AMRD_do_ex==0) AMRD_stop("require excision to be on","");

   PAMR_excision_on("chr",&CFE4D_fill_ex_mask,AMRD_ex,1);

   if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Sets all variables to their 'zero' values:
//=============================================================================
void CFE4D_AMRH_var_clear(void)
{
   int i,j,ind;

   ldptr();

   zero_f(phi1_n); 

   init_ghb_ads_(gb_tt_n,gb_tx_n,gb_ty_n,gb_xx_n,gb_xy_n,gb_yy_n,psi_n,
                 Hb_t_n,Hb_x_n,Hb_y_n,&AdS_L,x,y,chr,&AMRD_ex,&Nx,&Ny,&regtype);

   return;
}

//=============================================================================
// Initial data for free fields: (at tn=2) ... following vars also used in 
// t0_cnst_data
//=============================================================================
void CFE4D_free_data(void)
{
   int i;

   ldptr();

   CFE4D_AMRH_var_clear(); // constrained variables are set post-MG

   zero_f(phi1_t_n); // holds initial time derivatives for ID

//TEMPORARY//
   gauss2d_(eb_xx_n,&phi1_amp_1,&phi1_B_1,&phi1_r0_1,&phi1_delta_1,&phi1_x0_1[0],&phi1_x0_1[1],
            &phi1_ecc_1[0],&phi1_ecc_1[1],&AdS_L,x,y,&Nx,&Ny,&stype);

   gauss2d_(phi1_n,&phi1_amp_1,&phi1_B_1,&phi1_r0_1,&phi1_delta_1,&phi1_x0_1[0],&phi1_x0_1[1],
            &phi1_ecc_1[0],&phi1_ecc_1[1],&AdS_L,x,y,&Nx,&Ny,&stype);

   gauss2d_(w1,&phi1_amp_2,&phi1_B_2,&phi1_r0_2,&phi1_delta_2,&phi1_x0_2[0],&phi1_x0_2[1],
            &phi1_ecc_2[0],&phi1_ecc_2[1],&AdS_L,x,y,&Nx,&Ny,&stype);

   for (i=0; i<size; i++) phi1_n[i]+=w1[i]; 

   return;
}  

//=============================================================================
// Initialize any "elliptic_vars_t0" post construction of MGH, but before
// the start of vcycling.
//=============================================================================
void CFE4D_elliptic_vars_t0_init(void)
{
   // initializes dt, dx, dy
   ldptr_mg();

   // initializes zeta conformal factor to background AdS value zeta=1
   const_f(zeta,1);

}

//=============================================================================
// Initial constraint data --- called after each MG iteration.
//
// Here we also initialize past time level information if 
// AMRD_id_pl_method==3
//
// NOTE: not cleaning up memory allocation after reading in square black hole 
//       data
//
// NOTE: np1,n,nm1 variables are allocated only at the top level of the MG hierarchy,
//       so do an if(f_nm1){...}, for example, to make sure we are at the top level
//=============================================================================
void CFE4D_t0_cnst_data(void)
{
   int i,j,ind;

   ldptr_mg();

   // initialize time derivatives of gbars,hbars
   if (gb_xx_nm1)
   {
     init_ghbdot_(gb_tt_n,gb_tx_n,gb_ty_n,gb_xx_n,gb_xy_n,gb_yy_n,psi_n,
                  gb_tt_t_n,gb_tx_t_n,gb_ty_t_n,gb_xx_t_n,gb_xy_t_n,gb_yy_t_n,psi_t_n,
                  Hb_t_n,Hb_x_n,Hb_y_n,Hb_t_t_n,Hb_x_t_n,Hb_y_t_n,
                  &AdS_L,phys_bdy,x,y,&dt,chr,&AMRD_ex,&Nx,&Ny,&regtype);
   }

   // initialize gbars
   if ((background || skip_constraints) && ief_bh_r0==0)
   {
     init_ghb_ads_(gb_tt,gb_tx,gb_ty,gb_xx,gb_xy,gb_yy,psi,
                   Hb_t,Hb_x,Hb_y,&AdS_L,x,y,chr_mg,&AMRD_ex,&Nx,&Ny,&regtype);
   }
   else if (background || skip_constraints)
   {
     if (gb_xx_nm1) //"np1,n,nm1" variables only allocated on finest MG level
     {
       init_ads5d_bh_(&ief_bh_r0,&AdS_L,gb_tt,gb_tx,gb_ty,gb_xx,gb_xy,gb_yy,psi,
                      gb_tt_t_n,gb_tx_t_n,gb_ty_t_n,gb_xx_t_n,gb_xy_t_n,gb_yy_t_n,psi_t_n,
                      Hb_t,Hb_x,Hb_y,Hb_t_t_n,Hb_x_t_n,Hb_y_t_n,
                      phys_bdy,x,y,&dt,chr_mg,&AMRD_ex,&Nx,&Ny,&regtype);
     }
   }
   else
   {
     init_ghb_(zeta,
               gb_tt,gb_tx,gb_ty,gb_xx,gb_xy,gb_yy,psi,
               Hb_t,Hb_x,Hb_y,
               &AdS_L,mask_mg,phys_bdy,x,y,chr_mg,&AMRD_ex,&Nx,&Ny,&regtype,&rhoa,&rhob);
   }

   // initialize hbars and nm1, np1 time levels
   if (AMRD_id_pl_method==3 && gb_xx_nm1)
   {

//TEMPORARY//
//     init_hb_(gb_tt_np1,gb_tt_n,gb_tt_nm1,
//              gb_tx_np1,gb_tx_n,gb_tx_nm1,
//              gb_ty_np1,gb_ty_n,gb_ty_nm1,
//              gb_xx_np1,gb_xx_n,gb_xx_nm1,
//              gb_xy_np1,gb_xy_n,gb_xy_nm1,
//              gb_yy_np1,gb_yy_n,gb_yy_nm1,
//              psi_np1,psi_n,psi_nm1,
//              Hb_t_n,Hb_x_n,Hb_y_n,
//              &AdS_L,phys_bdy,x,y,&dt,chr,&AMRD_ex,&Nx,&Ny,&regtype);
//
//     init_nm1_(gb_tt_np1,gb_tt_n,gb_tt_nm1,gb_tt_t_n,
//               gb_tx_np1,gb_tx_n,gb_tx_nm1,gb_tx_t_n,
//               gb_ty_np1,gb_ty_n,gb_ty_nm1,gb_ty_t_n,
//               gb_xx_np1,gb_xx_n,gb_xx_nm1,gb_xx_t_n,
//               gb_xy_np1,gb_xy_n,gb_xy_nm1,gb_xy_t_n,
//               gb_yy_np1,gb_yy_n,gb_yy_nm1,gb_yy_t_n,
//               psi_np1,psi_n,psi_nm1,psi_t_n,
//               Hb_t_np1,Hb_t_n,Hb_t_nm1,Hb_t_t_n,
//               Hb_x_np1,Hb_x_n,Hb_x_nm1,Hb_x_t_n,
//               Hb_y_np1,Hb_y_n,Hb_y_nm1,Hb_y_t_n,
//               phi1_np1,phi1_n,phi1_nm1,phi1_t_n,tfunction,
//               &AdS_L,phys_bdy,x,y,&dt,chr,&AMRD_ex,&Nx,&Ny,&regtype);

     for (i=0; i<size; i++)
     {
       gb_tt_np1[i]=gb_tt_nm1[i]=gb_tt[i];
       gb_tx_np1[i]=gb_tx_nm1[i]=gb_tx[i];
       gb_ty_np1[i]=gb_ty_nm1[i]=gb_ty[i];
       gb_xx_np1[i]=gb_xx_nm1[i]=gb_xx[i];
       gb_xy_np1[i]=gb_xy_nm1[i]=gb_xy[i];
       gb_yy_np1[i]=gb_yy_nm1[i]=gb_yy[i];
       psi_np1[i]=psi_nm1[i]=psi[i];
     }

   }

   // store initial source functions 
   if (gb_xx_nm1)
   {
     for (i=0; i<size; i++)
     {
       Hb_t_0[i]=Hb_t[i];
       Hb_x_0[i]=Hb_x[i];
       Hb_y_0[i]=Hb_y[i];
     }
   }

   // fill test functions with conformal factor diagnostic profiles
   for (i=0; i<Nx; i++)
   {
      for (j=0; j<Ny; j++)
      {
         ind=i+j*Nx;
         test1[ind]=zeta[ind];
         test2[ind]=zeta_lop[ind];
         test3[ind]=zeta_rhs[ind];
         test4[ind]=zeta_res[ind];
      } 
   }
  
   return;
}

//=============================================================================
// Calculations prior to saving info to disk.
//
// NOTE: at this point, the time sequence is: n,nm1,np1
//=============================================================================
void CFE4D_pre_io_calc(void)
{
   ldptr();

   int i,j,ind;
   int j_shift;
   real ct,rho;

   dx=x[1]-x[0]; dy=y[1]-y[0];
   ct=PAMR_get_time(g_L);

   if (regtype==6) j_shift=3;
   if (regtype==5 || regtype==4 || regtype==3) j_shift=2;
   if (regtype==2 || regtype==1) j_shift=1;

   // output independent residual
   if (output_ires)
   {

      // compute independent residuals of the CFE4D system
      if (ct!=0)
      {
      //(NOTE: for t>t0, have cycled time sequence np1,n,nm1 to time sequence n,nm1,np1,
      // so here, time level n is the most advanced time level)
      ires_(efe_all_ires,
         efe_tt_ires,efe_tx_ires,efe_ty_ires,
         efe_xx_ires,efe_xy_ires,efe_yy_ires,
         efe_psi_ires,
         gb_tt_n,gb_tt_nm1,gb_tt_np1,
         gb_tx_n,gb_tx_nm1,gb_tx_np1,
         gb_ty_n,gb_ty_nm1,gb_ty_np1,
         gb_xx_n,gb_xx_nm1,gb_xx_np1,
         gb_xy_n,gb_xy_nm1,gb_xy_np1,
         gb_yy_n,gb_yy_nm1,gb_yy_np1,
         psi_n,psi_nm1,psi_np1,
         Hb_t_n,Hb_t_nm1,Hb_t_np1,
         Hb_x_n,Hb_x_nm1,Hb_x_np1,
         Hb_y_n,Hb_y_nm1,Hb_y_np1,
         phi1_n,phi1_nm1,phi1_np1,
         x,y,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,phys_bdy,ghost_width);
      }
      else
      {
      //(NOTE: for t=t0, have *not* cycled time sequence, so still np1,n,nm1,
      // so here, time level np1 is the most advanced time level)
      ires_(efe_all_ires,
         efe_tt_ires,efe_tx_ires,efe_ty_ires,
         efe_xx_ires,efe_xy_ires,efe_yy_ires,
         efe_psi_ires,
         gb_tt_np1,gb_tt_n,gb_tt_nm1,
         gb_tx_np1,gb_tx_n,gb_tx_nm1,
         gb_ty_np1,gb_ty_n,gb_ty_nm1,
         gb_xx_np1,gb_xx_n,gb_xx_nm1,
         gb_xy_np1,gb_xy_n,gb_xy_nm1,
         gb_yy_np1,gb_yy_n,gb_yy_nm1,
         psi_np1,psi_n,psi_nm1,
         Hb_t_np1,Hb_t_n,Hb_t_nm1,
         Hb_x_np1,Hb_x_n,Hb_x_nm1,
         Hb_y_np1,Hb_y_n,Hb_y_nm1,
         phi1_np1,phi1_n,phi1_nm1,
         x,y,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,phys_bdy,ghost_width);
      }

      // fill in independent residual evaluator test functions
      for (i=0; i<Nx; i++)
      {
         for (j=0; j<Ny; j++)
         {
            ind=i+j*Nx;
            rho=sqrt(x[i]*x[i]+y[j]*y[j]);

            // excise rho=1-1.5*dx pts (pure AdS diverges at rho=1, so cannot use these pts in difference stencils) 
            if (chr[ind]==AMRD_ex || 1-rho<1.5*dx_Lc)
            {
              iresall[ind]=0;  
              irestt[ind]=0;  
              irestx[ind]=0;  
              iresty[ind]=0;  
              iresxx[ind]=0;  
              iresxy[ind]=0;  
              iresyy[ind]=0;  
              irespsi[ind]=0;  
            }
            else
            {
              iresall[ind]=efe_all_ires[ind];
              irestt[ind]=efe_tt_ires[ind];
              irestx[ind]=efe_tx_ires[ind];
              iresty[ind]=efe_ty_ires[ind];
              iresxx[ind]=efe_xx_ires[ind];
              iresxy[ind]=efe_xy_ires[ind];
              iresyy[ind]=efe_yy_ires[ind];
              irespsi[ind]=efe_psi_ires[ind];
            }
         } 
      }

   }

   return;
}

// to be callable from fortran
void check_nan_(real *x, int *is_nan)
{
   if (isnan(*x)) *is_nan=1; else *is_nan=0;
}

//=============================================================================
// Returns some norm of the residual for the evolution variables ... just
// use the value from the most recent evolution step
//=============================================================================
#define LIN_ZERO_BND 1
int lin_zero_bnd_all=1;
real CFE4D_evo_residual(void)
{
   real l2norm=0,l2norm_phi1,l2norm_gb,l2norm_hb_t,l2norm_hb_i;
   int is_nan;

   ldptr();

   if (LIN_ZERO_BND) 
   {
      lin_zero_bnd_res_(phi1_res,phys_bdy,&lin_zero_bnd_all,&Nx,&Ny);
      lin_zero_bnd_res_(gb_res,phys_bdy,&lin_zero_bnd_all,&Nx,&Ny);
      lin_zero_bnd_res_(hb_t_res,phys_bdy,&lin_zero_bnd_all,&Nx,&Ny);
      lin_zero_bnd_res_(hb_i_res,phys_bdy,&lin_zero_bnd_all,&Nx,&Ny);
   }

   l2norm_phi1=norm_l2(phi1_res,mask,chr);
   l2norm_gb=norm_l2(gb_res,mask,chr);
   l2norm_hb_t=norm_l2(hb_t_res,mask,chr);
   l2norm_hb_i=norm_l2(hb_i_res,mask,chr);

   l2norm=l2norm_phi1;

   if (!background) l2norm+=(l2norm_gb+l2norm_hb_t+l2norm_hb_i);

   check_nan_(&l2norm,&is_nan);

//   if (is_nan)
//   {
//      printf("\nl2norm_phi1=%lf, l2norm_gb=%lf, l2norm_hb_t=%lf, l2norm_hb_i=%lf, g_norms[phi1_n_gfn-1]=%lf\n",
//              l2norm_phi1,l2norm_gb,l2norm_hb_t,l2norm_hb_i,g_norms[phi1_n_gfn-1]);
//      printf("[Nx,Ny]=[%i,%i],L=%i\n",Nx,Ny,g_L);
//      AMRD_stop("l2norm is nan ... stopping","");
//      l2norm=0;
//   }

   return l2norm;
}

//=============================================================================
//---lower frequency KO dissipation filter, called from CFE4D_evolve()
//   after pointers set up
//=============================================================================
void apply_diss_eps_k(void)
{
   int pbt_even[4]={PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_EVEN,PAMR_UNKNOWN};
   int pbt_unknown[4]={PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_UNKNOWN};
   int pbt_odd[4]={PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_ODD,PAMR_UNKNOWN};
   int even=PAMR_EVEN,odd=PAMR_ODD,i,j,ind,Nz;

   // define eps as an array, and dissipate each dimension one at a time 
   int do_ex=-1,ind_sweeps=1;

   real rho;
   dx=x[1]-x[0]; dy=y[1]-y[0];

   // using w2 for eps array, with new diss_eps_y_cutoff flag (off for diss_eps_y_cutoff=1) 
   for (i=0; i<Nx; i++)
      for (j=0; j<Ny; j++)
      {
         ind=i+j*Nx;
         rho=sqrt(x[i]*x[i]+y[j]*y[j]);
         if (rho>=1 || y[j]>diss_eps_y_cutoff) 
         {
           w2[ind]=0; 
         }
         else 
         {
           w2[ind]=pow(rho,diss_eps_k_cutoff_n)*diss_eps_k;
         }
      } 

   // only 2D Kreiss-Oliger dissipation
   Nz=1;

   // define effective excision mask w3 for dmdiss3d
   for (i=0; i<Nx; i++)
   {
      for (j=0; j<Ny; j++)
      {
         ind=i+j*Nx;
         rho=sqrt(x[i]*x[i]+y[j]*y[j]);
         w3[ind]=0;
         if (1-rho<i_shift*dx+dx/2) w3[ind]=AMRD_ex;
      }
   }

   // use effective excision mask w3 instead of excision mask chr with new diss_all flag
   if (diss_all==1)
   {
   dmdiss3d_ex_gen_(gb_tt_n,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
   dmdiss3d_ex_gen_(gb_tx_n,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
   dmdiss3d_ex_gen_(gb_ty_n,w1,w2,&diss_bdy_k,pbt_odd,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
   dmdiss3d_ex_gen_(gb_xx_n,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
   dmdiss3d_ex_gen_(gb_xy_n,w1,w2,&diss_bdy_k,pbt_odd,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
   dmdiss3d_ex_gen_(gb_yy_n,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
   dmdiss3d_ex_gen_(psi_n,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
   dmdiss3d_ex_gen_(Hb_t_n,w1,w2,&diss_bdy_k,pbt_unknown,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
   dmdiss3d_ex_gen_(Hb_x_n,w1,w2,&diss_bdy_k,pbt_unknown,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
   dmdiss3d_ex_gen_(Hb_y_n,w1,w2,&diss_bdy_k,pbt_unknown,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
   dmdiss3d_ex_gen_(phi1_n,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
   }
   else {
   dmdiss3d_ex_gen_(gb_xx_n,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
   }

   if (diss_all_past_k)
   {
      dmdiss3d_ex_gen_(gb_tt_nm1,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
      dmdiss3d_ex_gen_(gb_tx_nm1,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
      dmdiss3d_ex_gen_(gb_ty_nm1,w1,w2,&diss_bdy_k,pbt_odd,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
      dmdiss3d_ex_gen_(gb_xx_nm1,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
      dmdiss3d_ex_gen_(gb_xy_nm1,w1,w2,&diss_bdy_k,pbt_odd,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
      dmdiss3d_ex_gen_(gb_yy_nm1,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
      dmdiss3d_ex_gen_(psi_nm1,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
      dmdiss3d_ex_gen_(Hb_t_nm1,w1,w2,&diss_bdy_k,pbt_unknown,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
      dmdiss3d_ex_gen_(Hb_x_nm1,w1,w2,&diss_bdy_k,pbt_unknown,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
      dmdiss3d_ex_gen_(Hb_y_nm1,w1,w2,&diss_bdy_k,pbt_unknown,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
      dmdiss3d_ex_gen_(phi1_nm1,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
   }

   return;
}

//=============================================================================
// Performs 1 iteration of the evolution equations 
//
// NOTE: at this point, the time sequence is: np1,n,nm1
//=============================================================================
void CFE4D_evolve(int iter)
{
   int i,zero_i=0;
   int ltrace=0;
   real ct,zero=0;

   ldptr();

   ct=PAMR_get_time(g_L);

   if (ltrace) printf("CFE4D_evolve: iter=%i , time=%lf, lev=%i, rank=%i\n",iter,ct,g_L,my_rank);

   zero_f(hb_t_res);
   zero_f(hb_i_res);

   // when kmax nonzero, apply extra dissipation
   if (diss_kmax>=1 && diss_eps_k>0 && iter==1) apply_diss_eps_k();

   if (!background)
   {
      hb_t_evo_(hb_t_res,
                gb_tt_np1,gb_tt_n,gb_tt_nm1,
                gb_tx_np1,gb_tx_n,gb_tx_nm1,
                gb_ty_np1,gb_ty_n,gb_ty_nm1,
                gb_xx_np1,gb_xx_n,gb_xx_nm1,
                gb_xy_np1,gb_xy_n,gb_xy_nm1,
                gb_yy_np1,gb_yy_n,gb_yy_nm1,
                psi_np1,psi_n,psi_nm1,
                Hb_t_np1,Hb_t_n,Hb_t_nm1,
                Hb_x_np1,Hb_x_n,Hb_x_nm1,
                Hb_y_np1,Hb_y_n,Hb_y_nm1,
                phi1_np1,phi1_n,phi1_nm1,
                &AdS_L,x,y,&dt,chr,&AMRD_ex,
                phys_bdy,ghost_width,&Nx,&Ny,
                Hb_t_0,Hb_x_0,Hb_y_0,
                &gauge_t,&ct,&rho1_t,&rho2_t,&rho3_t,&rho4_t,&xi1_t,&xi2_t,
                &c1_t,&c2_t,&c3_t,&cbulk_t);

      hb_i_evo_(hb_i_res,
                gb_tt_np1,gb_tt_n,gb_tt_nm1,
                gb_tx_np1,gb_tx_n,gb_tx_nm1,
                gb_ty_np1,gb_ty_n,gb_ty_nm1,
                gb_xx_np1,gb_xx_n,gb_xx_nm1,
                gb_xy_np1,gb_xy_n,gb_xy_nm1,
                gb_yy_np1,gb_yy_n,gb_yy_nm1,
                psi_np1,psi_n,psi_nm1,
                Hb_t_np1,Hb_t_n,Hb_t_nm1,
                Hb_x_np1,Hb_x_n,Hb_x_nm1,
                Hb_y_np1,Hb_y_n,Hb_y_nm1,
                phi1_np1,phi1_n,phi1_nm1,
                &AdS_L,x,y,&dt,chr,&AMRD_ex,
                phys_bdy,ghost_width,&Nx,&Ny,
                Hb_t_0,Hb_x_0,Hb_y_0,
                &gauge_i,&ct,&rho1_i,&rho2_i,&rho3_i,&rho4_i,&xi1_i,&xi2_i,
                &c1_i,&c2_i,&c3_i,&cbulk_i);

     g_evo_opt_(eb_res,gb_res,phi1_res,cl_res,
                eb_xx_np1,eb_xx_n,eb_xx_nm1,
                eb_xy_np1,eb_xy_n,eb_xy_nm1,
                eb_xz_np1,eb_xz_n,eb_xz_nm1,
                eb_yy_np1,eb_yy_n,eb_yy_nm1,
                eb_yz_np1,eb_yz_n,eb_yz_nm1,
                eb_zz_np1,eb_zz_n,eb_zz_nm1,
                gb_tt_np1,gb_tt_n,gb_tt_nm1,
                gb_tx_np1,gb_tx_n,gb_tx_nm1,
                gb_ty_np1,gb_ty_n,gb_ty_nm1,
                gb_xx_np1,gb_xx_n,gb_xx_nm1,
                gb_xy_np1,gb_xy_n,gb_xy_nm1,     
                gb_yy_np1,gb_yy_n,gb_yy_nm1,
                psi_np1,psi_n,psi_nm1,
                Hb_t_np1,Hb_t_n,Hb_t_nm1,
                Hb_x_np1,Hb_x_n,Hb_x_nm1,
                Hb_y_np1,Hb_y_n,Hb_y_nm1,
                phi1_np1,phi1_n,phi1_nm1,
                &AdS_L,x,y,&dt,chr,&AMRD_ex,
                phys_bdy,ghost_width,&Nx,&Ny,
                &background,&kappa_cd,&rho_cd,
                &interptype,&i_shift,&regtype,
                &diss_kmax,tfunction);
   }

   return;
}

//=============================================================================
// sets excision mask (NO ITERATOR, SO DON'T LOAD POINTERS!!!)
// 
// outside rho=1 grid is also excised
//=============================================================================
void CFE4D_fill_ex_mask(real *mask, int dim, int *shape, real *bbox, real excised)
{
   int i,j,ind,l;
   real x,y,dx,dy,rho,xp,yp,ex_r_xp,ex_r_yp,r;

   dx=(bbox[1]-bbox[0])/(shape[0]-1);
   dy=(bbox[3]-bbox[2])/(shape[1]-1);

   for (i=0; i<shape[0]; i++)
   {
      x=bbox[0]+i*dx;
      for (j=0; j<shape[1]; j++)
      {
         y=bbox[2]+j*dy;
         rho=sqrt(x*x+y*y);
         ind=i+j*shape[0];

         if (rho>=(1-dx_Lc/2)) 
         {
            //excise outside rho>=1-dx
            mask[ind]=excised; 
         }
         else 
         {
            mask[ind]=excised-1;
            //when MAX_BHS.ne.0, excise inside a fraction ex_r of the horizon
            for (l=0; l<MAX_BHS; l++)
            {
               if (ex_r[l][0]>0)
               {
                 xp=(x-ex_xc[l][0]);
                 yp=(y-ex_xc[l][1]);
                 ex_r_xp=(ex_r[l][0]*(1-ex_rbuf[l]));
                 ex_r_yp=(ex_r[l][1]*(1-ex_rbuf[l]));
                 if ((r=sqrt(xp*xp/ex_r_xp/ex_r_xp+yp*yp/ex_r_yp/ex_r_yp))<1) mask[ind]=excised;
               }
            }
         }
      }
   }
}

//=============================================================================
//=============================================================================
void CFE4D_fill_bh_bboxes(real *bbox, int *num, int max_num)
{
   if (max_num<MAX_BHS) AMRD_stop("CFE4D_fill_bh_bboxes: error max_num too small\n","");

   *num=0;
}

//=============================================================================
// The following routine searches for AH's, manages excision, 
// and if t==0, determines past time level information using evolve
//
// NOTE: at this point, the time sequence is: n,nm1,np1 (unless t=0)
//=============================================================================
#define AH_RESET_AFTER_FAIL 0
int AH_count[MAX_BHS],found_AH[MAX_BHS],freq0[MAX_BHS];
int found_count_AH[MAX_BHS];
int pre_tstep_global_first=1,valid;

void CFE4D_pre_tstep(int L)
{
   char name[256];
   int AH[MAX_BHS];
   int AH_shape[1],got_an_AH,do_reinit_ex,do_repop;
   real M,J,c_equat,c_polar;
   real AH_bbox[2],AH_min_resid0,AH_min_resid1,AH_min_resid2;
   real ex_r0[2],ex_xc0[2],dt;

   real ct;
   real new_rbuf;
   real tol_save;     

   int omt;

   int n,i,j,k,l,Lf,Lc;

   int is,ie,js,je;

   ct=PAMR_get_time(L);

   Lf=PAMR_get_max_lev(PAMR_AMRH);
   Lc=PAMR_get_min_lev(PAMR_AMRH);  //if (PAMR_get_max_lev(PAMR_AMRH)>1) Lc=2; else Lc=1;

   if (AMRD_state!=AMRD_STATE_EVOLVE) return; // if disable, enable(?) reset_AH_shapes below

   if (pre_tstep_global_first)
   {
      for (l=0; l<MAX_BHS; l++) { AH_count[l]=found_AH[l]=found_count_AH[l]=0; freq0[l]=AH_freq[l]; }
      pre_tstep_global_first=0;
   }

   // search for AHs at t>0
   do_repop=do_reinit_ex=got_an_AH=0;

   for (l=0; l<MAX_BHS; l++)
   {
      real prev_AH_R[AH_Nchi[l]*AH_Nphi[l]];
      real prev_AH_xc[2];
      for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) {prev_AH_R[i]=AH_R[l][i];} 
      for (i=0; i<2; i++) {prev_AH_xc[i]=AH_xc[l][i];}
      c_AH=l;
      if (AH_max_iter[l]>0 && L==AH_Lmin[l] && ct>=AH_tmin[l])
      {
         if (AH_count[l]<0) { AH[l]=1; if (AMRD_state==AMRD_STATE_EVOLVE) M=J=0;}
         else if (!(AH_count[l] % freq0[l]) && !(c_AH==3 && ct==0)) // for fourth BH, do not search at t=0
         {
            omt=0; // over-max-tolerance
            AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid0); // AH finder

            if (AH[l]) { freq0[l]=AH_freq_aft[l]; found_AH[l]=1; got_an_AH=1; AH_tol[l]=AH_tol_aft[l]; found_count_AH[l]++; } // if this time found AH

            // if previously found but failed now
            if (found_AH[l] && !AH[l]) 
            {
               if (AH_reset_scale[l]>0)
               {
                  // expand old initial-guess surface
                  if (my_rank==0 && AMRD_evo_trace>=1)
                     printf("t=%lf ... lost AH[%i] ...\n" 
                            "expanding old initial-guess surface by a factor %lf\n",ct,l+1,AH_reset_scale[l]);
                  for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) AH_R[l][i]=prev_AH_R[i]*AH_reset_scale[l];

                  if (!(AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid1)))
                  { 
                     // shrink old initial-guess surface
                     if (my_rank==0 && AMRD_evo_trace>=1) printf("... still can't find one (min_resid=%lf)\n"
                         "... shrinking old initial-guess surface by a factor %lf\n",AH_min_resid1,1/AH_reset_scale[l]);
                     for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) AH_R[l][i]=prev_AH_R[i]/AH_reset_scale[l]; 

                     if (!(AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid2)))
                     {

                        // increase AH_tol to force an AH to be found, starting with shrunken initial-guess surface
                        if (AH_min_resid2<AH_min_resid1 && AH_min_resid2<AH_min_resid0) 
                        {
                           if (my_rank==0 && AMRD_evo_trace>=1)
                              printf("starting from shrunken initial-guess surface\n");
                           for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) AH_R[l][i]=prev_AH_R[i]/AH_reset_scale[l]; 
                           if (my_rank==0 && AMRD_evo_trace>=1) 
                              printf("and temporarily increasing tolerance to %lf; will then use the resulting surface\n"
                              ,AH_min_resid2*AH_omt_scale[l]);
                           omt=1;
                           tol_save=AH_tol[l];
                           AH_tol[l]=AH_min_resid2*AH_omt_scale[l];
                           if (!(AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid2)))
                           {
                              if (my_rank==0 && AMRD_evo_trace>=1) printf("BUG: couldn't find *same* AH\n");
                              if (AH_RESET_AFTER_FAIL) found_AH[l]=0;
                           }
                           if (my_rank==0 && AMRD_evo_trace>=1) printf("setting tolerance back to %lf\n",tol_save);
                           AH_tol[l]=tol_save;
                        }

                        // increase AH_tol to force an AH to be found, starting with expanded initial-guess surface
                        else if (AH_min_resid1<AH_min_resid0) 
                        {  
                           if (my_rank==0 && AMRD_evo_trace>=1)
                              printf("starting from expanded initial-guess surface\n");
                           for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) AH_R[l][i]=prev_AH_R[i]*AH_reset_scale[l]; 
                           if (my_rank==0 && AMRD_evo_trace>=1) 
                              printf("and temporarily increasing tolerance to %lf; use result surface\n",AH_min_resid1*AH_omt_scale[l]);
                           omt=1;
                           tol_save=AH_tol[l];
                           AH_tol[l]=AH_min_resid1*AH_omt_scale[l];
                           if (!(AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid1)))
                           {
                              if (my_rank==0 && AMRD_evo_trace>=1) printf("BUG: couldn't find *same* AH\n");
                              if (AH_RESET_AFTER_FAIL) found_AH[l]=0;
                           }
                           if (my_rank==0 && AMRD_evo_trace>=1) printf("setting tolerance back to %lf\n",tol_save);
                           AH_tol[l]=tol_save;
                        }

                        // increase AH_tol to force an AH to be found, starting with old initial-guess surface
                        else 
                        {
                           if (my_rank==0 && AMRD_evo_trace>=1)
                              printf("starting from old initial-guess surface\n");
                           for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) AH_R[l][i]=prev_AH_R[i]; 
                           if (my_rank==0 && AMRD_evo_trace>=1) 
                              printf("and temporarily increasing tolerance to %lf; use result surface\n",AH_min_resid0*AH_omt_scale[l]);
                           omt=1;
                           tol_save=AH_tol[l];
                           AH_tol[l]=AH_min_resid0*AH_omt_scale[l];
                           if (!(AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid0)))
                           {
                              if (my_rank==0 && AMRD_evo_trace>=1) printf("BUG: couldn't find *same* AH\n");
                              if (AH_RESET_AFTER_FAIL) found_AH[l]=0;
                           }
                           if (my_rank==0 && AMRD_evo_trace>=1) printf("setting tolerance back to %lf\n",tol_save);
                           AH_tol[l]=tol_save;
                        }

                     }                 
                  }
               }
            }

            // if still not found
            if (!(found_AH[l])) for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) {AH_R[l][i]=prev_AH_R[i];}

            // if never found AH
            if (!(found_AH[l])) for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) AH_R[l][i]=AH_r0[l];

            // save AH grid functions if this time found AH
            if (AH[l]) 
            {
               AH_shape[0]=AH_Nchi[l];
               AH_bbox[0]=0;
               if (AH_xc[l][1]<dy) {AH_bbox[1]=M_PI;} else {AH_bbox[1]=2*M_PI;}
               int rank=1;
               sprintf(name,"%sAH_R_%i",AMRD_save_tag,l);
               gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_R[l]);
               sprintf(name,"%sAH_theta_%i",AMRD_save_tag,l);
               gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_theta[l]);
            }

            // fill in excision parameters 
            // ( ex_xc0[0],ex_xc0[1],ex_xc0[2] are filled with coordinate center for excision
            //   and ex_r0[0],ex_r0[1],ex_r0[2] are filled with principle axis radii for excision )
            if (found_AH[l] && AH[l] && AMRD_do_ex)
            {
               fill_ex_params_(AH_R[l],AH_xc[l],ex_r0,ex_xc0,&AH_Nchi[l],&AH_Nphi[l],&dx,&dy,&axisym);

               if (no_AH_intersect(ex_r0,ex_xc0,l))
               {
                 do_reinit_ex=1;
                 do_repop=1;
                 // saves local ex_r0,ex_xc0 to global ex_r, ex_xc
                 ex_r[l][0]=ex_r0[0]; //excision ellipse x-semiaxis
                 ex_r[l][1]=ex_r0[1]; //excision ellipse y-semiaxis
                 ex_xc[l][0]=ex_xc0[0]; //excision ellipse x-coordinate center
                 ex_xc[l][1]=ex_xc0[1]; //excision ellipse y-coordinate center 
               }
            }

         }
         AH_count[l]++;
      }
   }

   // repopulate if needed
   int repop_n=1;
   if (do_repop) AMRD_repopulate(repop_n,ex_repop_io); //TMP FOR REPOP: THIS IS REPOPING OUTER BDY 

   do_reinit_ex=1; //REMOVE THIS LATER (but: not doing PAMR_excision_on at all causes problems)

   // re-initialize mask function
   if (do_reinit_ex)
   {
     remove_redundant_AH();
     PAMR_excision_on("chr",&CFE4D_fill_ex_mask,AMRD_ex,1);
   }

   return;
}

//=============================================================================
// The following routine prints diagnostic quantities
//
// NOTE: at this point, the time sequence is: n,nm1,np1 (unless t=0)
//=============================================================================
void CFE4D_post_tstep(int L)
{
   int itrace=1,valid;
   static int local_first = 1;

   real ct;
   int n,i,j,Lf,Lc;

   int is,ie,js,je;

   ct = PAMR_get_time(L);

   Lf=PAMR_get_max_lev(PAMR_AMRH);
   Lc=PAMR_get_min_lev(PAMR_AMRH);  //if (PAMR_get_max_lev(PAMR_AMRH)>1) Lc=2; else Lc=1;

   if (AMRD_state!=AMRD_STATE_EVOLVE) return; // if disable, enable(?) reset_AH_shapes below

   return;
}


#define ACTION_RELAX 1
#define ACTION_LOP 3
#define ACTION_RESIDUAL 2
//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's)
//=============================================================================
real CFE4D_MG_residual(void)
{
   int action=ACTION_RESIDUAL;
   real norm;

   ldptr_mg();

   if (background || skip_constraints)
   {
      zero_f(zeta_res);
      return 0;
   }

   // solves for zeta conformal factor at t=0; residual 
   mg_sup_(&action,zeta,zeta_rhs,zeta_lop,zeta_res,phi1,
           &AdS_L,mask_mg,phys_bdy,chr_mg,&AMRD_ex, 
           x,y,&norm,&Nx,&Ny);

   return norm;
}


//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual.
//=============================================================================
real CFE4D_MG_relax(void)
{
   int action=ACTION_RELAX;
   real norm;

   ldptr_mg();

   if (background || skip_constraints)
   {
      const_f(zeta,1);
      return 0;
    }

   // solves for zeta conformal factor at t=0; relaxation 
   mg_sup_(&action,zeta,zeta_rhs,zeta_lop,zeta_res,phi1,
           &AdS_L,mask_mg,phys_bdy,chr_mg,&AMRD_ex,
           x,y,&norm,&Nx,&Ny);

   return norm;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f"
//=============================================================================
void CFE4D_L_op(void)
{
   int action=ACTION_LOP;
   real norm;

   ldptr_mg();

   if (background || skip_constraints) 
   {
      zero_f(zeta_lop);
      return;
   }

   // solves for zeta conformal factor at t=0; elliptic operator
   mg_sup_(&action,zeta,zeta_rhs,zeta_lop,zeta_res,phi1,
           &AdS_L,mask_mg,phys_bdy,chr_mg,&AMRD_ex,
           x,y,&norm,&Nx,&Ny);

   return;
}

//=============================================================================
//=============================================================================
void CFE4D_scale_tre(void)
{
}

//=============================================================================
// post-regrid initialization of constant functions
//=============================================================================
void CFE4D_post_regrid(void)
{
   int i;

   if (!background || skip_constraints) return;

   ldptr();

   init_ghb_ads_(gb_tt_n,gb_tx_n,gb_ty_n,gb_xx_n,gb_xy_n,gb_yy_n,psi_n,
                 Hb_t_n,Hb_x_n,Hb_y_n,&AdS_L,x,y,chr,&AMRD_ex,&Nx,&Ny,&regtype);

   for (i=0; i<size; i++)
   {
      gb_tt_nm1[i]=gb_tt_np1[i]=gb_tt_n[i];
      gb_tx_nm1[i]=gb_tx_np1[i]=gb_tx_n[i];
      gb_ty_nm1[i]=gb_ty_np1[i]=gb_ty_n[i];
      gb_xx_nm1[i]=gb_xx_np1[i]=gb_xx_n[i];
      gb_xy_nm1[i]=gb_xy_np1[i]=gb_xy_n[i];
      gb_yy_nm1[i]=gb_yy_np1[i]=gb_yy_n[i];
      psi_nm1[i]=psi_np1[i]=psi_n[i];
   }
}

//=============================================================================
//check-pointing
//=============================================================================
#define CP_DATA_SIZE 50000
void CFE4D_copy_block(char *p, char **q, int n, int dir, int *tot_size)
{
   char *p0;
   int n0;

   if (n==0) return;

   if ((*tot_size+n) > (CP_DATA_SIZE))
      AMRD_stop("CFE4D_copy_block: error ... CP_DATA_SIZE too small\n","");
   *tot_size+=n;

   n0=n;
   p0=p;

   if (dir==AMRD_CP_SAVE) while(n0--) *(*q)++=*p0++;
   else while(n0--) *p0++=*(*q)++;

   return;
}

void CFE4D_cp(int dir, char *data)
{
   int size=0;
   
   if (dir==AMRD_CP_SAVE)
   {
       cp_version=ADS5D_CP_VERSION;
       CFE4D_copy_block((char *)&cp_version,&data,sizeof(int),dir,&size);
   }
}

//=============================================================================
int main(int argc, char **argv)
{
   amrd_set_app_user_cp_hook(CFE4D_cp,CP_DATA_SIZE);
   amrd_set_app_pre_tstep_hook(CFE4D_pre_tstep);
   amrd_set_elliptic_vars_t0_init(CFE4D_elliptic_vars_t0_init);
   amrd(argc,argv,&CFE4D_id,&CFE4D_var_pre_init,
        &CFE4D_var_post_init, &CFE4D_AMRH_var_clear,
        &CFE4D_free_data, &CFE4D_t0_cnst_data,
        &CFE4D_evo_residual, &CFE4D_MG_residual,
        &CFE4D_evolve, &CFE4D_MG_relax, &CFE4D_L_op, 
        &CFE4D_pre_io_calc, &CFE4D_scale_tre, 
        &CFE4D_post_regrid, &CFE4D_post_tstep,
        &CFE4D_fill_ex_mask, &CFE4D_fill_bh_bboxes);
}

