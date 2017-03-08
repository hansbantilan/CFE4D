//============================================================================
// in polar coordinates t,x,y,z = t,x,theta,phi for x in [0,1]
// using r=2*x/(1-x^2) compactification
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

//=============================================================================
// set in fparam for now
//=============================================================================
real AdS_L;

//=============================================================================
// Carsten's constraint-damping parameters
//=============================================================================
real kappa_cd,rho_cd;

//=============================================================================
// scalar field parameters
//=============================================================================
int stype;

//=============================================================================
// id parameters
//=============================================================================

// gaussians
real phi1_amp_1,phi1_r0_1,phi1_delta_1,phi1_x0_1[1],phi1_width_1[1],phi1_amp2_1,phi1_r02_1,phi1_delta2_1,phi1_x02_1[1],phi1_width2_1[1];
int phi1_jj_1,phi1_jj2_1;
real phi1_amp_2,phi1_r0_2,phi1_delta_2,phi1_x0_2[1],phi1_width_2[1],phi1_amp2_2,phi1_r02_2,phi1_delta2_2,phi1_x02_2[1],phi1_width2_2[1];
int phi1_jj_2,phi1_jj2_2;
real phi1_amp_3,phi1_r0_3,phi1_delta_3,phi1_x0_3[1],phi1_width_3[1],phi1_amp2_3,phi1_r02_3,phi1_delta2_3,phi1_x02_3[1],phi1_width2_3[1];
int phi1_jj_3,phi1_jj2_3;

int background,skip_constraints;

// if > 0, initialize with exact BH
real ief_bh_r0;

// new parameters in rtfile
int interptype,i_shift;

//gauge parameters
int gauge_t;
int gauge_i;
real rho1_t,rho2_t,rho3_t,rho4_t,xi1_t,xi2_t,cbulk_t;
real rho1_i,rho2_i,rho3_i,rho4_i,xi1_i,xi2_i,cbulk_i;
real rhoa,rhob;

int cp_version; 

//=============================================================================
// some convenient, "local" global variables
//=============================================================================

real *cl_res;

real *phi1,*phi1_n,*phi1_np1,*phi1_nm1; // MGH, AMRH n/np1/nm1
real *phi1_t,*phi1_t_n;
real *kg_res;

real *gb_tt,*gb_tt_n,*gb_tt_np1,*gb_tt_nm1;
real *gb_tx,*gb_tx_n,*gb_tx_np1,*gb_tx_nm1;
real *gb_xx,*gb_xx_n,*gb_xx_np1,*gb_xx_nm1;
real *psi,*psi_n,*psi_np1,*psi_nm1;
real *gb_tt_t,*gb_tt_t_n;
real *gb_tx_t,*gb_tx_t_n;
real *gb_xx_t,*gb_xx_t_n;
real *psi_t,*psi_t_n;
real *gb_res;

real *db_txtx,*db_txtx_n,*db_txtx_np1,*db_txtx_nm1;
real *db_txty,*db_txty_n,*db_txty_np1,*db_txty_nm1;
real *db_txtz,*db_txtz_n,*db_txtz_np1,*db_txtz_nm1;
real *db_txxy,*db_txxy_n,*db_txxy_np1,*db_txxy_nm1;
real *db_txxz,*db_txxz_n,*db_txxz_np1,*db_txxz_nm1;
real *db_txyz,*db_txyz_n,*db_txyz_np1,*db_txyz_nm1;
real *db_tyty,*db_tyty_n,*db_tyty_np1,*db_tyty_nm1;
real *db_tytz,*db_tytz_n,*db_tytz_np1,*db_tytz_nm1;
real *db_tyxy,*db_tyxy_n,*db_tyxy_np1,*db_tyxy_nm1;
real *db_tyxz,*db_tyxz_n,*db_tyxz_np1,*db_tyxz_nm1;
real *db_tyyz,*db_tyyz_n,*db_tyyz_np1,*db_tyyz_nm1;
real *db_tztz,*db_tztz_n,*db_tztz_np1,*db_tztz_nm1;
real *db_tzxy,*db_tzxy_n,*db_tzxy_np1,*db_tzxy_nm1;
real *db_tzxz,*db_tzxz_n,*db_tzxz_np1,*db_tzxz_nm1;
real *db_tzyz,*db_tzyz_n,*db_tzyz_np1,*db_tzyz_nm1;
real *db_xyxy,*db_xyxy_n,*db_xyxy_np1,*db_xyxy_nm1;
real *db_xyxz,*db_xyxz_n,*db_xyxz_np1,*db_xyxz_nm1;
real *db_xyyz,*db_xyyz_n,*db_xyyz_np1,*db_xyyz_nm1;
real *db_xzxz,*db_xzxz_n,*db_xzxz_np1,*db_xzxz_nm1;
real *db_xzyz,*db_xzyz_n,*db_xzyz_np1,*db_xzyz_nm1;
real *db_yzyz,*db_yzyz_n,*db_yzyz_np1,*db_yzyz_nm1;
real *db_txtx_t,*db_txtx_t_n;
real *db_txty_t,*db_txty_t_n;
real *db_txtz_t,*db_txtz_t_n;
real *db_txxy_t,*db_txxy_t_n;
real *db_txxz_t,*db_txxz_t_n;
real *db_txyz_t,*db_txyz_t_n;
real *db_tyty_t,*db_tyty_t_n;
real *db_tytz_t,*db_tytz_t_n;
real *db_tyxy_t,*db_tyxy_t_n;
real *db_tyxz_t,*db_tyxz_t_n;
real *db_tyyz_t,*db_tyyz_t_n;
real *db_tztz_t,*db_tztz_t_n;
real *db_tzxy_t,*db_tzxy_t_n;
real *db_tzxz_t,*db_tzxz_t_n;
real *db_tzyz_t,*db_tzyz_t_n;
real *db_xyxy_t,*db_xyxy_t_n;
real *db_xyxz_t,*db_xyxz_t_n;
real *db_xyyz_t,*db_xyyz_t_n;
real *db_xzxz_t,*db_xzxz_t_n;
real *db_xzyz_t,*db_xzyz_t_n;
real *db_yzyz_t,*db_yzyz_t_n;
real *db_res;

real *Hb_t,*Hb_t_n,*Hb_t_np1,*Hb_t_nm1;
real *Hb_x,*Hb_x_n,*Hb_x_np1,*Hb_x_nm1;
real *Hb_t_t,*Hb_t_t_n;
real *Hb_x_t,*Hb_x_t_n;
real *hb_t_res,*hb_i_res;

real *Hb_t_0,*Hb_x_0;

real *zetab,*zetab_res,*zetab_lop,*zetab_rhs;

real *w1,*mg_w1;
real *w2,*mg_w2;
real *w3,*mg_w3;
real *w4,*mg_w4;

real *mask,*mask_mg,*chr,*chr_mg;

real *efe_all_ires,*iresall;
real *phj_ires,*phj;
real *pij_ires,*pij;
real *alphaq_ires,*alphaq;
real *thetap_ires,*thetap;

real qs_tt,qs_psi,qs_mass;
	
real *g_norms;

real *x;
int shape[1],ghost_width[2],Nx,phys_bdy[2],size,g_rank;
real base_bbox[2],bbox[2],dx,dt,dx_Lc;
int g_L;

int cl_res_gfn;

int phi1_gfn,phi1_n_gfn,phi1_np1_gfn,phi1_nm1_gfn; 
int phi1_t_gfn,phi1_t_n_gfn;
int kg_res_gfn;

int gb_tt_gfn,gb_tt_n_gfn,gb_tt_np1_gfn,gb_tt_nm1_gfn;
int gb_tx_gfn,gb_tx_n_gfn,gb_tx_np1_gfn,gb_tx_nm1_gfn;
int gb_xx_gfn,gb_xx_n_gfn,gb_xx_np1_gfn,gb_xx_nm1_gfn;
int psi_gfn,psi_n_gfn,psi_np1_gfn,psi_nm1_gfn;
int gb_tt_t_gfn,gb_tt_t_n_gfn;
int gb_tx_t_gfn,gb_tx_t_n_gfn;
int gb_xx_t_gfn,gb_xx_t_n_gfn;
int psi_t_gfn,psi_t_n_gfn;
int gb_res_gfn;

int db_txtx_gfn,db_txtx_n_gfn,db_txtx_np1_gfn,db_txtx_nm1_gfn;
int db_txty_gfn,db_txty_n_gfn,db_txty_np1_gfn,db_txty_nm1_gfn;
int db_txtz_gfn,db_txtz_n_gfn,db_txtz_np1_gfn,db_txtz_nm1_gfn;
int db_txxy_gfn,db_txxy_n_gfn,db_txxy_np1_gfn,db_txxy_nm1_gfn;
int db_txxz_gfn,db_txxz_n_gfn,db_txxz_np1_gfn,db_txxz_nm1_gfn;
int db_txyz_gfn,db_txyz_n_gfn,db_txyz_np1_gfn,db_txyz_nm1_gfn;
int db_tyty_gfn,db_tyty_n_gfn,db_tyty_np1_gfn,db_tyty_nm1_gfn;
int db_tytz_gfn,db_tytz_n_gfn,db_tytz_np1_gfn,db_tytz_nm1_gfn;
int db_tyxy_gfn,db_tyxy_n_gfn,db_tyxy_np1_gfn,db_tyxy_nm1_gfn;
int db_tyxz_gfn,db_tyxz_n_gfn,db_tyxz_np1_gfn,db_tyxz_nm1_gfn;
int db_tyyz_gfn,db_tyyz_n_gfn,db_tyyz_np1_gfn,db_tyyz_nm1_gfn;
int db_tztz_gfn,db_tztz_n_gfn,db_tztz_np1_gfn,db_tztz_nm1_gfn;
int db_tzxy_gfn,db_tzxy_n_gfn,db_tzxy_np1_gfn,db_tzxy_nm1_gfn;
int db_tzxz_gfn,db_tzxz_n_gfn,db_tzxz_np1_gfn,db_tzxz_nm1_gfn;
int db_tzyz_gfn,db_tzyz_n_gfn,db_tzyz_np1_gfn,db_tzyz_nm1_gfn;
int db_xyxy_gfn,db_xyxy_n_gfn,db_xyxy_np1_gfn,db_xyxy_nm1_gfn;
int db_xyxz_gfn,db_xyxz_n_gfn,db_xyxz_np1_gfn,db_xyxz_nm1_gfn;
int db_xyyz_gfn,db_xyyz_n_gfn,db_xyyz_np1_gfn,db_xyyz_nm1_gfn;
int db_xzxz_gfn,db_xzxz_n_gfn,db_xzxz_np1_gfn,db_xzxz_nm1_gfn;
int db_xzyz_gfn,db_xzyz_n_gfn,db_xzyz_np1_gfn,db_xzyz_nm1_gfn;
int db_yzyz_gfn,db_yzyz_n_gfn,db_yzyz_np1_gfn,db_yzyz_nm1_gfn;
int db_txtx_t_gfn,db_txtx_t_n_gfn;
int db_txty_t_gfn,db_txty_t_n_gfn;
int db_txtz_t_gfn,db_txtz_t_n_gfn;
int db_txxy_t_gfn,db_txxy_t_n_gfn;
int db_txxz_t_gfn,db_txxz_t_n_gfn;
int db_txyz_t_gfn,db_txyz_t_n_gfn;
int db_tyty_t_gfn,db_tyty_t_n_gfn;
int db_tytz_t_gfn,db_tytz_t_n_gfn;
int db_tyxy_t_gfn,db_tyxy_t_n_gfn;
int db_tyxz_t_gfn,db_tyxz_t_n_gfn;
int db_tyyz_t_gfn,db_tyyz_t_n_gfn;
int db_tztz_t_gfn,db_tztz_t_n_gfn;
int db_tzxy_t_gfn,db_tzxy_t_n_gfn;
int db_tzxz_t_gfn,db_tzxz_t_n_gfn;
int db_tzyz_t_gfn,db_tzyz_t_n_gfn;
int db_xyxy_t_gfn,db_xyxy_t_n_gfn;
int db_xyxz_t_gfn,db_xyxz_t_n_gfn;
int db_xyyz_t_gfn,db_xyyz_t_n_gfn;
int db_xzxz_t_gfn,db_xzxz_t_n_gfn;
int db_xzyz_t_gfn,db_xzyz_t_n_gfn;
int db_yzyz_t_gfn,db_yzyz_t_n_gfn;
int db_res_gfn;

int Hb_t_gfn,Hb_t_n_gfn,Hb_t_np1_gfn,Hb_t_nm1_gfn;
int Hb_x_gfn,Hb_x_n_gfn,Hb_x_np1_gfn,Hb_x_nm1_gfn;
int Hb_t_t_gfn,Hb_t_t_n_gfn;
int Hb_x_t_gfn,Hb_x_t_n_gfn;
int hb_t_res_gfn,hb_i_res_gfn;

int Hb_t_0_gfn,Hb_x_0_gfn;	

int zetab_gfn,zetab_res_gfn,zetab_lop_gfn,zetab_rhs_gfn;

int w1_gfn,mg_w1_gfn;
int w2_gfn,mg_w2_gfn;
int w3_gfn,mg_w3_gfn;
int w4_gfn,mg_w4_gfn;

int mask_gfn,mask_mg_gfn,chr_gfn,chr_mg_gfn;

int efe_all_ires_gfn,iresall_gfn;
int phj_ires_gfn,phj_gfn;
int pij_ires_gfn,pij_gfn;
int alphaq_ires_gfn,alphaq_gfn;
int thetap_ires_gfn,thetap_gfn;

//=============================================================================
// arrays holding quasiset of CFT on ESU at outer bdy
//=============================================================================
real *lqstt0,*lqspsi0,*lqsmass0;
real *qstt0,*qspsi0,*qsmass0;

//=============================================================================
// excision parameters
//=============================================================================
real lxh;
real xh;
real ex_rbuf;

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
    if ((phi1_t_gfn   = PAMR_get_gfn("phi1_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((phi1_t_n_gfn = PAMR_get_gfn("phi1_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((kg_res_gfn   = PAMR_get_gfn("kg_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_tt_gfn     = PAMR_get_gfn("gb_tt",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_nm1_gfn = PAMR_get_gfn("gb_tt",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_n_gfn   = PAMR_get_gfn("gb_tt",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_np1_gfn = PAMR_get_gfn("gb_tt",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_gfn     = PAMR_get_gfn("gb_tx",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_nm1_gfn = PAMR_get_gfn("gb_tx",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_n_gfn   = PAMR_get_gfn("gb_tx",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_np1_gfn = PAMR_get_gfn("gb_tx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_gfn     = PAMR_get_gfn("gb_xx",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_nm1_gfn = PAMR_get_gfn("gb_xx",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_n_gfn   = PAMR_get_gfn("gb_xx",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_np1_gfn = PAMR_get_gfn("gb_xx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_gfn     = PAMR_get_gfn("psi",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_nm1_gfn = PAMR_get_gfn("psi",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_n_gfn   = PAMR_get_gfn("psi",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_np1_gfn = PAMR_get_gfn("psi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_t_gfn = PAMR_get_gfn("gb_tt_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_t_n_gfn = PAMR_get_gfn("gb_tt_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_t_gfn = PAMR_get_gfn("gb_tx_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_t_n_gfn = PAMR_get_gfn("gb_tx_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_t_gfn = PAMR_get_gfn("gb_xx_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_t_n_gfn = PAMR_get_gfn("gb_xx_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_t_gfn = PAMR_get_gfn("psi_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_t_n_gfn = PAMR_get_gfn("psi_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_res_gfn    = PAMR_get_gfn("gb_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((db_txtx_gfn     = PAMR_get_gfn("db_txtx",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txtx_nm1_gfn = PAMR_get_gfn("db_txtx",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txtx_n_gfn   = PAMR_get_gfn("db_txtx",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txtx_np1_gfn = PAMR_get_gfn("db_txtx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txty_gfn     = PAMR_get_gfn("db_txty",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txty_nm1_gfn = PAMR_get_gfn("db_txty",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txty_n_gfn   = PAMR_get_gfn("db_txty",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txty_np1_gfn = PAMR_get_gfn("db_txty",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txtz_gfn     = PAMR_get_gfn("db_txtz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txtz_nm1_gfn = PAMR_get_gfn("db_txtz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txtz_n_gfn   = PAMR_get_gfn("db_txtz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txtz_np1_gfn = PAMR_get_gfn("db_txtz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txxy_gfn     = PAMR_get_gfn("db_txxy",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txxy_nm1_gfn = PAMR_get_gfn("db_txxy",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txxy_n_gfn   = PAMR_get_gfn("db_txxy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txxy_np1_gfn = PAMR_get_gfn("db_txxy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txxz_gfn     = PAMR_get_gfn("db_txxz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txxz_nm1_gfn = PAMR_get_gfn("db_txxz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txxz_n_gfn   = PAMR_get_gfn("db_txxz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txxz_np1_gfn = PAMR_get_gfn("db_txxz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txyz_gfn     = PAMR_get_gfn("db_txyz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txyz_nm1_gfn = PAMR_get_gfn("db_txyz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txyz_n_gfn   = PAMR_get_gfn("db_txyz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txyz_np1_gfn = PAMR_get_gfn("db_txyz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyty_gfn     = PAMR_get_gfn("db_tyty",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyty_nm1_gfn = PAMR_get_gfn("db_tyty",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyty_n_gfn   = PAMR_get_gfn("db_tyty",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyty_np1_gfn = PAMR_get_gfn("db_tyty",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tytz_gfn     = PAMR_get_gfn("db_tytz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tytz_nm1_gfn = PAMR_get_gfn("db_tytz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tytz_n_gfn   = PAMR_get_gfn("db_tytz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tytz_np1_gfn = PAMR_get_gfn("db_tytz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyxy_gfn     = PAMR_get_gfn("db_tyxy",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyxy_nm1_gfn = PAMR_get_gfn("db_tyxy",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyxy_n_gfn   = PAMR_get_gfn("db_tyxy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyxy_np1_gfn = PAMR_get_gfn("db_tyxy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyxz_gfn     = PAMR_get_gfn("db_tyxz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyxz_nm1_gfn = PAMR_get_gfn("db_tyxz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyxz_n_gfn   = PAMR_get_gfn("db_tyxz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyxz_np1_gfn = PAMR_get_gfn("db_tyxz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyyz_gfn     = PAMR_get_gfn("db_tyyz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyyz_nm1_gfn = PAMR_get_gfn("db_tyyz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyyz_n_gfn   = PAMR_get_gfn("db_tyyz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyyz_np1_gfn = PAMR_get_gfn("db_tyyz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tztz_gfn     = PAMR_get_gfn("db_tztz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tztz_nm1_gfn = PAMR_get_gfn("db_tztz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tztz_n_gfn   = PAMR_get_gfn("db_tztz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tztz_np1_gfn = PAMR_get_gfn("db_tztz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzxy_gfn     = PAMR_get_gfn("db_tzxy",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzxy_nm1_gfn = PAMR_get_gfn("db_tzxy",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzxy_n_gfn   = PAMR_get_gfn("db_tzxy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzxy_np1_gfn = PAMR_get_gfn("db_tzxy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzxz_gfn     = PAMR_get_gfn("db_tzxz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzxz_nm1_gfn = PAMR_get_gfn("db_tzxz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzxz_n_gfn   = PAMR_get_gfn("db_tzxz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzxz_np1_gfn = PAMR_get_gfn("db_tzxz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzyz_gfn     = PAMR_get_gfn("db_tzyz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzyz_nm1_gfn = PAMR_get_gfn("db_tzyz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzyz_n_gfn   = PAMR_get_gfn("db_tzyz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzyz_np1_gfn = PAMR_get_gfn("db_tzyz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyxy_gfn     = PAMR_get_gfn("db_xyxy",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyxy_nm1_gfn = PAMR_get_gfn("db_xyxy",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyxy_n_gfn   = PAMR_get_gfn("db_xyxy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyxy_np1_gfn = PAMR_get_gfn("db_xyxy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyxz_gfn     = PAMR_get_gfn("db_xyxz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyxz_nm1_gfn = PAMR_get_gfn("db_xyxz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyxz_n_gfn   = PAMR_get_gfn("db_xyxz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyxz_np1_gfn = PAMR_get_gfn("db_xyxz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyyz_gfn     = PAMR_get_gfn("db_xyyz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyyz_nm1_gfn = PAMR_get_gfn("db_xyyz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyyz_n_gfn   = PAMR_get_gfn("db_xyyz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyyz_np1_gfn = PAMR_get_gfn("db_xyyz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xzxz_gfn     = PAMR_get_gfn("db_xzxz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xzxz_nm1_gfn = PAMR_get_gfn("db_xzxz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xzxz_n_gfn   = PAMR_get_gfn("db_xzxz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xzxz_np1_gfn = PAMR_get_gfn("db_xzxz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xzyz_gfn     = PAMR_get_gfn("db_xzyz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xzyz_nm1_gfn = PAMR_get_gfn("db_xzyz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xzyz_n_gfn   = PAMR_get_gfn("db_xzyz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xzyz_np1_gfn = PAMR_get_gfn("db_xzyz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_yzyz_gfn     = PAMR_get_gfn("db_yzyz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_yzyz_nm1_gfn = PAMR_get_gfn("db_yzyz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((db_yzyz_n_gfn   = PAMR_get_gfn("db_yzyz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_yzyz_np1_gfn = PAMR_get_gfn("db_yzyz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txtx_t_gfn     = PAMR_get_gfn("db_txtx_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txtx_t_n_gfn   = PAMR_get_gfn("db_txtx_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txty_t_gfn     = PAMR_get_gfn("db_txty_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txty_t_n_gfn   = PAMR_get_gfn("db_txty_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txtz_t_gfn     = PAMR_get_gfn("db_txtz_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txtz_t_n_gfn   = PAMR_get_gfn("db_txtz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txxy_t_gfn     = PAMR_get_gfn("db_txxy_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txxy_t_n_gfn   = PAMR_get_gfn("db_txxy_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txxz_t_gfn     = PAMR_get_gfn("db_txxz_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txxz_t_n_gfn   = PAMR_get_gfn("db_txxz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txyz_t_gfn     = PAMR_get_gfn("db_txyz_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_txyz_t_n_gfn   = PAMR_get_gfn("db_txyz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyty_t_gfn     = PAMR_get_gfn("db_tyty_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyty_t_n_gfn   = PAMR_get_gfn("db_tyty_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tytz_t_gfn     = PAMR_get_gfn("db_tytz_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tytz_t_n_gfn   = PAMR_get_gfn("db_tytz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyxy_t_gfn     = PAMR_get_gfn("db_tyxy_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyxy_t_n_gfn   = PAMR_get_gfn("db_tyxy_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyxz_t_gfn     = PAMR_get_gfn("db_tyxz_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyxz_t_n_gfn   = PAMR_get_gfn("db_tyxz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyyz_t_gfn     = PAMR_get_gfn("db_tyyz_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tyyz_t_n_gfn   = PAMR_get_gfn("db_tyyz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tztz_t_gfn     = PAMR_get_gfn("db_tztz_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tztz_t_n_gfn   = PAMR_get_gfn("db_tztz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzxy_t_gfn     = PAMR_get_gfn("db_tzxy_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzxy_t_n_gfn   = PAMR_get_gfn("db_tzxy_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzxz_t_gfn     = PAMR_get_gfn("db_tzxz_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzxz_t_n_gfn   = PAMR_get_gfn("db_tzxz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzyz_t_gfn     = PAMR_get_gfn("db_tzyz_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_tzyz_t_n_gfn   = PAMR_get_gfn("db_tzyz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyxy_t_gfn     = PAMR_get_gfn("db_xyxy_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyxy_t_n_gfn   = PAMR_get_gfn("db_xyxy_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyxz_t_gfn     = PAMR_get_gfn("db_xyxz_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyxz_t_n_gfn   = PAMR_get_gfn("db_xyxz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyyz_t_gfn     = PAMR_get_gfn("db_xyyz_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xyyz_t_n_gfn   = PAMR_get_gfn("db_xyyz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xzxz_t_gfn     = PAMR_get_gfn("db_xzxz_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xzxz_t_n_gfn   = PAMR_get_gfn("db_xzxz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xzyz_t_gfn     = PAMR_get_gfn("db_xzyz_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_xzyz_t_n_gfn   = PAMR_get_gfn("db_xzyz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_yzyz_t_gfn     = PAMR_get_gfn("db_yzyz_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((db_yzyz_t_n_gfn   = PAMR_get_gfn("db_yzyz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((db_res_gfn    = PAMR_get_gfn("db_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((Hb_t_gfn      = PAMR_get_gfn("Hb_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_nm1_gfn  = PAMR_get_gfn("Hb_t",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_n_gfn    = PAMR_get_gfn("Hb_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_np1_gfn  = PAMR_get_gfn("Hb_t",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_gfn      = PAMR_get_gfn("Hb_x",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_nm1_gfn  = PAMR_get_gfn("Hb_x",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_n_gfn    = PAMR_get_gfn("Hb_x",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_np1_gfn  = PAMR_get_gfn("Hb_x",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_t_gfn  = PAMR_get_gfn("Hb_t_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_t_n_gfn  = PAMR_get_gfn("Hb_t_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_t_gfn  = PAMR_get_gfn("Hb_x_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_t_n_gfn  = PAMR_get_gfn("Hb_x_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((hb_t_res_gfn  = PAMR_get_gfn("hb_t_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((hb_i_res_gfn  = PAMR_get_gfn("hb_i_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((Hb_t_0_gfn  = PAMR_get_gfn("Hb_t_0",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_0_gfn  = PAMR_get_gfn("Hb_x_0",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((zetab_gfn     = PAMR_get_gfn("zetab",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((zetab_res_gfn = PAMR_get_gfn("zetab_res",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((zetab_lop_gfn = PAMR_get_gfn("zetab_lop",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((zetab_rhs_gfn = PAMR_get_gfn("zetab_rhs",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    if ((w1_gfn   = PAMR_get_gfn("w1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((w2_gfn   = PAMR_get_gfn("w2",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((w3_gfn   = PAMR_get_gfn("w3",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((w4_gfn   = PAMR_get_gfn("w4",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((mg_w1_gfn   = PAMR_get_gfn("mg_w1",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mg_w2_gfn   = PAMR_get_gfn("mg_w2",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mg_w3_gfn   = PAMR_get_gfn("mg_w3",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mg_w4_gfn   = PAMR_get_gfn("mg_w4",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    if ((efe_all_ires_gfn= PAMR_get_gfn("efe_all_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((iresall_gfn  = PAMR_get_gfn("iresall",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((phj_ires_gfn= PAMR_get_gfn("phj_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((phj_gfn  = PAMR_get_gfn("phj",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((pij_ires_gfn  = PAMR_get_gfn("pij_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((pij_gfn  = PAMR_get_gfn("pij",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((alphaq_ires_gfn  = PAMR_get_gfn("alphaq_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((alphaq_gfn  = PAMR_get_gfn("alphaq",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((thetap_ires_gfn  = PAMR_get_gfn("thetap_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((thetap_gfn  = PAMR_get_gfn("thetap",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((mask_mg_gfn = PAMR_get_gfn("cmask",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_gfn    = PAMR_get_gfn("cmask",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((chr_gfn     = PAMR_get_gfn("chr",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((chr_mg_gfn  = PAMR_get_gfn("chr",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    g_norms=AMRD_get_global_norms();
}

//=============================================================================
// call with valid iter to set up globals:
//=============================================================================
void ldptr_bbox(void)
{
   real dx0[1];
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

   if ((bbox[0]-base_bbox[0])<dx/2) phys_bdy[0]=1; else phys_bdy[0]=0;
   if ((base_bbox[1]-bbox[1])<dx/2) phys_bdy[1]=1; else phys_bdy[1]=0;

   Nx=shape[0];

   size=Nx;
}

void ldptr(void)
{
   real *x0[1],*gfs[PAMR_MAX_GFNS];

   ldptr_bbox();

   PAMR_get_g_x(x0);

   x=x0[0];

   PAMR_get_g_gfs(gfs);

   cl_res   = gfs[cl_res_gfn-1];

   phi1     = gfs[phi1_gfn-1];
   phi1_n   = gfs[phi1_n_gfn-1];
   phi1_np1 = gfs[phi1_np1_gfn-1];
   phi1_nm1 = gfs[phi1_nm1_gfn-1];
   phi1_t   = gfs[phi1_t_gfn-1];
   phi1_t_n = gfs[phi1_t_n_gfn-1];
   kg_res   = gfs[kg_res_gfn-1];

   gb_tt     = gfs[gb_tt_gfn-1];
   gb_tt_n   = gfs[gb_tt_n_gfn-1];
   gb_tt_np1 = gfs[gb_tt_np1_gfn-1];
   gb_tt_nm1 = gfs[gb_tt_nm1_gfn-1];
   gb_tx     = gfs[gb_tx_gfn-1];
   gb_tx_n   = gfs[gb_tx_n_gfn-1];
   gb_tx_np1 = gfs[gb_tx_np1_gfn-1];
   gb_tx_nm1 = gfs[gb_tx_nm1_gfn-1];
   gb_xx     = gfs[gb_xx_gfn-1];
   gb_xx_n   = gfs[gb_xx_n_gfn-1];
   gb_xx_np1 = gfs[gb_xx_np1_gfn-1];
   gb_xx_nm1 = gfs[gb_xx_nm1_gfn-1];
   psi     = gfs[psi_gfn-1];
   psi_n   = gfs[psi_n_gfn-1];
   psi_np1 = gfs[psi_np1_gfn-1];
   psi_nm1 = gfs[psi_nm1_gfn-1];
   gb_tt_t = gfs[gb_tt_t_gfn-1];
   gb_tt_t_n = gfs[gb_tt_t_n_gfn-1];
   gb_tx_t = gfs[gb_tx_t_gfn-1];
   gb_tx_t_n = gfs[gb_tx_t_n_gfn-1];
   gb_xx_t = gfs[gb_xx_t_gfn-1];
   gb_xx_t_n = gfs[gb_xx_t_n_gfn-1];
   psi_t = gfs[psi_t_gfn-1];
   psi_t_n = gfs[psi_t_n_gfn-1];
   gb_res    = gfs[gb_res_gfn-1];

   db_txtx     = gfs[db_txtx_gfn-1];
   db_txtx_n   = gfs[db_txtx_n_gfn-1];
   db_txtx_np1 = gfs[db_txtx_np1_gfn-1];
   db_txtx_nm1 = gfs[db_txtx_nm1_gfn-1];
   db_txty     = gfs[db_txty_gfn-1];
   db_txty_n   = gfs[db_txty_n_gfn-1];
   db_txty_np1 = gfs[db_txty_np1_gfn-1];
   db_txty_nm1 = gfs[db_txty_nm1_gfn-1];
   db_txtz     = gfs[db_txtz_gfn-1];
   db_txtz_n   = gfs[db_txtz_n_gfn-1];
   db_txtz_np1 = gfs[db_txtz_np1_gfn-1];
   db_txtz_nm1 = gfs[db_txtz_nm1_gfn-1];
   db_txxy     = gfs[db_txxy_gfn-1];
   db_txxy_n   = gfs[db_txxy_n_gfn-1];
   db_txxy_np1 = gfs[db_txxy_np1_gfn-1];
   db_txxy_nm1 = gfs[db_txxy_nm1_gfn-1];
   db_txxz     = gfs[db_txxz_gfn-1];
   db_txxz_n   = gfs[db_txxz_n_gfn-1];
   db_txxz_np1 = gfs[db_txxz_np1_gfn-1];
   db_txxz_nm1 = gfs[db_txxz_nm1_gfn-1];
   db_txyz     = gfs[db_txyz_gfn-1];
   db_txyz_n   = gfs[db_txyz_n_gfn-1];
   db_txyz_np1 = gfs[db_txyz_np1_gfn-1];
   db_txyz_nm1 = gfs[db_txyz_nm1_gfn-1];
   db_tyty     = gfs[db_tyty_gfn-1];
   db_tyty_n   = gfs[db_tyty_n_gfn-1];
   db_tyty_np1 = gfs[db_tyty_np1_gfn-1];
   db_tyty_nm1 = gfs[db_tyty_nm1_gfn-1];
   db_tytz     = gfs[db_tytz_gfn-1];
   db_tytz_n   = gfs[db_tytz_n_gfn-1];
   db_tytz_np1 = gfs[db_tytz_np1_gfn-1];
   db_tytz_nm1 = gfs[db_tytz_nm1_gfn-1];
   db_tyxy     = gfs[db_tyxy_gfn-1];
   db_tyxy_n   = gfs[db_tyxy_n_gfn-1];
   db_tyxy_np1 = gfs[db_tyxy_np1_gfn-1];
   db_tyxy_nm1 = gfs[db_tyxy_nm1_gfn-1];
   db_tyxz     = gfs[db_tyxz_gfn-1];
   db_tyxz_n   = gfs[db_tyxz_n_gfn-1];
   db_tyxz_np1 = gfs[db_tyxz_np1_gfn-1];
   db_tyxz_nm1 = gfs[db_tyxz_nm1_gfn-1];
   db_tyyz     = gfs[db_tyyz_gfn-1];
   db_tyyz_n   = gfs[db_tyyz_n_gfn-1];
   db_tyyz_np1 = gfs[db_tyyz_np1_gfn-1];
   db_tyyz_nm1 = gfs[db_tyyz_nm1_gfn-1];
   db_tztz     = gfs[db_tztz_gfn-1];
   db_tztz_n   = gfs[db_tztz_n_gfn-1];
   db_tztz_np1 = gfs[db_tztz_np1_gfn-1];
   db_tztz_nm1 = gfs[db_tztz_nm1_gfn-1];
   db_tzxy     = gfs[db_tzxy_gfn-1];
   db_tzxy_n   = gfs[db_tzxy_n_gfn-1];
   db_tzxy_np1 = gfs[db_tzxy_np1_gfn-1];
   db_tzxy_nm1 = gfs[db_tzxy_nm1_gfn-1];
   db_tzxz     = gfs[db_tzxz_gfn-1];
   db_tzxz_n   = gfs[db_tzxz_n_gfn-1];
   db_tzxz_np1 = gfs[db_tzxz_np1_gfn-1];
   db_tzxz_nm1 = gfs[db_tzxz_nm1_gfn-1];
   db_tzyz     = gfs[db_tzyz_gfn-1];
   db_tzyz_n   = gfs[db_tzyz_n_gfn-1];
   db_tzyz_np1 = gfs[db_tzyz_np1_gfn-1];
   db_tzyz_nm1 = gfs[db_tzyz_nm1_gfn-1];
   db_xyxy     = gfs[db_xyxy_gfn-1];
   db_xyxy_n   = gfs[db_xyxy_n_gfn-1];
   db_xyxy_np1 = gfs[db_xyxy_np1_gfn-1];
   db_xyxy_nm1 = gfs[db_xyxy_nm1_gfn-1];
   db_xyxz     = gfs[db_xyxz_gfn-1];
   db_xyxz_n   = gfs[db_xyxz_n_gfn-1];
   db_xyxz_np1 = gfs[db_xyxz_np1_gfn-1];
   db_xyxz_nm1 = gfs[db_xyxz_nm1_gfn-1];
   db_xyyz     = gfs[db_xyyz_gfn-1];
   db_xyyz_n   = gfs[db_xyyz_n_gfn-1];
   db_xyyz_np1 = gfs[db_xyyz_np1_gfn-1];
   db_xyyz_nm1 = gfs[db_xyyz_nm1_gfn-1];
   db_xzxz     = gfs[db_xzxz_gfn-1];
   db_xzxz_n   = gfs[db_xzxz_n_gfn-1];
   db_xzxz_np1 = gfs[db_xzxz_np1_gfn-1];
   db_xzxz_nm1 = gfs[db_xzxz_nm1_gfn-1];
   db_xzyz     = gfs[db_xzyz_gfn-1];
   db_xzyz_n   = gfs[db_xzyz_n_gfn-1];
   db_xzyz_np1 = gfs[db_xzyz_np1_gfn-1];
   db_xzyz_nm1 = gfs[db_xzyz_nm1_gfn-1];
   db_yzyz     = gfs[db_yzyz_gfn-1];
   db_yzyz_n   = gfs[db_yzyz_n_gfn-1];
   db_yzyz_np1 = gfs[db_yzyz_np1_gfn-1];
   db_yzyz_nm1 = gfs[db_yzyz_nm1_gfn-1];
   db_txtx_t     = gfs[db_txtx_t_gfn-1];
   db_txtx_t_n   = gfs[db_txtx_t_n_gfn-1];
   db_txty_t     = gfs[db_txty_t_gfn-1];
   db_txty_t_n   = gfs[db_txty_t_n_gfn-1];
   db_txtz_t     = gfs[db_txtz_t_gfn-1];
   db_txtz_t_n   = gfs[db_txtz_t_n_gfn-1];
   db_txxy_t     = gfs[db_txxy_t_gfn-1];
   db_txxy_t_n   = gfs[db_txxy_t_n_gfn-1];
   db_txxz_t     = gfs[db_txxz_t_gfn-1];
   db_txxz_t_n   = gfs[db_txxz_t_n_gfn-1];
   db_txyz_t     = gfs[db_txyz_t_gfn-1];
   db_txyz_t_n   = gfs[db_txyz_t_n_gfn-1];
   db_tyty_t     = gfs[db_tyty_t_gfn-1];
   db_tyty_t_n   = gfs[db_tyty_t_n_gfn-1];
   db_tytz_t     = gfs[db_tytz_t_gfn-1];
   db_tytz_t_n   = gfs[db_tytz_t_n_gfn-1];
   db_tyxy_t     = gfs[db_tyxy_t_gfn-1];
   db_tyxy_t_n   = gfs[db_tyxy_t_n_gfn-1];
   db_tyxz_t     = gfs[db_tyxz_t_gfn-1];
   db_tyxz_t_n   = gfs[db_tyxz_t_n_gfn-1];
   db_tyyz_t     = gfs[db_tyyz_t_gfn-1];
   db_tyyz_t_n   = gfs[db_tyyz_t_n_gfn-1];
   db_tztz_t     = gfs[db_tztz_t_gfn-1];
   db_tztz_t_n   = gfs[db_tztz_t_n_gfn-1];
   db_tzxy_t     = gfs[db_tzxy_t_gfn-1];
   db_tzxy_t_n   = gfs[db_tzxy_t_n_gfn-1];
   db_tzxz_t     = gfs[db_tzxz_t_gfn-1];
   db_tzxz_t_n   = gfs[db_tzxz_t_n_gfn-1];
   db_tzyz_t     = gfs[db_tzyz_t_gfn-1];
   db_tzyz_t_n   = gfs[db_tzyz_t_n_gfn-1];
   db_xyxy_t     = gfs[db_xyxy_t_gfn-1];
   db_xyxy_t_n   = gfs[db_xyxy_t_n_gfn-1];
   db_xyxz_t     = gfs[db_xyxz_t_gfn-1];
   db_xyxz_t_n   = gfs[db_xyxz_t_n_gfn-1];
   db_xyyz_t     = gfs[db_xyyz_t_gfn-1];
   db_xyyz_t_n   = gfs[db_xyyz_t_n_gfn-1];
   db_xzxz_t     = gfs[db_xzxz_t_gfn-1];
   db_xzxz_t_n   = gfs[db_xzxz_t_n_gfn-1];
   db_xzyz_t     = gfs[db_xzyz_t_gfn-1];
   db_xzyz_t_n   = gfs[db_xzyz_t_n_gfn-1];
   db_yzyz_t     = gfs[db_yzyz_t_gfn-1];
   db_yzyz_t_n   = gfs[db_yzyz_t_n_gfn-1];
   db_res    = gfs[db_res_gfn-1];

   Hb_t      = gfs[Hb_t_gfn-1];
   Hb_t_n    = gfs[Hb_t_n_gfn-1];
   Hb_t_nm1  = gfs[Hb_t_nm1_gfn-1];
   Hb_t_np1  = gfs[Hb_t_np1_gfn-1];
   Hb_x      = gfs[Hb_x_gfn-1];
   Hb_x_n    = gfs[Hb_x_n_gfn-1];
   Hb_x_nm1  = gfs[Hb_x_nm1_gfn-1];
   Hb_x_np1  = gfs[Hb_x_np1_gfn-1];
   Hb_t_t  = gfs[Hb_t_t_gfn-1];
   Hb_t_t_n  = gfs[Hb_t_t_n_gfn-1];
   Hb_x_t  = gfs[Hb_x_t_gfn-1];
   Hb_x_t_n  = gfs[Hb_x_t_n_gfn-1];
   hb_t_res  = gfs[hb_t_res_gfn-1];
   hb_i_res  = gfs[hb_i_res_gfn-1];

   Hb_t_0  = gfs[Hb_t_0_gfn-1];
   Hb_x_0  = gfs[Hb_x_0_gfn-1];

   zetab     = gfs[zetab_gfn-1];
   zetab_lop = gfs[zetab_lop_gfn-1];
   zetab_res = gfs[zetab_res_gfn-1];
   zetab_rhs = gfs[zetab_rhs_gfn-1];

   w1 = gfs[w1_gfn-1];
   w2 = gfs[w2_gfn-1];
   w3 = gfs[w3_gfn-1];
   w4 = gfs[w4_gfn-1];

   mg_w1    =gfs[mg_w1_gfn-1]; 
   mg_w2    =gfs[mg_w2_gfn-1]; 
   mg_w3    =gfs[mg_w3_gfn-1]; 
   mg_w4    =gfs[mg_w4_gfn-1]; 

   efe_all_ires = gfs[efe_all_ires_gfn-1];
   iresall = gfs[iresall_gfn-1];
   phj_ires = gfs[phj_ires_gfn-1];
   phj = gfs[phj_gfn-1];
   pij_ires = gfs[pij_ires_gfn-1];
   pij = gfs[pij_gfn-1];
   alphaq_ires = gfs[alphaq_ires_gfn-1];
   alphaq = gfs[alphaq_gfn-1];
   thetap_ires = gfs[thetap_ires_gfn-1];
   thetap = gfs[thetap_gfn-1];

   mask    = gfs[mask_gfn-1];
   mask_mg = gfs[mask_mg_gfn-1];
   chr = gfs[chr_gfn-1]; 
   chr_mg = gfs[chr_mg_gfn-1]; 
}

//=============================================================================
// PAMR_get_dxdt() only works with AMR hierarchy levels ... here we use
// lambda for dt, but this only works if rhosp=rhotm
//=============================================================================
void ldptr_mg(void)
{
   real lambda;

   ldptr();

   dx=x[1]-x[0];
   PAMR_get_lambda(&lambda);
   dt=lambda*dx;
}

//=============================================================================
// utility routines
//=============================================================================
void const_f(real *f, real c)
{
   int i;

   for (i=0; i<Nx; i++) f[i]=c;
}

void zero_f(real *f)
{
   const_f(f,0);
}

void zero_f_ex(real *f, real *chr)
{
   int i;

   for (i=0; i<Nx; i++) if (chr[i]==AMRD_ex) f[i]=0;
}

real norm_l2(real *f, real *cmask, real *chr)
{
   int i;
   real norm=0;
   int sum=0;

   for (i=0; i<Nx; i++) 
      if (cmask[i]==AMRD_CMASK_ON && (chr[i]!=AMRD_ex)) { sum++; norm+=f[i]*f[i]; }

   if (!sum) sum=1;
   return (sqrt(norm/sum));
}

real ex_ahr(real xh)
{
   return xh*(1-ex_rbuf);
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

   cp_version=ADS5DP_CP_VERSION;
   AMRD_int_param(pfile,"cp_version",&cp_version,1);
   
   AdS_L=1.0; AMRD_real_param(pfile,"AdS_L",&AdS_L,1);

   interptype=0; AMRD_int_param(pfile,"interptype",&interptype,1);
   i_shift=0; AMRD_int_param(pfile,"i_shift",&i_shift,1);

   lxh=0;
   xh=0;
   ex_rbuf=0.1; AMRD_real_param(pfile,"ex_rbuf",&ex_rbuf,1);

   return;
}

void CFE4D_var_post_init(char *pfile)
{
   if (my_rank==0)
   {
      printf("===================================================================\n");
      printf("Reading CFE4D parameters:\n\n");
   }

   phi1_amp_1=phi1_r0_1=phi1_x0_1[0]=phi1_amp2_1=phi1_r02_1=phi1_x02_1[0]=0;
   phi1_amp_2=phi1_r0_2=phi1_x0_2[0]=phi1_amp2_2=phi1_r02_2=phi1_x02_2[0]=0;
   phi1_amp_3=phi1_r0_3=phi1_x0_3[0]=phi1_amp2_3=phi1_r02_3=phi1_x02_3[0]=0;
   phi1_width_1[0]=phi1_width2_1[0]=1;
   phi1_width_2[0]=phi1_width2_2[0]=1;
   phi1_width_3[0]=phi1_width2_3[0]=1;
   phi1_jj_1=phi1_jj2_1=0;
   phi1_jj_2=phi1_jj2_2=0;
   phi1_jj_3=phi1_jj2_3=0;

   AMRD_real_param(pfile,"phi1_amp_1",&phi1_amp_1,1);
   AMRD_real_param(pfile,"phi1_r0_1",&phi1_r0_1,1);
   AMRD_real_param(pfile,"phi1_delta_1",&phi1_delta_1,1);
   AMRD_real_param(pfile,"phi1_x0_1",phi1_x0_1,AMRD_dim);
   AMRD_real_param(pfile,"phi1_width_1",phi1_width_1,AMRD_dim);
   AMRD_int_param(pfile,"phi1_jj_1",&phi1_jj_1,1);
   AMRD_real_param(pfile,"phi1_amp2_1",&phi1_amp2_1,1);
   AMRD_real_param(pfile,"phi1_r02_1",&phi1_r02_1,1);
   AMRD_real_param(pfile,"phi1_delta2_1",&phi1_delta2_1,1);
   AMRD_real_param(pfile,"phi1_x02_1",phi1_x02_1,AMRD_dim);
   AMRD_real_param(pfile,"phi1_width2_1",phi1_width2_1,AMRD_dim);
   AMRD_int_param(pfile,"phi1_jj2_1",&phi1_jj2_1,1);

   AMRD_real_param(pfile,"phi1_amp_2",&phi1_amp_2,1);
   AMRD_real_param(pfile,"phi1_r0_2",&phi1_r0_2,1);
   AMRD_real_param(pfile,"phi1_delta_2",&phi1_delta_2,1);
   AMRD_real_param(pfile,"phi1_x0_2",phi1_x0_2,AMRD_dim);
   AMRD_real_param(pfile,"phi1_width_2",phi1_width_2,AMRD_dim);
   AMRD_int_param(pfile,"phi1_jj_2",&phi1_jj_2,1);
   AMRD_real_param(pfile,"phi1_amp2_2",&phi1_amp2_2,1);
   AMRD_real_param(pfile,"phi1_r02_2",&phi1_r02_2,1);
   AMRD_real_param(pfile,"phi1_delta2_2",&phi1_delta2_2,1);
   AMRD_real_param(pfile,"phi1_x02_2",phi1_x02_2,AMRD_dim);
   AMRD_real_param(pfile,"phi1_width2_2",phi1_width2_2,AMRD_dim);
   AMRD_int_param(pfile,"phi1_jj2_2",&phi1_jj2_2,1);

   AMRD_real_param(pfile,"phi1_amp_3",&phi1_amp_3,1);
   AMRD_real_param(pfile,"phi1_r0_3",&phi1_r0_3,1);
   AMRD_real_param(pfile,"phi1_delta_3",&phi1_delta_3,1);
   AMRD_real_param(pfile,"phi1_x0_3",phi1_x0_3,AMRD_dim);
   AMRD_real_param(pfile,"phi1_width_3",phi1_width_3,AMRD_dim);
   AMRD_int_param(pfile,"phi1_jj_3",&phi1_jj_3,1);
   AMRD_real_param(pfile,"phi1_amp2_3",&phi1_amp2_3,1);
   AMRD_real_param(pfile,"phi1_r02_3",&phi1_r02_3,1);
   AMRD_real_param(pfile,"phi1_delta2_3",&phi1_delta2_3,1);
   AMRD_real_param(pfile,"phi1_x02_3",phi1_x02_3,AMRD_dim);
   AMRD_real_param(pfile,"phi1_width2_3",phi1_width2_3,AMRD_dim);
   AMRD_int_param(pfile,"phi1_jj2_3",&phi1_jj2_3,1);

   kappa_cd=0; AMRD_real_param(pfile,"kappa_cd",&kappa_cd,1);
   rho_cd=0; AMRD_real_param(pfile,"rho_cd",&rho_cd,1);

   stype=0; AMRD_int_param(pfile,"stype",&stype,1);

   background=0; AMRD_int_param(pfile,"background",&background,1);
   skip_constraints=0; AMRD_int_param(pfile,"skip_constraints",&skip_constraints,1);

   gauge_t=0; AMRD_int_param(pfile,"gauge_t",&gauge_t,1);
   gauge_i=0; AMRD_int_param(pfile,"gauge_i",&gauge_i,1);
   rho1_t=1; AMRD_real_param(pfile,"rho1_t",&rho1_t,1);
   rho2_t=1; AMRD_real_param(pfile,"rho2_t",&rho2_t,1);
   rho3_t=1; AMRD_real_param(pfile,"rho3_t",&rho3_t,1);
   rho4_t=1; AMRD_real_param(pfile,"rho4_t",&rho4_t,1);
   xi1_t=1; AMRD_real_param(pfile,"xi1_t",&xi1_t,1);
   xi2_t=1; AMRD_real_param(pfile,"xi2_t",&xi2_t,1);
   cbulk_t=0; AMRD_real_param(pfile,"cbulk_t",&cbulk_t,1);
   rho1_i=1; AMRD_real_param(pfile,"rho1_i",&rho1_i,1);
   rho2_i=1; AMRD_real_param(pfile,"rho2_i",&rho2_i,1);
   rho3_i=1; AMRD_real_param(pfile,"rho3_i",&rho3_i,1);
   rho4_i=1; AMRD_real_param(pfile,"rho4_i",&rho4_i,1);
   xi1_i=1; AMRD_real_param(pfile,"xi1_i",&xi1_i,1);
   xi2_i=1; AMRD_real_param(pfile,"xi2_i",&xi2_i,1);
   cbulk_i=0; AMRD_real_param(pfile,"cbulk_i",&cbulk_i,1);

   rhoa=1; AMRD_real_param(pfile,"rhoa",&rhoa,1);
   rhob=1; AMRD_real_param(pfile,"rhob",&rhob,1);

   lqstt0   = malloc((AMRD_steps/AMRD_save_ivec0[3]+1)*sizeof(real));
   lqspsi0  = malloc((AMRD_steps/AMRD_save_ivec0[3]+1)*sizeof(real));
   lqsmass0 = malloc((AMRD_steps/AMRD_save_ivec0[3]+1)*sizeof(real));
   qstt0    = malloc((AMRD_steps/AMRD_save_ivec0[3]+1)*sizeof(real));
   qspsi0   = malloc((AMRD_steps/AMRD_save_ivec0[3]+1)*sizeof(real));
   qsmass0  = malloc((AMRD_steps/AMRD_save_ivec0[3]+1)*sizeof(real));

   ief_bh_r0=0; AMRD_real_param(pfile,"ief_bh_r0",&ief_bh_r0,1);

   real Mh,rh,xh;
   Mh=ief_bh_r0/2;
   rh=(pow(2,1.0/3.0)*pow(9*pow(AdS_L,2)*ief_bh_r0+sqrt(12*pow(AdS_L,6)+81*pow(AdS_L,4)*pow(ief_bh_r0,2)),2.0/3.0)-2*pow(3,1.0/3.0)*pow(AdS_L,2))/pow(6,2.0/3.0)/pow(9*pow(AdS_L,2)*ief_bh_r0+sqrt(12*pow(AdS_L,6)+81*pow(AdS_L,4)*pow(ief_bh_r0,2)),1.0/3.0);
   if (my_rank==0 && ief_bh_r0!=0) printf("\nanalytic BH initial data\n"
         "BH parameter r0=%lf\n"
         "BH mass M=%lf\n"
         "BH radius xh=%lf (areal radius rh=%lf in uncompactified)\n"
         "excision radius xh*(1-ex_rbuf)=%lf\n",ief_bh_r0,Mh,(sqrt(1+4*rh*rh)-1)/(2*rh),rh,(sqrt(1+4*rh*rh)-1)/(2*rh)*(1-ex_rbuf));

   if (AMRD_do_ex==0) AMRD_stop("require excision to be on","");

   if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Sets all variables to their 'zero' values:
//=============================================================================
void CFE4D_AMRH_var_clear(void)
{
   ldptr();

   zero_f(phi1_n); 

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

   zero_f(phi1_t_n); // sets initial time derivatives for ID
   zero_f(gb_tt_t_n);
   zero_f(gb_tx_t_n);
   zero_f(gb_xx_t_n);
   zero_f(psi_t_n);
   zero_f(Hb_t_t_n);
   zero_f(Hb_x_t_n);

   gauss1d_(phi1_n,
            &phi1_amp_1,&phi1_r0_1,&phi1_delta_1,&phi1_x0_1[0],&phi1_width_1[0],
            &phi1_amp2_1,&phi1_r02_1,&phi1_delta2_1,&phi1_x02_1[0],&phi1_width2_1[0],
            &AdS_L,x,&Nx);

   gauss1d_(w1,
            &phi1_amp_2,&phi1_r0_2,&phi1_delta_2,&phi1_x0_2[0],&phi1_width_2[0],
            &phi1_amp2_2,&phi1_r02_2,&phi1_delta2_2,&phi1_x02_2[0],&phi1_width2_2[0],
            &AdS_L,x,&Nx);

   gauss1d_(w2,
            &phi1_amp_3,&phi1_r0_3,&phi1_delta_3,&phi1_x0_3[0],&phi1_width_3[0],
            &phi1_amp2_3,&phi1_r02_3,&phi1_delta2_3,&phi1_x02_3[0],&phi1_width2_3[0],
            &AdS_L,x,&Nx);

   for (i=0; i<size; i++) phi1_n[i]+=w1[i]+w2[i]; 

   return;
}  

//=============================================================================
// Initialize any "elliptic_vars_t0" post construction of MGH, but before
// the start of vcycling.
//=============================================================================
void CFE4D_elliptic_vars_t0_init(void)
{
   // initializes dt, dx
   ldptr_mg();

   // initializes zetab conformal factor variable to zero
   const_f(zetab,0);
}

//=============================================================================
// Initial constraint data --- called after each MG iteration.
//
// Here we also initialize past time level information if
// AMRD_id_pl_method==3
//
// NOTE: np1,n,nm1 variables are allocated only at the top level of the MG hierarchy, 
//       so do an if(f_nm1){...}, for example, to make sure we are at the top level
//=============================================================================
void CFE4D_t0_cnst_data(void)
{
   int i;

   ldptr_mg();

   // initialize gbars
   if (skip_constraints)
   {
     if (ief_bh_r0!=0)
     {
       init_sch_(gb_tt,gb_tx,gb_xx,psi,&ief_bh_r0,
                 &AdS_L,phys_bdy,x,chr_mg,&AMRD_ex,&Nx); 
     }
     else
     {
       init_ads_(gb_tt,gb_tx,gb_xx,psi,
                 &AdS_L,phys_bdy,x,chr_mg,&AMRD_ex,&Nx);
     }
   }
   else
   {
     init_ghb_(zetab,phi1,gb_tt,gb_tx,gb_xx,psi,
               &AdS_L,phys_bdy,x,chr_mg,&AMRD_ex,&Nx,&rhoa,&rhob);
   }

   // initialize nm1,np1 time levels and hbars
   if (AMRD_id_pl_method==3 && phi1_nm1)
   {

     init_hb_(gb_tt_np1,gb_tt_n,gb_tt_nm1,
              gb_tx_np1,gb_tx_n,gb_tx_nm1,
              gb_xx_np1,gb_xx_n,gb_xx_nm1,
              psi_np1,psi_n,psi_nm1,
              Hb_t_n,Hb_x_n,
              &AdS_L,phys_bdy,x,&dt,chr,&AMRD_ex,&Nx);

     init_nm1_(gb_tt_np1,gb_tt_n,gb_tt_nm1,gb_tt_t_n,
               gb_tx_np1,gb_tx_n,gb_tx_nm1,gb_tx_t_n,
               gb_xx_np1,gb_xx_n,gb_xx_nm1,gb_xx_t_n,
               psi_np1,psi_n,psi_nm1,psi_t_n,
               Hb_t_np1,Hb_t_n,Hb_t_nm1,Hb_t_t_n,
               Hb_x_np1,Hb_x_n,Hb_x_nm1,Hb_x_t_n,
               phi1_np1,phi1_n,phi1_nm1,phi1_t_n,
               &AdS_L,phys_bdy,x,&dt,chr,&AMRD_ex,&Nx);

     // store initial source functions, metric components
     for (i=0; i<size; i++)
     {
       Hb_t_0[i]=Hb_t[i];
       Hb_x_0[i]=Hb_x[i];
     }

   }

   return;
}

//=============================================================================
// Calculations prior to saving info to disk.
//
// NOTE: at this point, the time sequence is: n,nm1,np1 (unless t=0)
//=============================================================================
void CFE4D_pre_io_calc(void)
{
   ldptr();

   int i;
   real ct;

   dx=x[1]-x[0];
   ct=PAMR_get_time(g_L);

   // compute independent residuals of the AdS4D system
   if (ct!=0)
   {
     //(NOTE: for t>0, have cycled time sequence np1,n,nm1 to time sequence n,nm1,np1,
     // so here, time level n is the most advanced time level)
     ires_(efe_all_ires,phj_ires,
        pij_ires,alphaq_ires,thetap_ires,
        gb_tt_n,gb_tt_nm1,gb_tt_np1,
        gb_tx_n,gb_tx_nm1,gb_tx_np1,
        gb_xx_n,gb_xx_nm1,gb_xx_np1,
        psi_n,psi_nm1,psi_np1,
        Hb_t_np1,Hb_t_n,Hb_t_nm1,
        Hb_x_np1,Hb_x_n,Hb_x_nm1,
        phi1_n,phi1_nm1,phi1_np1,
        x,&dt,chr,&AdS_L,&AMRD_ex,&Nx,phys_bdy,ghost_width);
   }
   else
   {
     //(NOTE: for t=0, have *not* cycled time sequence, so still np1,n,nm1,
     // so here, time level np1 is the most advanced time level)
     ires_(efe_all_ires,phj_ires,
        pij_ires,alphaq_ires,thetap_ires,
        gb_tt_np1,gb_tt_n,gb_tt_nm1,
        gb_tx_np1,gb_tx_n,gb_tx_nm1,
        gb_xx_np1,gb_xx_n,gb_xx_nm1,
        psi_np1,psi_n,psi_nm1,
        Hb_t_np1,Hb_t_n,Hb_t_nm1,
        Hb_x_np1,Hb_x_n,Hb_x_nm1,
        phi1_np1,phi1_n,phi1_nm1,
        x,&dt,chr,&AdS_L,&AMRD_ex,&Nx,phys_bdy,ghost_width);
   }

   // fill in ires arrays with independent residuals
   for (i=0; i<Nx; i++)
   {
     // excise i=Nx-1 pts (pure AdS diverges at i=Nx, so cannot use these pts in difference stencils) 
     if (chr[i]==AMRD_ex 
         || x[i]<0.5*dx_Lc || 1-x[i]<0.5*dx_Lc) 
     {
       phj[i]=0;  
       iresall[i]=0;
       pij[i]=0;
       alphaq[i]=0;
       thetap[i]=0;
     }
     else
     {
       phj[i]=phj_ires[i];
       iresall[i]=efe_all_ires[i];
       pij[i]=pij_ires[i];
       alphaq[i]=alphaq_ires[i];
       thetap[i]=thetap_ires[i];
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
   real l2norm=0,l2norm_gb,l2norm_hb_t,l2norm_hb_i,l2norm_phi1;
   int is_nan;

   ldptr();

   if (LIN_ZERO_BND) 
   {
      lin_zero_bnd_res_(gb_res,phys_bdy,&lin_zero_bnd_all,&Nx);
      lin_zero_bnd_res_(hb_t_res,phys_bdy,&lin_zero_bnd_all,&Nx);
      lin_zero_bnd_res_(hb_i_res,phys_bdy,&lin_zero_bnd_all,&Nx);
      lin_zero_bnd_res_(kg_res,phys_bdy,&lin_zero_bnd_all,&Nx);
   }

   l2norm_gb  =norm_l2(gb_res,mask,chr);
   l2norm_hb_t=norm_l2(hb_t_res,mask,chr);
   l2norm_hb_i=norm_l2(hb_i_res,mask,chr);
   l2norm_phi1=norm_l2(kg_res,mask,chr);

   l2norm=l2norm_gb+l2norm_hb_t+l2norm_hb_i+l2norm_phi1;

   check_nan_(&l2norm,&is_nan);

   if (is_nan)
   {
      printf("\nl2norm_gb=%lf, l2norm_hb_t=%lf, l2norm_hb_i=%lf, l2norm_phi1=%lf\n",
              l2norm_gb,l2norm_hb_t,l2norm_hb_i,l2norm_phi1);
      printf("Nx=%i,L=%i\n",Nx,g_L);
      AMRD_stop("l2norm is nan ... stopping","");
      l2norm=0;
   }

   return l2norm;
}

//=============================================================================
// Performs 1 iteration of the evolution equations 
//
// NOTE: at this point, the time sequence is: np1,n,nm1 
//=============================================================================
void CFE4D_evolve(int iter)
{
   int i;
   int zero_i=0;
   int ltrace=0;
   real ct,zero=0;

   ldptr();

   ct=PAMR_get_time(g_L);

   if (ltrace) printf("CFE4D_evolve: iter=%i , time=%lf, lev=%i, rank=%i\n",iter,ct,g_L,my_rank);

   hb_t_evo_(hb_t_res,
             gb_tt_np1,gb_tt_n,gb_tt_nm1,
             gb_tx_np1,gb_tx_n,gb_tx_nm1,
             gb_xx_np1,gb_xx_n,gb_xx_nm1,
             psi_np1,psi_n,psi_nm1,
             Hb_t_np1,Hb_t_n,Hb_t_nm1,
             Hb_x_np1,Hb_x_n,Hb_x_nm1,
             phi1_np1,phi1_n,phi1_nm1,
             &AdS_L,x,&dt,chr,&AMRD_ex,
             phys_bdy,ghost_width,&Nx,
             Hb_t_0,Hb_x_0,
             &gauge_t,&ct,&rho1_t,&rho2_t,&rho3_t,&rho4_t,&xi1_t,&xi2_t,&cbulk_t);

   hb_i_evo_(hb_i_res,
             gb_tt_np1,gb_tt_n,gb_tt_nm1,
             gb_tx_np1,gb_tx_n,gb_tx_nm1,
             gb_xx_np1,gb_xx_n,gb_xx_nm1,
             psi_np1,psi_n,psi_nm1,
             Hb_t_np1,Hb_t_n,Hb_t_nm1,
             Hb_x_np1,Hb_x_n,Hb_x_nm1,
             phi1_np1,phi1_n,phi1_nm1,
             &AdS_L,x,&dt,chr,&AMRD_ex,
             phys_bdy,ghost_width,&Nx,
             Hb_t_0,Hb_x_0,
             &gauge_i,&ct,&rho1_i,&rho2_i,&rho3_i,&rho4_i,&xi1_i,&xi2_i,&cbulk_i);

   g_evo_opt_(gb_res,kg_res,cl_res,
              gb_tt_np1,gb_tt_n,gb_tt_nm1,
              gb_tx_np1,gb_tx_n,gb_tx_nm1,
              gb_xx_np1,gb_xx_n,gb_xx_nm1,
              psi_np1,psi_n,psi_nm1,
              Hb_t_np1,Hb_t_n,Hb_t_nm1,
              Hb_x_np1,Hb_x_n,Hb_x_nm1,
              phi1_np1,phi1_n,phi1_nm1,
              &AdS_L,x,&dt,chr,&AMRD_ex,
              phys_bdy,ghost_width,&Nx,
              &background,&kappa_cd,&rho_cd);

   return;
}

//=============================================================================
// sets excision mask (NO ITERATOR, SO DON'T LOAD POINTERS!!!)
// 
// in polar code, only excised regions are those inside horizons
//=============================================================================
void CFE4D_fill_ex_mask(real *mask, int dim, int *shape, real *bbox, real excised)
{
   int i;
   int l,n,a,b;
   int axiregpts;
   int rho_sp[PAMR_MAX_LEVS],rho_tm[PAMR_MAX_LEVS];
   int Lf,rho_sp0,rho_tm0;
   real x,dx,ct;

   PAMR_get_rho(rho_sp,rho_tm,PAMR_MAX_LEVS);
   Lf=PAMR_get_max_lev(PAMR_AMRH);
   rho_sp0=rho_sp[Lf]; rho_tm0=rho_tm[Lf];

   dx=(bbox[1]-bbox[0])/(shape[0]-1);

   for (i=0; i<shape[0]; i++)
   {
      x=bbox[0]+i*dx;
      mask[i]=excised-1;
      if (xh!=0) // only activates when AH found
      {
         if (x/ex_ahr(xh)<1) mask[i]=excised;
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
int AH_count[MAX_BHS],found_AH[MAX_BHS],freq0[MAX_BHS];
int pre_tstep_global_first=1,valid;

void CFE4D_pre_tstep(int L)
{
   char name[256];
   int AH[MAX_BHS];
   int AH_shape[1],got_an_AH,do_reinit_ex,do_repop;
   real AH_bbox[2],AH_min_resid0;

   real ct;
   real new_rbuf;
   int n,i,Lf,Lc;

   int is,ie;
  
   ct=PAMR_get_time(L);

   Lf=PAMR_get_max_lev(PAMR_AMRH);
   Lc=PAMR_get_min_lev(PAMR_AMRH);  //if (PAMR_get_max_lev(PAMR_AMRH)>1) Lc=2; else Lc=1;

   // horizon quantities, when at coarsest level L=Lc
   if (L==Lc)
   {
     valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
     while(valid)
     {
       ldptr();

       // find radius of outermost trapped surface
       if (ct!=0)
       {
         //(NOTE: for t>0, have cycled time sequence np1,n,nm1 to time sequence n,nm1,np1,
         // so here, time level n is the most advanced time level)
         ires_(efe_all_ires,phj_ires,
            pij_ires,alphaq_ires,thetap_ires,
            gb_tt_n,gb_tt_nm1,gb_tt_np1,
            gb_tx_n,gb_tx_nm1,gb_tx_np1,
            gb_xx_n,gb_xx_nm1,gb_xx_np1,
            psi_n,psi_nm1,psi_np1,
            Hb_t_np1,Hb_t_n,Hb_t_nm1,
            Hb_x_np1,Hb_x_n,Hb_x_nm1,
            phi1_n,phi1_nm1,phi1_np1,
            x,&dt,chr,&AdS_L,&AMRD_ex,&Nx,phys_bdy,ghost_width);
       }
       else
       {
         //(NOTE: for t=0, have *not* cycled time sequence, so still np1,n,nm1,
         // so here, time level np1 is the most advanced time level)
         ires_(efe_all_ires,phj_ires,
            pij_ires,alphaq_ires,thetap_ires,
            gb_tt_np1,gb_tt_n,gb_tt_nm1,
            gb_tx_np1,gb_tx_n,gb_tx_nm1,
            gb_xx_np1,gb_xx_n,gb_xx_nm1,
            psi_np1,psi_n,psi_nm1,
            Hb_t_np1,Hb_t_n,Hb_t_nm1,
            Hb_x_np1,Hb_x_n,Hb_x_nm1,
            phi1_np1,phi1_n,phi1_nm1,
            x,&dt,chr,&AdS_L,&AMRD_ex,&Nx,phys_bdy,ghost_width);
       }
       for (i=0; i<Nx; i++)
       {
         if(thetap_ires[i]<0 && x[i]>lxh) {lxh=x[i];}
       }

       valid=PAMR_next_g();
     }
     MPI_Allreduce(&lxh,&xh,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

     // re-initialize mask function
     if (xh!=0) {PAMR_excision_on("chr",&CFE4D_fill_ex_mask,AMRD_ex,1);}

     // calculates horizon mass and horizon areal radius
     int ih;
     real xmin,xmax;
     real lM0,M0,lR0,R0;
     real larea,lgchch,lgthth;
     xmin=x[0]; xmax=x[Nx-1];
     lM0=0; M0=0; lR0=0; R0=0;
     larea=0; lgchch=0; lgthth=0;
     if (xh>xmin && xh<xmax) //only considers the processor whose grid contain the horizon, so ih index calculation makes sense
     {
        ih=(xh-xmin)/(x[1]-x[0]);

        valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
        while(valid)
        {
//           lgchch=x[ih]*x[ih]/pow(1-x[ih]*x[ih],2)+psi_n[ih]*x[ih]*x[ih];
//           lgthth=x[ih]*x[ih]/pow(1-x[ih]*x[ih],2)+psi_n[ih]*x[ih]*x[ih];
//           larea=4*M_PI*sqrt(lgchch*lgthth);
//           lR0=sqrt(larea/4/M_PI);
   
           lR0=x[ih]*sqrt(1/(1-x[ih]*x[ih])/(1-x[ih]*x[ih])+psi_n[ih]); //equivalent to above four lines
   
           lM0=lR0*(1+lR0*lR0/AdS_L/AdS_L)/2;
   
           valid=PAMR_next_g();
        }
     } 
     MPI_Allreduce(&lM0,&M0,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
     MPI_Allreduce(&lR0,&R0,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
     
     // displays horizon mass and horizon areal radius
     if (my_rank==0)
     {
        if (xh!=0) 
        {
           printf("\n ... found an AH ... \n");
           printf("     horizon mass: %5.3lf, from horizon radius xh=%5.3lf (areal radius rh=%5.3lf)\n:",M0,xh,R0);
        }
     }
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
   double lE4,gE4,liphi4,giphi4,cE4;
   char out1_name[256];
   static FILE *out1;
   static int local_first = 1;

   real ct;
   int n,i,Lf,Lc;

   int is,ie;

   ct=PAMR_get_time(L);

   Lf=PAMR_get_max_lev(PAMR_AMRH);  
   Lc=PAMR_get_min_lev(PAMR_AMRH);  //if (PAMR_get_max_lev(PAMR_AMRH)>1) Lc=2; else Lc=1;

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

   // solves for zetab conformal factor at t=0; residual
   mg_sup_(&action,zetab,zetab_rhs,zetab_lop,zetab_res,
           phi1,&AdS_L,mask_mg,phys_bdy,chr_mg,&AMRD_ex,x,&norm,&Nx);

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

   // solves for zetab conformal factor at t=0; residual
   mg_sup_(&action,zetab,zetab_rhs,zetab_lop,zetab_res,
           phi1,&AdS_L,mask_mg,phys_bdy,chr_mg,&AMRD_ex,x,&norm,&Nx);

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

   // solves for zetab conformal factor at t=0; residual
   mg_sup_(&action,zetab,zetab_rhs,zetab_lop,zetab_res,
           phi1,&AdS_L,mask_mg,phys_bdy,chr_mg,&AMRD_ex,x,&norm,&Nx);

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
}

//=============================================================================
//check-pointing
//=============================================================================
#define CP_DATA_SIZE 50000
void CFE4D_copy_block(char *r, char **q, int n, int dir, int *tot_size)
{
   char *r0;
   int n0;

   if (n==0) return;

   if ((*tot_size+n) > (CP_DATA_SIZE))
      AMRD_stop("CFE4D_copy_block: error ... CP_DATA_SIZE too small\n","");
   *tot_size+=n;

   n0=n;
   r0=r;

   if (dir==AMRD_CP_SAVE) while(n0--) *(*q)++=*r0++;
   else while(n0--) *r0++=*(*q)++;

   return;
}

void CFE4D_cp(int dir, char *data)
{
   int size=0;
   
   if (dir==AMRD_CP_SAVE)
   {
       cp_version=ADS5DP_CP_VERSION;
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

