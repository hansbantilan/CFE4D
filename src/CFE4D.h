#ifndef _ADS5D_H
#define _ADS5D_H

/*============================================================================= */
/* global variables and prototypes for CFE4D                                    */
/*============================================================================= */

extern int axisym;
extern real z_axysim;

extern real AdS_L;

extern real kappa_cd,rho_cd;

extern int post_tstep_global_first;

extern real diss_eps_k,diss_eps_y_cutoff;
extern int diss_kmax,diss_eps_k_cutoff_n,diss_bdy_k,diss_all_past_k,diss_all;

/*============================================================================= */
/* id parameters                                                                */
/*============================================================================= */

//TEST//
#define MAX_BHS 4

extern real phi1_amp_1,phi1_r0_1,phi1_delta_1,phi1_x0_1[2],phi1_ecc_1[2];
extern real phi1_amp_2,phi1_r0_2,phi1_delta_2,phi1_x0_2[2],phi1_ecc_2[2];

extern real phi4_qb_amp_1,phi4_qb_x0_1[2],phi4_qb_w0_1,phi4_qb_k_1[2],phi4_qb_delta_1;
extern real phi4_qb_amp_2,phi4_qb_x0_2[2],phi4_qb_w0_2,phi4_qb_k_2[2],phi4_qb_delta_2;

extern real ief_bh_r0;
extern real ex_rbuf[MAX_BHS];
extern int ex_reset_rbuf;
extern real ex_r[MAX_BHS][2],ex_xc[MAX_BHS][2];

extern int bh_background,background,skip_constraints;
extern int use_bg_inc,output_ires,output_quasiset;
extern int which_gb_ii_res,which_gb_ii_delta_ads; 

extern int interptype;
extern int i_shift;
extern int regtype;
extern int stype;
 
extern int harmonize;

/* gauge parameters: */
extern int gauge_t;
extern int gauge_i;
extern real c1_t,c2_t,c3_t;
extern real c1_i,c2_i,c3_i;
extern real rho1_t,rho2_t,xi1_t,xi2_t;
extern real rho1_i,rho2_i,xi1_i,xi2_i;



/* TTEST */
/* apph parameters  */
extern int AH_Nchi[MAX_BHS],AH_Nphi[MAX_BHS],AH_Lmin[MAX_BHS],AH_Lmax[MAX_BHS],AH_find_best_fit[MAX_BHS];
extern int AH_max_iter[MAX_BHS],AH_freq[MAX_BHS],AH_freq_aft[MAX_BHS],AH_rsteps[MAX_BHS],AH_maxinc[MAX_BHS];
extern real AH_tol[MAX_BHS],AH_tol_aft[MAX_BHS],AH_r0[MAX_BHS],AH_lambda[MAX_BHS],AH_lambda_min[MAX_BHS];
extern real AH_eps[MAX_BHS],AH_r1[MAX_BHS],AH_tol_scale[MAX_BHS],AH_reset_scale[MAX_BHS];
extern real AH_xc[MAX_BHS][2],AH_max_tol_inc[MAX_BHS],AH_tmin[MAX_BHS],AH_omt_scale[MAX_BHS];
extern int use_AH_new_smooth,use_AH_new;
extern int c_AH;

extern real *AH_theta[MAX_BHS],*AH_R[MAX_BHS],*AH_w1[MAX_BHS],*AH_w2[MAX_BHS],*AH_w3[MAX_BHS];
extern real *AH_theta_ads[MAX_BHS];
extern int *AH_lev[MAX_BHS],*AH_own[MAX_BHS];

extern int skip_ires;

extern int out1_freq;

extern int AH_count[MAX_BHS],found_AH[MAX_BHS],freq0[MAX_BHS];



#define ADS5D_CP_VERSION 1
extern int cp_version;

extern real *cl_res;

extern real *phi1,*phi1_n,*phi1_np1,*phi1_nm1;
extern real *phi4_r,*phi4_r_n,*phi4_r_np1,*phi4_r_nm1;
extern real *phi4_i,*phi4_i_n,*phi4_i_np1,*phi4_i_nm1;

extern real *phi4_r,*phi4_r_n,*phi4_r_np1,*phi4_r_nm1,*phi4_r_t;
extern real *phi4_i,*phi4_i_n,*phi4_i_np1,*phi4_i_nm1,*phi4_i_t;

extern real *gb_tt,*gb_tt_n,*gb_tt_np1,*gb_tt_nm1; 
extern real *gb_tx,*gb_tx_n,*gb_tx_np1,*gb_tx_nm1; 
extern real *gb_ty,*gb_ty_n,*gb_ty_np1,*gb_ty_nm1; 
extern real *gb_xx,*gb_xx_n,*gb_xx_np1,*gb_xx_nm1; 
extern real *gb_xy,*gb_xy_n,*gb_xy_np1,*gb_xy_nm1; 
extern real *gb_yy,*gb_yy_n,*gb_yy_np1,*gb_yy_nm1; 
extern real *psi,*psi_n,*psi_np1,*psi_nm1; 

extern real *Hb_t,*Hb_t_n,*Hb_t_np1,*Hb_t_nm1;
extern real *Hb_x,*Hb_x_n,*Hb_x_np1,*Hb_x_nm1;
extern real *Hb_y,*Hb_y_n,*Hb_y_np1,*Hb_y_nm1;

extern real *phi1_t,*phi1_t_n;
extern real *phi4_r_t,*phi4_r_t_n;
extern real *phi4_i_t,*phi4_i_t_n;
extern real *gb_tt_t,*gb_tt_t_n;
extern real *gb_tx_t,*gb_tx_t_n;
extern real *gb_ty_t,*gb_ty_t_n;
extern real *gb_xx_t,*gb_xx_t_n;
extern real *gb_xy_t,*gb_xy_t_n;
extern real *gb_yy_t,*gb_yy_t_n;
extern real *psi_t,*psi_t_n;
extern real *Hb_t_t,*Hb_t_t_n;
extern real *Hb_x_t,*Hb_x_t_n;
extern real *Hb_y_t,*Hb_y_t_n;

extern real *w1,*mg_w1;
extern real *w2,*mg_w2;
extern real *w3,*mg_w3;
extern real *w4,*mg_w4;

extern real *mask,*mask_mg,*chr,*chr_mg,AMRD_ex;
extern real *gu_tt,*gu_tx,*gu_ty,*gu_xx,*gu_xy,*gu_yy,*m_g_det;
extern real *kg_ires,*alpha,*theta,*f,*K;

extern real *phi1_res,*phi4_res,*gb_res,*gb_ii_res;
extern real *efe_tt_ires,*efe_tx_ires,*efe_ty_ires;
extern real *efe_xx_ires,*efe_xy_ires,*efe_yy_ires;
extern real *efe_psi_ires;
extern real *quasiset_tt,*quasiset_tx,*quasiset_ty;
extern real *quasiset_xx,*quasiset_xy,*quasiset_yy;
extern real *quasiset_psi;
extern real *quasiset_mass;

extern real *einstein_ll, *set_ll, *gamma_ull, *ricci_ll, *ricci;
extern real *phi1_x, *phi4_r_x, *phi4_i_x, *g_ll, *g_uu;
extern real *gfull_uu;

extern real *dm,*dmtest,*rho_in;

extern real *tfunction;
extern real *chrfunction;
extern real *irestt;
extern real *irestx;
extern real *iresty;
extern real *iresxx;
extern real *iresxy;
extern real *iresyy;
extern real *irespsi;
extern real *qstt;
extern real *qstx;
extern real *qsty;
extern real *qsxx;
extern real *qsxy;
extern real *qsyy;
extern real *qspsi;

extern real *zeta,*zeta_res,*zeta_lop,*zeta_rhs;

extern real *hb_t_res,*hb_i_res;
extern real *Hb_t_0,*Hb_x_0,*Hb_y_0;

extern real *g_norms;

extern real *x,*y,*z;
extern int shape[2],ghost_width[4],Nx,Ny,Nz,phys_bdy[4],size,g_rank;
extern real base_bbox[4],bbox[4],dx,dy,dz,dt,dx_Lc;
extern int g_L;

extern int cl_res_gfn;

extern int phi1_gfn,phi1_n_gfn,phi1_np1_gfn,phi1_nm1_gfn; 
extern int phi4_r_gfn,phi4_r_n_gfn,phi4_r_np1_gfn,phi4_r_nm1_gfn; 
extern int phi4_i_gfn,phi4_i_n_gfn,phi4_i_np1_gfn,phi4_i_nm1_gfn; 

extern int gb_tt_gfn,gb_tt_n_gfn,gb_tt_np1_gfn,gb_tt_nm1_gfn; 
extern int gb_tx_gfn,gb_tx_n_gfn,gb_tx_np1_gfn,gb_tx_nm1_gfn; 
extern int gb_ty_gfn,gb_ty_n_gfn,gb_ty_np1_gfn,gb_ty_nm1_gfn; 
extern int gb_xx_gfn,gb_xx_n_gfn,gb_xx_np1_gfn,gb_xx_nm1_gfn; 
extern int gb_xy_gfn,gb_xy_n_gfn,gb_xy_np1_gfn,gb_xy_nm1_gfn; 
extern int gb_yy_gfn,gb_yy_n_gfn,gb_yy_np1_gfn,gb_yy_nm1_gfn; 
extern int psi_gfn,psi_n_gfn,psi_np1_gfn,psi_nm1_gfn; 

extern int Hb_t_gfn,Hb_t_n_gfn,Hb_t_np1_gfn,Hb_t_nm1_gfn;
extern int Hb_x_gfn,Hb_x_n_gfn,Hb_x_np1_gfn,Hb_x_nm1_gfn;
extern int Hb_y_gfn,Hb_y_n_gfn,Hb_y_np1_gfn,Hb_y_nm1_gfn;

extern int mask_gfn,mask_mg_gfn,chr_gfn,chr_mg_gfn;
extern int kg_ires_gfn,alpha_gfn,theta_gfn;

extern int phi1_res_gfn,phi4_res_gfn,gb_res_gfn,gb_ii_res_gfn;
extern int phi4_r_res_gfn,phi4_i_res_gfn;
extern int efe_tt_ires_gfn,efe_xx_ires_gfn,efe_xy_ires_gfn,efe_yy_ires_gfn;
extern int efe_psi_ires_gfn;

extern int tfunction_gfn;
extern int chrfunction_gfn;
extern int irestt_gfn;
extern int irestx_gfn;
extern int iresty_gfn;
extern int iresxx_gfn;
extern int iresxy_gfn;
extern int iresyy_gfn;
extern int irespsi_gfn;
extern int qstt_gfn;
extern int qstx_gfn;
extern int qsty_gfn;
extern int qsxx_gfn;
extern int qsxy_gfn;
extern int qsyy_gfn;
extern int qspsi_gfn;

extern int zeta_gfn,zeta_res_gfn,zeta_lop_gfn,zeta_rhs_gfn;

extern int hb_t_res_gfn,hb_i_res_gfn;

extern int w1_gfn,mg_w1_gfn;
extern int w2_gfn,mg_w2_gfn;
extern int w3_gfn,mg_w3_gfn;
extern int w4_gfn,mg_w4_gfn;

#define MAX_N 2049
extern real csr[MAX_N],snr[MAX_N];

extern int skip_ires;
extern int skip_exp;

void set_gfns(void);
void ldptr_bbox(void);
void ldptr(void);
void ldptr_mg(void);
void const_f(real *f, real c);
void zero_f(real *f);
real norm_l2(real *f, real *cmask, real *chr);
void calc_gbu(void);
void calc_gbu_nm1(void);

/* prototypes for the various fortran functions we use: */
void g_evo_opt_(real *gb_res, real *phi1_res, real *cl_res,
                real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1,
                real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1,
                real *gb_ty_np1, real *gb_ty_n, real *gb_ty_nm1,
                real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1,
                real *gb_xy_np1, real *gb_xy_n, real *gb_xy_nm1,
                real *gb_yy_np1, real *gb_yy_n, real *gb_yy_nm1,
                real *psi_np1, real *psi_n, real *psi_nm1,
                real *Hb_t_np1, real *Hb_t_n, real *Hb_t_nm1,
                real *Hb_x_np1, real *Hb_x_n, real *Hb_x_nm1,
                real *Hb_y_np1, real *Hb_y_n, real *Hb_y_nm1,
                real *phi1_np1, real *phi1_n, real *phi1_nm1,
                real *AdS_L, real *x, real *y, real *dt, real *chr, real *ex,
                int *phys_bdy, int *ghost_width, int *Nx, int *Ny,
                int *background, real *kappa_cd, real *rho_cd,
                int *interptype, int *i_shift, int *regtype,
                int *diss_kmax, real *tfunction);

void init_ads5d_bh_(real *ief_bh_r0, real *AdS_L, real *gb_tt, real *gb_tx, 
                  real *gb_ty, real *gb_xx, real *gb_xy, real *gb_yy, real *psi,
                  real *gb_tt_t, real *gb_tx_t, real *gb_ty_t, 
                  real *gb_xx_t, real *gb_xy_t, real *gb_yy_t, real *psi_t,
                  real *Hb_t, real *Hb_x, real *Hb_y,
                  real *Hb_t_t, real *Hb_x_t, real *Hb_y_t,
                  int *phys_bdy, real *x, real *y, real *dt, real *chr, 
                  real *ex, int *Nx, int *Ny, int *regtype);

void init_ghb_ads_(real *gb_tt, real *gb_tx, real *gb_ty,
                   real *gb_xx, real *gb_xy, real *gb_yy, real *psi,
                   real *Hb_t, real *Hb_x, real *Hb_y, real *AdS_L, 
                   real *x, real *y, real *chr, real *ex, int *Nx, int *Ny,
                   int *regtype);

void hb_t_evo_(real *res,
               real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1,
               real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1,
               real *gb_ty_np1, real *gb_ty_n, real *gb_ty_nm1,
               real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1,
               real *gb_xy_np1, real *gb_xy_n, real *gb_xy_nm1,
               real *gb_yy_np1, real *gb_yy_n, real *gb_yy_nm1,
               real *psi_np1, real *psi_n, real *psi_nm1,
               real *Hb_t_np1, real *Hb_t_n, real *Hb_t_nm1,
               real *Hb_x_np1, real *Hb_x_n, real *Hb_x_nm1,
               real *Hb_y_np1, real *Hb_y_n, real *Hb_y_nm1,
               real *phi1_np1, real *phi1_n, real *phi1_nm1,
               real *AdS_L, real *x, real *y, real *dt, real *chr, real *ex, 
               int *phys_bdy, int *ghost_width, int *Nx, int *Ny, 
               real *Hb_t_0, real *Hb_x_0, real *Hb_y_0,
               int *gauge, real *t_n, real *rho1_t, real *rho2_t, real *rho3, real *rho4, real *xi1_t, real *xi2_t, 
               real *c1, real *c2, real *c3, real *cbulk);

void hb_i_evo_(real *res,
               real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1,
               real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1,
               real *gb_ty_np1, real *gb_ty_n, real *gb_ty_nm1,
               real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1,
               real *gb_xy_np1, real *gb_xy_n, real *gb_xy_nm1,
               real *gb_yy_np1, real *gb_yy_n, real *gb_yy_nm1,
               real *psi_np1, real *psi_n, real *psi_nm1,
               real *Hb_t_np1, real *Hb_t_n, real *Hb_t_nm1,
               real *Hb_x_np1, real *Hb_x_n, real *Hb_x_nm1,
               real *Hb_y_np1, real *Hb_y_n, real *Hb_y_nm1,
               real *phi1_np1, real *phi1_n, real *phi1_nm1,
               real *AdS_L, real *x, real *y, real *dt, real *chr, real *ex, 
               int *phys_bdy, int *ghost_width, int *Nx, int *Ny, 
               real *Hb_t_0, real *Hb_x_0, real *Hb_y_0,
               int *gauge, real *t_n, real *rho1_i, real *rho2_i, real *rho3, real *rho4, real *xi1_i, real *xi2_i, 
               real *c1, real *c2, real *c3, real *cbulk);

void ires_(real *efe_all_ires,
           real *efe_tt_ires, real *efe_tx_ires, real *efe_ty_ires,
           real *efe_xx_ires, real *efe_xy_ires, real *efe_yy_ires,
           real *efe_psi_ires,
           real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1,
           real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1,
           real *gb_ty_np1, real *gb_ty_n, real *gb_ty_nm1,
           real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1,
           real *gb_xy_np1, real *gb_xy_n, real *gb_xy_nm1,
           real *gb_yy_np1, real *gb_yy_n, real *gb_yy_nm1,
           real *psi_np1, real *psi_n, real *psi_nm1,
           real *Hb_t_np1, real *Hb_t_n, real *Hb_t_nm1,
           real *Hb_x_np1, real *Hb_x_n, real *Hb_x_nm1,
           real *Hb_y_np1, real *Hb_y_n, real *Hb_y_nm1,
           real *phi1_np1, real *phi1_n, real *phi1_nm1,
           real *x, real *y, real *dt, real *chr, 
           real *AdS_L, real *ex, int *Nx, int *Ny, int *phys_bdy, int *ghost_width);

void gu_calc_(real *gb_tt, real *gb_tx, real *gb_ty, real *gb_xx, 
               real *gb_xy, real *gb_yy, real *gb_psi, 
              real *gu_tt, real *gu_tx, real *gu_ty, real *gu_xx,
              real *gu_xy, real *gu_yy, 
              real *m_g_det, real *chr, real *ex, int *Nx, int *Ny);

void mg_sup_(int *action, real *zeta, real *zeta_rhs, real *zeta_lop, real *zeta_res, real *phi1, 
             real *AdS_L, real *cmask, int *phys_bdy, real *chr, real *ex, 
             real *x, real *y, real *norm, int *Nx, int *Ny);

void init_ghb_(real *zeta,
               real *gb_tt, real *gb_tx, real *gb_ty, real *gb_xy,
               real *gb_xx, real *gb_yy, real *psi, 
               real *Hb_t, real *Hb_x, real *Hb_y, 
               real *AdS_L, real *cmask, int *phys_bdy, real *x, real *y, 
               real *chr, real *ex, int *Nx, int *Ny, int *regtype,
               real* rhoa, real* rhob);

void gauss2d_(real *f, real *A, real *B, real *r0, real *delta, real *xu0, real *yu0, real *ex, real *ey,
              real *AdS_L, real *x, real *y, int *Nx, int *Ny, int *stype);

void CFE4D_fill_ex_mask(real *mask, int dim, int *shape, real *bbox, real excised);

void init_nm1_(real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1, real *gb_tt_t_n,
               real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1, real *gb_tx_t_n,
               real *gb_ty_np1, real *gb_ty_n, real *gb_ty_nm1, real *gb_ty_t_n,
               real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1, real *gb_xx_t_n,
               real *gb_xy_np1, real *gb_xy_n, real *gb_xy_nm1, real *gb_xy_t_n,
               real *gb_yy_np1, real *gb_yy_n, real *gb_yy_nm1, real *gb_yy_t_n,
               real *psi_np1, real *psi_n, real *psi_nm1, real *psi_t_n,
               real *Hb_t_np1, real *Hb_t_n, real *Hb_t_nm1, real *Hb_t_t_n,
               real *Hb_x_np1, real *Hb_x_n, real *Hb_x_nm1, real *Hb_x_t_n,
               real *Hb_y_np1, real *Hb_y_n, real *Hb_y_nm1, real *Hb_y_t_n,
               real *phi1_np1, real *phi1_n, real *phi1_nm1, real *phi1_t_n, real *tfunction,
               real *AdS_L, int *phys_bdy, real *x, real *y, real *dt,
               real *chr, real *ex, int *Nx, int *Ny, int *regtype);

void init_hb_(real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1,
              real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1,
              real *gb_ty_np1, real *gb_ty_n, real *gb_ty_nm1,
              real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1,
              real *gb_xy_np1, real *gb_xy_n, real *gb_xy_nm1,
              real *gb_yy_np1, real *gb_yy_n, real *gb_yy_nm1,
              real *psi_np1, real *psi_n, real *psi_nm1,
              real *Hb_t_n, real *Hb_x_n, real *Hb_y_n, 
              real *AdS_L, int *phys_bdy, real *x, real *y, real *dt, real *chr, real *ex,
              int *Nx, int *Ny, int *regtype);

void init_ghbdot_(real *gb_tt_n, real *gb_tx_n, real *gb_ty_n,
                  real *gb_xx_n, real *gb_xy_n, real *gb_yy_n, real *psi_n,
                  real *gb_tt_t_n, real *gb_tx_t_n, real *gb_ty_t_n,
                  real *gb_xx_t_n, real *gb_xy_t_n, real *gb_yy_t_n, real *psi_t_n,
                  real *Hb_t_n, real *Hb_x_n, real *Hb_y_n,
                  real *Hb_t_t_n, real *Hb_x_t_n, real *Hb_y_t_n,
                  real *AdS_L, int *phys_bdy, real *x, real *y, real *dt, real *chr, real *ex,
                  int *Nx, int *Ny, int *regtype);

void lin_zero_bnd_res_(real *f, int *phys_bdy, int *all, int *Nx, int *Ny);

void approx_qb_(real *phi_r, real *phi_i, real *phi_r_dot, real *phi_i_dot, 
                real *amp, real *x0, real *delta, real *omega, real *v0, real *AdS_L,
                real *x, real *y, int *Nx, int *Ny);

void dmdiss3d_ex_gen_(real *f,real *work,real *eps,int *do_bdy,int *phys_bdy_type, int *even,
                      int *odd,int *nx,int *ny,
                      int *nz, real *chr, real *ex, int *do_ex, int *ind_sweeps, int *kmax);

#endif
