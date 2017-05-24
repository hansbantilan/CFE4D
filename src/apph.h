#ifndef _APPH_H
#define _APPH_H

//-----------------------------------------------------------------------------
// global variables and prototypes for AH-finder
//-----------------------------------------------------------------------------

extern int AH_Nchi[MAX_BHS],AH_Nphi[MAX_BHS],AH_Lmin[MAX_BHS],AH_Lmax[MAX_BHS];
extern int AH_max_iter[MAX_BHS],AH_freq[MAX_BHS],AH_freq_aft[MAX_BHS],AH_rsteps[MAX_BHS],AH_maxinc[MAX_BHS];
extern real AH_tol[MAX_BHS],AH_tol_aft[MAX_BHS],AH_r0[MAX_BHS],AH_r1[MAX_BHS],AH_lambda[MAX_BHS],AH_lambda_min[MAX_BHS];
extern real AH_eps[MAX_BHS];
extern real AH_tmin[MAX_BHS];
extern real AH_xc[MAX_BHS][2];
extern int c_AH;

extern real *AH_theta[MAX_BHS],*AH_R[MAX_BHS],*AH_w1[MAX_BHS],*AH_w2[MAX_BHS],*AH_w3[MAX_BHS];
extern int *AH_lev[MAX_BHS],*AH_own[MAX_BHS];

extern int AH_count[MAX_BHS],found_AH[MAX_BHS],freq0[MAX_BHS];
extern int found_count_AH[MAX_BHS];

//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------

int find_apph(real *M, real *J, real *c_equat, real *c_polar, int use_R_ic, real *AH_min_resid);

//-----------------------------------------------------------------------------
// related fotran routines
//-----------------------------------------------------------------------------

void is_inside_(int *a_in_b, int *a_int_b, real *r_a, real *xc_a, real *r_b, real *xc_b, int *dim);

void ah_is_int_(int *is_int, real *AH_R, real *AH_xc, int *i, int *j, real *bbox, real *dx, real *dy, 
               int *AH_Nchi, int *AH_Nphi, int *axisym);

void ah_fill_own_(real *AH_R, real *AH_xc, int *AH_own, int *AH_lev, real *bbox, real *dx, real *dy, 
                 int *rank, int *AdS_L, int *AH_Nchi, int *AH_Nphi, int *axisym);

void ah_fill_f_(real *AH_R, real *AH_xc, real *f, int *is, int *ie, int *js, int *je, 
                real *x, real *y, real *z, int *AH_Nchi, int *AH_Nphi, int *Nx, int *Ny, int *Nz, int *axisym);

void calc_exp0_(real *AH_R, real *AH_xc, real *AH_theta, int *i0, int *j0, int *AH_Nchi,int *AH_Nphi, real *theta, real *f, real *da, real *d_ceq, real *d_cp,
                real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1,
                real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1,
                real *gb_ty_np1, real *gb_ty_n, real *gb_ty_nm1,
                real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1,
                real *gb_xy_np1, real *gb_xy_n, real *gb_xy_nm1,
                real *gb_yy_np1, real *gb_yy_n, real *gb_yy_nm1,
                real *psi_np1, real *psi_n, real *psi_nm1,
                real *AdS_L, real *x, real *y, real *z, real *dt, real *chr, 
                real *ex, int *do_ex, int *Nx, int *Ny, int *Nz, int *axisym);

void smooth_ah_r_(real *AH_R,real *AH_w1,real *AH_eps, int *AH_Nchi, int *AH_Nphi);

void smooth_ah_r_b_(real *AH_R,real *AH_w1,real *AH_eps, int *AH_Nchi, int *AH_Nphi);

void reg_ah_r_(real *AH_R, int *AH_Nchi, int *AH_Nphi);

void adjust_ah_xc_(real *AH_R, real *AH_xc, int *AH_Nchi, int *AH_Nphi, real *dx, real *dy, real *dz, int *axisym);

void fill_ex_params_(real *AH_R, real *AH_xc, real *ex_r, real *ex_xc, int *AH_Nchi, int *AH_Nphi, real *dx, real *dy, int *axisym);

#endif
