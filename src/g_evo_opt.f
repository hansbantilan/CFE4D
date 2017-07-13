c----------------------------------------------------------------------
c in cartesian coordinates t,x,y for x in [-1,1], y in [0,1]
c
c An experimental evolution routine for the gb,phi1, computing
c the residual at the time just prior to updated
c
c choosing theta=pi/2
c
c L below is the AdS length scale
c----------------------------------------------------------------------
        subroutine g_evo_opt(eb_res,gb_res,kg_res,cl_res,
     &                       eb_xx_np1,eb_xx_n,eb_xx_nm1,
     &                       eb_xy_np1,eb_xy_n,eb_xy_nm1,
     &                       eb_xz_np1,eb_xz_n,eb_xz_nm1,
     &                       eb_yy_np1,eb_yy_n,eb_yy_nm1,
     &                       eb_yz_np1,eb_yz_n,eb_yz_nm1,
     &                       eb_zz_np1,eb_zz_n,eb_zz_nm1,
     &                       gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                       gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                       gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                       gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                       gb_xy_np1,gb_xy_n,gb_xy_nm1,     
     &                       gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                       psi_np1,psi_n,psi_nm1,
     &                       Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                       Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                       Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                       phi1_np1,phi1_n,phi1_nm1,
     &                       L,x,y,dt,chr,ex,
     &                       phys_bdy,ghost_width,Nx,Ny,
     &                       background,kappa_cd,rho_cd,
     &                       interptype,i_shift,regtype,
     &                       diss_kmax,tfunction)
        implicit none
        integer Nx,Ny
        integer phys_bdy(4),ghost_width(4)
        integer background
        integer interptype
        integer regtype
        integer i_shift
        integer diss_kmax
        integer max_ghost_width
        real*8 kappa_cd,rho_cd
        real*8 eb_res(Nx,Ny),gb_res(Nx,Ny),kg_res(Nx,Ny),cl_res(Nx,Ny)
        real*8 eb_xx_np1(Nx,Ny),eb_xx_n(Nx,Ny),eb_xx_nm1(Nx,Ny)
        real*8 eb_xy_np1(Nx,Ny),eb_xy_n(Nx,Ny),eb_xy_nm1(Nx,Ny)
        real*8 eb_xz_np1(Nx,Ny),eb_xz_n(Nx,Ny),eb_xz_nm1(Nx,Ny)
        real*8 eb_yy_np1(Nx,Ny),eb_yy_n(Nx,Ny),eb_yy_nm1(Nx,Ny)
        real*8 eb_yz_np1(Nx,Ny),eb_yz_n(Nx,Ny),eb_yz_nm1(Nx,Ny)
        real*8 eb_zz_np1(Nx,Ny),eb_zz_n(Nx,Ny),eb_zz_nm1(Nx,Ny)
        real*8 gb_tt_np1(Nx,Ny),gb_tt_n(Nx,Ny),gb_tt_nm1(Nx,Ny)
        real*8 gb_tx_np1(Nx,Ny),gb_tx_n(Nx,Ny),gb_tx_nm1(Nx,Ny)
        real*8 gb_ty_np1(Nx,Ny),gb_ty_n(Nx,Ny),gb_ty_nm1(Nx,Ny)
        real*8 gb_xx_np1(Nx,Ny),gb_xx_n(Nx,Ny),gb_xx_nm1(Nx,Ny)
        real*8 gb_xy_np1(Nx,Ny),gb_xy_n(Nx,Ny),gb_xy_nm1(Nx,Ny)
        real*8 gb_yy_np1(Nx,Ny),gb_yy_n(Nx,Ny),gb_yy_nm1(Nx,Ny)
        real*8 psi_np1(Nx,Ny),psi_n(Nx,Ny),psi_nm1(Nx,Ny)
        real*8 phi1_np1(Nx,Ny),phi1_n(Nx,Ny),phi1_nm1(Nx,Ny)
        real*8 Hb_t_np1(Nx,Ny),Hb_t_n(Nx,Ny),Hb_t_nm1(Nx,Ny)
        real*8 Hb_x_np1(Nx,Ny),Hb_x_n(Nx,Ny),Hb_x_nm1(Nx,Ny)
        real*8 Hb_y_np1(Nx,Ny),Hb_y_n(Nx,Ny),Hb_y_nm1(Nx,Ny)
        real*8 tfunction(Nx,Ny)
        real*8 L
        real*8 x(Nx),y(Ny),dt,chr(Nx,Ny),ex
        real*8 chr2(Nx,Ny)

        integer a,b,c,d,e,f,m,n,p,q
        integer rb,i,j
        integer is,ie,js,je,is_a_nan

        real*8 dx,dy
        real*8 x0,y0,rho0

        real*8 phi1_res,phi1_J

        real*8 cfe(4,4),cfe_J(4,4)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 g0_tt_t, g0_tt_x, g0_tt_y, g0_tt_tt
        real*8 g0_tt_xx,g0_tt_yy,g0_tt_tx,g0_tt_ty
        real*8 g0_tt_xy

        real*8 g0_tx_t, g0_tx_x, g0_tx_y, g0_tx_tt
        real*8 g0_tx_xx,g0_tx_yy,g0_tx_tx,g0_tx_ty
        real*8 g0_tx_xy

        real*8 g0_ty_t, g0_ty_x, g0_ty_y, g0_ty_tt
        real*8 g0_ty_xx,g0_ty_yy,g0_ty_tx,g0_ty_ty
        real*8 g0_ty_xy

        real*8 g0_xx_t, g0_xx_x, g0_xx_y, g0_xx_tt
        real*8 g0_xx_xx,g0_xx_yy,g0_xx_tx,g0_xx_ty
        real*8 g0_xx_xy

        real*8 g0_xy_t, g0_xy_x, g0_xy_y, g0_xy_tt
        real*8 g0_xy_xx,g0_xy_yy,g0_xy_tx,g0_xy_ty
        real*8 g0_xy_xy

        real*8 g0_yy_t, g0_yy_x, g0_yy_y, g0_yy_tt
        real*8 g0_yy_xx,g0_yy_yy,g0_yy_tx,g0_yy_ty
        real*8 g0_yy_xy

        real*8 g0_psi_t, g0_psi_x, g0_psi_y, g0_psi_tt
        real*8 g0_psi_xx,g0_psi_yy,g0_psi_tx,g0_psi_ty
        real*8 g0_psi_xy

        real*8 g0u_tt0,g0_tt0
        real*8 g0u_tx0,g0_tx0
        real*8 g0u_ty0,g0_ty0
        real*8 g0u_xx0,g0_xx0
        real*8 g0u_xy0,g0_xy0
        real*8 g0u_yy0,g0_yy0
        real*8 g0_psi0

        real*8 H0_t0,H0_x0,H0_y0

        real*8 C_t,C_x,C_y
        real*8 C_t_tt_J,C_t_tx_J,C_t_ty_J,C_t_xx_J
        real*8 C_t_xy_J,C_t_yy_J,C_t_psi_J
        real*8 C_x_tt_J,C_x_tx_J,C_x_ty_J,C_x_xx_J
        real*8 C_x_xy_J,C_x_yy_J,C_x_psi_J
        real*8 C_y_tt_J,C_y_tx_J,C_y_ty_J,C_y_xx_J
        real*8 C_y_xy_J,C_y_yy_J,C_y_psi_J
        real*8 nu_t,nu_x,nu_y,nl_t,nl_x,nl_y
        real*8 d_gb_tt_res,d_gb_tx_res,d_gb_ty_res
        real*8 d_gb_xx_res,d_gb_xy_res,d_gb_yy_res
        real*8 d_psi_res
        real*8 d_gb_tt_J,d_gb_tx_J,d_gb_ty_J
        real*8 d_gb_xx_J,d_gb_xy_J,d_gb_yy_J
        real*8 d_psi_J

        logical ltrace,is_nan,dump,first_nan
        logical first_evolved_pt
        parameter (ltrace=.false.)
        data first_nan/.true./

        integer i2,j2

        !--------------------------------------------------------------
        ! using g0_ab = gads_ab + h0_ab
        !       g0^ab = gads^ab + h0^ab
        !       H0_a  = Hads_a + A_a
        ! where h0_ab, h0^ab are not inverses of each other
        !
        ! g0_ll(a,b)           = g0_ab           !use for diagnostics  !
        ! g0_uu(a,b)           = g0^ab           !eg: compare g0_ll    !
        ! g0_ll_x(a,b,c)       = g0_ab_,c        !    with gads_ll+h_ll!  
        ! g0_ll_xx(a,b,c,d)    = g0_ab_,cd
        ! g0_uu_x(a,b,c)       = g0^ab_,c
        !                      = -g0^ad g0^be g0_de_,c 
        !
        ! gads_ll(a,b)         = gads_ab
        ! gads_uu(a,b)         = gads^ab
        ! gads_ll_x(a,b,c)     = gads_ab_,c
        ! gads_ll_xx(a,b,c,d)  = gads_ab_,cd
        ! gads_uu_x(a,b,c)     = gads^ab_,c
        !                      = -gads^ad gads^be gads_de_,c 
        ! 
        ! h0_ll(a,b)           = (1-rho0**2)*gb_ab
        ! h0_uu(a,b)           = inverse(gads+(1-rho0**2)*gb)^ab - gads^ab
        ! h0_ll_x(a,b,c)       = ((1-rho0**2)*gb_ab)_,c
        ! h0_ll_xx(a,b,c,d)    = ((1-rho0**2)*gb_ab)_,cd
        ! h0_uu_x(a,b,c)       = g^ab_,c - gads^ab_,c
        !
        ! gammagg(a,b,c) = 0.5d0*gads_uu(a,d) 
        !                  *(gads_ll_x(c,d,b)-gads_ll_x(b,c,d)+gads_ll_x(d,b,c))
        ! gammahh(a,b,c) = 0.5d0*h0_uu(a,d) 
        !                  *(h0_ll_x(c,d,b)-h0_ll_x(b,c,d)+h0_ll_x(d,b,c))
        ! gammagh(a,b,c) = 0.5d0*gads_uu(a,d) 
        !                  *(h0_ll_x(c,d,b)-h0_ll_x(b,c,d)+h0_ll_x(d,b,c))
        ! gammahg(a,b,c) = 0.5d0*h0_uu(a,d) 
        !                  *(gads_ll_x(c,d,b)-gads_ll_x(b,c,d)+gads_ll_x(d,b,c))
        !
        ! cuuuu(a,b,c,d) = gads_uu(a,b)*gads_uu(c,d)
        !                   +h0_uu(a,b)*h0_uu(c,d) 
        !                   +gads_uu(a,b)*h0_uu(c,d) 
        !                   +h0_uu(a,b)*gads_uu(c,d) 
        !
        ! dlll(a,b,c)    = g0_ll_x(b,c,a)-g0_ll_x(a,b,c)+g0_ll_x(c,a,b)
        !
        ! A_l(a)     = Hb_a
        ! A_l_x(a,b) = Hb_a,b
        ! 
        ! phi10_x(a)= phi1_,a
        !
        ! grad_phi1_sq = g^cd*phi1_,c*phi1_,d
        !
        ! set_ab = 2*phi1_,a*phi1_,b - g_ab*grad_phi1_sq
        ! tr_set= g^cd*se_cd 
        !
        ! efe(a,b) = residual ... hardcoded expressions (see below)
        !
        ! t,x,y=1,2,3
        ! 
        ! NOTE: g0_ll_xx,gads_ll_xx,h0_ll_xx,efe,efe_J
        !       do *NOT* symmetric components filled in
        !
        !--------------------------------------------------------------
        real*8 efe(4,4),efe_J(4,4)
        real*8 term1(4,4),term2(4,4),term3(4,4),term4(4,4)
        real*8 term5(4,4),term6(4,4),term7(4,4),term8(4,4)
        real*8 gammagg(4,4,4),gammahh(4,4,4)
        real*8 gammagh(4,4,4),gammahg(4,4,4) 
        real*8 cuuuu(4,4,4,4),dlll(4,4,4)
        real*8 dphi1(4)

        real*8 cd_ll(4,4),cd_J_ll(4,4)

        real*8 ndotc,n_l(4),n_u(4),c_l(4),c_J_l(4)
        real*8 n_l_x(4,4) 

        real*8 tr_set,grad_phi1_sq
        
        real*8 g0u_tt_ads0,g0u_xx_ads0,g0u_xy_ads0,g0u_yy_ads0
        real*8 g0u_psi_ads0

        real*8 H0_t_ads0,H0_x_ads0,H0_y_ads0

        real*8 dgb_J,ddgb_J,ddgb_J_tx,ddgb_J_ty
        real*8 dc_J

        real*8 lambda4

        real*8 h0_tt0
        real*8 h0_tx0
        real*8 h0_ty0
        real*8 h0_xx0
        real*8 h0_xy0
        real*8 h0_yy0
        real*8 h0_psi0

        real*8 h0_tt_t, h0_tt_x, h0_tt_y, h0_tt_tt
        real*8 h0_tt_xx,h0_tt_yy,h0_tt_tx,h0_tt_ty
        real*8 h0_tt_xy

        real*8 h0_tx_t, h0_tx_x, h0_tx_y, h0_tx_tt
        real*8 h0_tx_xx,h0_tx_yy,h0_tx_tx,h0_tx_ty
        real*8 h0_tx_xy

        real*8 h0_ty_t, h0_ty_x, h0_ty_y, h0_ty_tt
        real*8 h0_ty_xx,h0_ty_yy,h0_ty_tx,h0_ty_ty
        real*8 h0_ty_xy

        real*8 h0_xx_t, h0_xx_x, h0_xx_y, h0_xx_tt
        real*8 h0_xx_xx,h0_xx_yy,h0_xx_tx,h0_xx_ty
        real*8 h0_xx_xy

        real*8 h0_xy_t, h0_xy_x, h0_xy_y, h0_xy_tt
        real*8 h0_xy_xx,h0_xy_yy,h0_xy_tx,h0_xy_ty
        real*8 h0_xy_xy

        real*8 h0_yy_t, h0_yy_x, h0_yy_y, h0_yy_tt
        real*8 h0_yy_xx,h0_yy_yy,h0_yy_tx,h0_yy_ty
        real*8 h0_yy_xy

        real*8 h0_psi_t, h0_psi_x, h0_psi_y, h0_psi_tt
        real*8 h0_psi_xx,h0_psi_yy,h0_psi_tx,h0_psi_ty
        real*8 h0_psi_xy

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,w,y,z)
        !--------------------------------------------------------------
        real*8 e0_ll(4,4),e0_ll_x(4,4,4),e0_ll_xx(4,4,4,4)

        real*8 g0_ll(4,4),g0_uu(4,4)
        real*8 g0_ll_x(4,4,4),g0_uu_x(4,4,4),g0_ll_xx(4,4,4,4)
        real*8 gads_ll(4,4),gads_uu(4,4)
        real*8 gads_ll_x(4,4,4),gads_uu_x(4,4,4),gads_ll_xx(4,4,4,4)
        real*8 h0_ll(4,4),h0_uu(4,4)
        real*8 h0_ll_x(4,4,4),h0_uu_x(4,4,4),h0_ll_xx(4,4,4,4)
        real*8 gamma_ull(4,4,4),gamma_ull_x(4,4,4,4)
        real*8 riemann_ulll(4,4,4,4)
        real*8 ricci_ll(4,4),ricci_lu(4,4),ricci
        real*8 einstein_ll(4,4),set_ll(4,4)
        real*8 Hads_l(4),A_l(4),A_l_x(4,4)
        real*8 phi10_x(4),phi10_xx(4,4)

        !--------------------------------------------------------------
        ! local variables for tensor manipulations 
        !--------------------------------------------------------------
        real*8 weyl(4,4,4,4),weyl_x(4,4,4,4,4)
        !TEMPORARY
        real*8 b0_ll(4,4),b0_ll_x(4,4,4),b0_ll_xx(4,4,4,4) !move to tensor_init()
        real*8 omega !move to tensor_init()

        !--------------------------------------------------------------
        ! initialize fixed-size variables 
        !--------------------------------------------------------------
        data lambda4/0.0/

        data g0u_tt_ads0,g0u_xx_ads0/0.0,0.0/
        data g0u_xy_ads0,g0u_yy_ads0/0.0,0.0/
        data g0u_psi_ads0/0.0/

        data H0_t_ads0,H0_x_ads0,H0_y_ads0/0.0,0.0,0.0/

        data dgb_J,ddgb_J,ddgb_J_tx,ddgb_J_ty/0.0,0.0,0.0,0.0/
        data dc_J/0.0/

        data dlll/64*0.0/
        data cuuuu/256*0.0/

        data term1,term2/16*0.0,16*0.0/
        data term3,term4/16*0.0,16*0.0/
        data term5,term6/16*0.0,16*0.0/
        data term7,term8/16*0.0,16*0.0/

        data cfe,cfe_J/16*0.0,16*0.0/

        data efe,efe_J/16*0.0,16*0.0/
        data cd_ll,cd_J_ll/16*0.0,16*0.0/

        data n_l,n_u,c_l,c_J_l/4*0.0,4*0.0,4*0.0,4*0.0/
        data n_l_x/16*0.0/

        data rb,i,j/0,0,0/
        data i2,j2/0,0/
        data is,ie,js,je,is_a_nan/0,0,0,0,0/
        data a,b,c,d,e/0,0,0,0,0/

        data dx,dy/0.0,0.0/

        data g0_tt_t, g0_tt_x, g0_tt_y, g0_tt_tt/0.0,0.0,0.0,0.0/
        data g0_tt_xx,g0_tt_yy,g0_tt_tx,g0_tt_ty/0.0,0.0,0.0,0.0/
        data g0_tt_xy/0.0/

        data g0_tx_t, g0_tx_x, g0_tx_y, g0_tx_tt/0.0,0.0,0.0,0.0/
        data g0_tx_xx,g0_tx_yy,g0_tx_tx,g0_tx_ty/0.0,0.0,0.0,0.0/
        data g0_tx_xy/0.0/

        data g0_ty_t, g0_ty_x, g0_ty_y, g0_ty_tt/0.0,0.0,0.0,0.0/
        data g0_ty_xx,g0_ty_yy,g0_ty_tx,g0_ty_ty/0.0,0.0,0.0,0.0/
        data g0_ty_xy/0.0/

        data g0_xx_t, g0_xx_x, g0_xx_y, g0_xx_tt/0.0,0.0,0.0,0.0/
        data g0_xx_xx,g0_xx_yy,g0_xx_tx,g0_xx_ty/0.0,0.0,0.0,0.0/
        data g0_xx_xy/0.0/

        data g0_xy_t, g0_xy_x, g0_xy_y, g0_xy_tt/0.0,0.0,0.0,0.0/
        data g0_xy_xx,g0_xy_yy,g0_xy_tx,g0_xy_ty/0.0,0.0,0.0,0.0/
        data g0_xy_xy/0.0/

        data g0_yy_t, g0_yy_x, g0_yy_y, g0_yy_tt/0.0,0.0,0.0,0.0/
        data g0_yy_xx,g0_yy_yy,g0_yy_tx,g0_yy_ty/0.0,0.0,0.0,0.0/
        data g0_yy_xy/0.0/

        data g0_psi_t, g0_psi_x, g0_psi_y, g0_psi_tt/0.0,0.0,0.0,0.0/
        data g0_psi_xx,g0_psi_yy,g0_psi_tx,g0_psi_ty/0.0,0.0,0.0,0.0/
        data g0_psi_xy/0.0/

        data g0u_tt0,g0_tt0/0.0,0.0/
        data g0u_tx0,g0_tx0/0.0,0.0/
        data g0u_ty0,g0_ty0/0.0,0.0/
        data g0u_xx0,g0_xx0/0.0,0.0/
        data g0u_xy0,g0_xy0/0.0,0.0/
        data g0u_yy0,g0_yy0/0.0,0.0/
        data g0_psi0/0.0/

        data H0_t0,H0_x0,H0_y0/0.0,0.0,0.0/

        data x0,y0,rho0/0.0,0.0,0.0/

        data C_t,C_x,C_y/0.0,0.0,0.0/
        data C_t_tt_J,C_t_tx_J,C_t_ty_J,C_t_xx_J/0.0,0.0,0.0,0.0/
        data C_t_xy_J,C_t_yy_J,C_t_psi_J/0.0,0.0,0.0/
        data C_x_tt_J,C_x_tx_J,C_x_ty_J,C_x_xx_J/0.0,0.0,0.0,0.0/
        data C_x_xy_J,C_x_yy_J,C_x_psi_J/0.0,0.0,0.0/
        data C_y_tt_J,C_y_tx_J,C_y_ty_J,C_y_xx_J/0.0,0.0,0.0,0.0/
        data C_y_xy_J,C_y_yy_J,C_y_psi_J/0.0,0.0,0.0/
        data nu_t,nu_x,nu_y,nl_t,nl_x,nl_y/0.0,0.0,0.0,0.0,0.0,0.0/
        data d_gb_tt_res,d_gb_tx_res,d_gb_ty_res/0.0,0.0,0.0/
        data d_gb_xx_res,d_gb_xy_res,d_gb_yy_res/0.0,0.0,0.0/
        data d_psi_res/0.0/
        data d_gb_tt_J,d_gb_tx_J,d_gb_ty_J/0.0,0.0,0.0/
        data d_gb_xx_J,d_gb_xy_J,d_gb_yy_J/0.0,0.0,0.0/
        data d_psi_J/0.0/

        data grad_phi1_sq/1*0.0/
        data Hads_l,A_l,dphi1/4*0.0,4*0.0,4*0.0/
        data A_l_x/16*0.0/

        data e0_ll,e0_ll_x,e0_ll_xx/16*0.0,64*0.0,256*0.0/
        data b0_ll,b0_ll_x,b0_ll_xx/16*0.0,64*0.0,256*0.0/

        data weyl,weyl_x/256*0.0,1024*0/
        data omega/0.0/

        data g0_ll,g0_uu,gads_ll/16*0.0,16*0.0,16*0.0/
        data gads_uu,h0_ll,h0_uu/16*0.0,16*0.0,16*0.0/
        data gammagg,gammahh/64*0.0,64*0.0/
        data gammagh,gammahg/64*0.0,64*0.0/
        data g0_ll_x,g0_uu_x/64*0.0,64*0.0/
        data gads_ll_x,gads_uu_x/64*0.0,64*0.0/
        data h0_ll_x,h0_uu_x/64*0.0,64*0.0/
        data g0_ll_xx/256*0.0/
        data gads_ll_xx,h0_ll_xx/256*0.0,256*0.0/

        data gamma_ull/64*0.0/
        data gamma_ull_x/256*0.0/
        data riemann_ulll/256*0.0/
        data ricci/0.0/
        data ricci_ll,ricci_lu/16*0.0,16*0.0/
        data einstein_ll,set_ll/16*0.0,16*0.0/

        data phi10_x/4*0.0/
        data phi10_xx/16*0.0/

        !--------------------------------------------------------------
        if (ltrace) write(*,*) 'gb_psi_evo ... N=',Nx,Ny

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        if (abs((y(2)-y(1))/dx-1).gt.1.0d-8) then
           write(*,*) 'error ... g_evo_opt not updated for dx!=dy!=dz'
           stop
        end if

        ! CFE4D cosmological constant
        !(lambda4=-(n-1)(n-2)/2/L^2) for n=4 dimensional AdS)
        lambda4=-3/L/L

        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1

        !(nearest-to-axis points are not evolved, according to regtype choice) 
        if (regtype.eq.7 .or. regtype.eq.6) then
          if (abs(y(1)).lt.dy/2) js=4
        else if (regtype.eq.5 .or. regtype.eq.4 .or. regtype.eq.3) then
          if (abs(y(1)).lt.dy/2) js=3
        else
          if (abs(y(1)).lt.dy/2) js=2
        endif

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)

        ! check kmax value against ghost_width
        max_ghost_width=max(ghost_width(1),ghost_width(2),
     &                      ghost_width(3),ghost_width(4))
        if (max_ghost_width.lt.2*diss_kmax) then
          write(*,*) 'WARNING ... ghost_width < 2*diss_kmax'
          write(*,*) 'max{ghost_width}=',max_ghost_width
          write(*,*) 'diss_kmax=',diss_kmax
          stop
        endif

        ! zero all outer boundary points
        if (phys_bdy(1).eq.1) then 
          do j=1,Ny
            gb_tt_np1(1,j) = 0
            gb_tx_np1(1,j) = 0
            gb_ty_np1(1,j) = 0
            gb_xx_np1(1,j) = 0
            gb_xy_np1(1,j) = 0
            gb_yy_np1(1,j) = 0
            psi_np1(1,j) = 0
            phi1_np1(1,j) = 0
          end do
        end if
        if (phys_bdy(2).eq.1) then 
          do j=1,Ny
            gb_tt_np1(Nx,j) = 0
            gb_tx_np1(Nx,j) = 0
            gb_ty_np1(Nx,j) = 0
            gb_xx_np1(Nx,j) = 0
            gb_xy_np1(Nx,j) = 0
            gb_yy_np1(Nx,j) = 0
            psi_np1(Nx,j) = 0
            phi1_np1(Nx,j) = 0
          end do
        end if
        if (phys_bdy(4).eq.1) then 
          do i=1,Nx
            gb_tt_np1(i,Ny) = 0
            gb_tx_np1(i,Ny) = 0
            gb_ty_np1(i,Ny) = 0
            gb_xx_np1(i,Ny) = 0
            gb_xy_np1(i,Ny) = 0
            gb_yy_np1(i,Ny) = 0
            psi_np1(i,Ny) = 0
            phi1_np1(i,Ny) = 0
          end do
        end if

        ! define chr2
        do i=is,ie
          do j=js,je
            if (chr(i,j).ne.ex.and.
     &          sqrt(x(i)**2+y(j)**2).ge.1.0d0-3*dx/2.and.
     &          (chr(i-1,j).eq.ex.or.chr(i+1,j).eq.ex.or.
     &           chr(i,j-1).eq.ex.or.chr(i,j+1).eq.ex)) then
              chr2(i,j)=ex
            else 
              chr2(i,j)=ex-1
            end if
          end do
        end do

        !(MAIN LOOP) red-black loop through spacetime points x(i),y(j)  
        do rb=1,2

          do j=js,je
            do i=is+mod(j+rb,2),ie,2
              x0=x(i)
              y0=y(j)
              rho0=sqrt(x0**2+y0**2)

              dump=.false.

              if (ltrace) write(*,*) 'i,j:',i,j

              first_evolved_pt=.false.

              ! define first_evolved_pt
              if (chr(i,j).ne.ex.and.chr2(i,j).ne.ex.and.
     &            (chr2(i-1,j).eq.ex.or.chr2(i+1,j).eq.ex.or.
     &             chr2(i,j-1).eq.ex.or.chr2(i,j+1).eq.ex)) then
                first_evolved_pt=.true.
              end if

              !(REGION) interior not one-point-away-from-ads-bdy points; evolve 
              if (chr(i,j).ne.ex .and. chr2(i,j).ne.ex) then

              ! computes tensors at point i,j
              call tensor_init(
     &                eb_xx_np1,eb_xx_n,eb_xx_nm1,
     &                eb_xy_np1,eb_xy_n,eb_xy_nm1,
     &                eb_xz_np1,eb_xz_n,eb_xz_nm1,
     &                eb_yy_np1,eb_yy_n,eb_yy_nm1,
     &                eb_yz_np1,eb_yz_n,eb_yz_nm1,
     &                eb_zz_np1,eb_zz_n,eb_zz_nm1,
     &                gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                psi_np1,psi_n,psi_nm1,
     &                Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                phi1_np1,phi1_n,phi1_nm1,
     &                e0_ll,e0_ll_x,e0_ll_xx,
     &                g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &                gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &                h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                A_l,A_l_x,Hads_l,
     &                gamma_ull,gamma_ull_x,
     &                riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                einstein_ll,set_ll,
     &                phi10_x,phi10_xx,
     &                x,y,dt,chr,L,ex,Nx,Ny,i,j)

                do a=1,4
                  do b=1,4
                    do c=1,4
                      dlll(a,b,c)=
     &                    g0_ll_x(b,c,a)-g0_ll_x(a,b,c)+g0_ll_x(c,a,b)
                      gammagg(a,b,c)=0
                      gammahh(a,b,c)=0
                      gammagh(a,b,c)=0
                      gammahg(a,b,c)=0
                      do d=1,4
                        gammagg(a,b,c)=gammagg(a,b,c)
     &                                 +0.5d0*gads_uu(a,d)
     &                                      *(gads_ll_x(c,d,b)
     &                                       -gads_ll_x(b,c,d)
     &                                       +gads_ll_x(d,b,c))
                        gammahh(a,b,c)=gammahh(a,b,c)
     &                                 +0.5d0*h0_uu(a,d)
     &                                      *(h0_ll_x(c,d,b)
     &                                       -h0_ll_x(b,c,d)
     &                                       +h0_ll_x(d,b,c))
                        gammagh(a,b,c)=gammagh(a,b,c)
     &                                 +0.5d0*gads_uu(a,d)
     &                                      *(h0_ll_x(c,d,b)
     &                                       -h0_ll_x(b,c,d)
     &                                       +h0_ll_x(d,b,c))
                        gammahg(a,b,c)=gammahg(a,b,c)
     &                                 +0.5d0*h0_uu(a,d)
     &                                      *(gads_ll_x(c,d,b)
     &                                       -gads_ll_x(b,c,d)
     &                                       +gads_ll_x(d,b,c))
                        cuuuu(a,b,c,d)=gads_uu(a,b)*gads_uu(c,d)+
     &                                 h0_uu(a,b)*h0_uu(c,d)+
     &                                 gads_uu(a,b)*h0_uu(c,d)+
     &                                 h0_uu(a,b)*gads_uu(c,d)
                      end do
                    end do
                  end do
                end do

                !---------------------------------------------------------------- 
                ! Analytically remove the pure AdS terms in the C_l,
                ! denoted ->, so 
                ! C_a = (Hads_a+A_a) - boxx_a
                ! for
                ! 
                ! Hads_a - boxx_a =  (1/2 gads^cd gads_cd,a - gads^cd
                ! gads_ca,d )
                !                   -(1/2 g^cd g_cd,a - g^cd g_ca,d )
                !
                !                 ->-(1/2 h^cd h_cd,a - h^cd h_ca,d )
                !                   -(1/2 gads^cd h_cd,a - gads^cd
                !                   h_ca,d )
                !                   -(1/2 h^cd gads_cd,a - h^cd
                !                   gads_ca,d )
                !----------------------------------------------------------------
                do a=1,4
                  c_l(a)=A_l(a)
     &                   -( 0.5d0*( h0_uu(1,1)*h0_ll_x(1,1,a)+
     &                              h0_uu(2,2)*h0_ll_x(2,2,a)+
     &                              h0_uu(3,3)*h0_ll_x(3,3,a)+
     &                              h0_uu(4,4)*h0_ll_x(4,4,a)+
     &                           2*(h0_uu(1,2)*h0_ll_x(1,2,a)+
     &                              h0_uu(1,3)*h0_ll_x(1,3,a)+
     &                              h0_uu(1,4)*h0_ll_x(1,4,a)+
     &                              h0_uu(2,3)*h0_ll_x(2,3,a)+
     &                              h0_uu(2,4)*h0_ll_x(2,4,a)+
     &                              h0_uu(3,4)*h0_ll_x(3,4,a)) )
     &                     -1.0d0*( h0_uu(1,1)*h0_ll_x(1,a,1)+
     &                              h0_uu(2,2)*h0_ll_x(2,a,2)+
     &                              h0_uu(3,3)*h0_ll_x(3,a,3)+
     &                              h0_uu(4,4)*h0_ll_x(4,a,4)+
     &                             (h0_uu(1,2)*h0_ll_x(1,a,2)+
     &                              h0_uu(2,1)*h0_ll_x(2,a,1)+
     &                              h0_uu(1,3)*h0_ll_x(1,a,3)+
     &                              h0_uu(3,1)*h0_ll_x(3,a,1)+
     &                              h0_uu(1,4)*h0_ll_x(1,a,4)+
     &                              h0_uu(4,1)*h0_ll_x(4,a,1)+
     &                              h0_uu(2,3)*h0_ll_x(2,a,3)+
     &                              h0_uu(3,2)*h0_ll_x(3,a,2)+
     &                              h0_uu(2,4)*h0_ll_x(2,a,4)+
     &                              h0_uu(4,2)*h0_ll_x(4,a,2)+
     &                              h0_uu(3,4)*h0_ll_x(3,a,4)+
     &                              h0_uu(4,3)*h0_ll_x(4,a,3)) ) )
     &                   -( 0.5d0*( gads_uu(1,1)*h0_ll_x(1,1,a)+
     &                              gads_uu(2,2)*h0_ll_x(2,2,a)+
     &                              gads_uu(3,3)*h0_ll_x(3,3,a)+
     &                              gads_uu(4,4)*h0_ll_x(4,4,a)+
     &                           2*(gads_uu(1,2)*h0_ll_x(1,2,a)+
     &                              gads_uu(1,3)*h0_ll_x(1,3,a)+
     &                              gads_uu(1,4)*h0_ll_x(1,4,a)+
     &                              gads_uu(2,3)*h0_ll_x(2,3,a)+
     &                              gads_uu(2,4)*h0_ll_x(2,4,a)+
     &                              gads_uu(3,4)*h0_ll_x(3,4,a)) )
     &                     -1.0d0*( gads_uu(1,1)*h0_ll_x(1,a,1)+
     &                              gads_uu(2,2)*h0_ll_x(2,a,2)+
     &                              gads_uu(3,3)*h0_ll_x(3,a,3)+
     &                              gads_uu(4,4)*h0_ll_x(4,a,4)+
     &                             (gads_uu(1,2)*h0_ll_x(1,a,2)+
     &                              gads_uu(2,1)*h0_ll_x(2,a,1)+
     &                              gads_uu(1,3)*h0_ll_x(1,a,3)+
     &                              gads_uu(3,1)*h0_ll_x(3,a,1)+
     &                              gads_uu(1,4)*h0_ll_x(1,a,4)+
     &                              gads_uu(4,1)*h0_ll_x(4,a,1)+
     &                              gads_uu(2,3)*h0_ll_x(2,a,3)+
     &                              gads_uu(3,2)*h0_ll_x(3,a,2)+
     &                              gads_uu(2,4)*h0_ll_x(2,a,4)+
     &                              gads_uu(4,2)*h0_ll_x(4,a,2)+
     &                              gads_uu(3,4)*h0_ll_x(3,a,4)+
     &                              gads_uu(4,3)*h0_ll_x(4,a,3)) ) )
     &                   -( 0.5d0*( h0_uu(1,1)*gads_ll_x(1,1,a)+
     &                              h0_uu(2,2)*gads_ll_x(2,2,a)+
     &                              h0_uu(3,3)*gads_ll_x(3,3,a)+
     &                              h0_uu(4,4)*gads_ll_x(4,4,a)+
     &                           2*(h0_uu(1,2)*gads_ll_x(1,2,a)+
     &                              h0_uu(1,3)*gads_ll_x(1,3,a)+
     &                              h0_uu(1,4)*gads_ll_x(1,4,a)+
     &                              h0_uu(2,3)*gads_ll_x(2,3,a)+
     &                              h0_uu(2,4)*gads_ll_x(2,4,a)+
     &                              h0_uu(3,4)*gads_ll_x(3,4,a)) )
     &                     -1.0d0*( h0_uu(1,1)*gads_ll_x(1,a,1)+
     &                              h0_uu(2,2)*gads_ll_x(2,a,2)+
     &                              h0_uu(3,3)*gads_ll_x(3,a,3)+
     &                              h0_uu(4,4)*gads_ll_x(4,a,4)+
     &                             (h0_uu(1,2)*gads_ll_x(1,a,2)+
     &                              h0_uu(2,1)*gads_ll_x(2,a,1)+
     &                              h0_uu(1,3)*gads_ll_x(1,a,3)+
     &                              h0_uu(3,1)*gads_ll_x(3,a,1)+
     &                              h0_uu(1,4)*gads_ll_x(1,a,4)+
     &                              h0_uu(4,1)*gads_ll_x(4,a,1)+
     &                              h0_uu(2,3)*gads_ll_x(2,a,3)+
     &                              h0_uu(3,2)*gads_ll_x(3,a,2)+
     &                              h0_uu(2,4)*gads_ll_x(2,a,4)+
     &                              h0_uu(4,2)*gads_ll_x(4,a,2)+
     &                              h0_uu(3,4)*gads_ll_x(3,a,4)+
     &                              h0_uu(4,3)*gads_ll_x(4,a,3)) ) )
                end do

                n_l(1)=-1/sqrt(-g0_uu(1,1))
                do a=1,4
                  n_u(a)=n_l(1)*g0_uu(a,1)+
     &                   n_l(2)*g0_uu(a,2)+
     &                   n_l(3)*g0_uu(a,3)+
     &                   n_l(4)*g0_uu(a,4)
                end do
                do b=1,4
                  n_l_x(1,b)=-1/2.0d0/sqrt(-g0_uu(1,1))**3
     &                        *g0_uu_x(1,1,b)
                end do

                ndotc  =n_u(1)*c_l(1)+
     &                  n_u(2)*c_l(2)+
     &                  n_u(3)*c_l(3)+
     &                  n_u(4)*c_l(4)
 
                grad_phi1_sq=phi10_x(1)*phi10_x(1)*g0_uu(1,1)+
     &                       phi10_x(2)*phi10_x(2)*g0_uu(2,2)+
     &                       phi10_x(3)*phi10_x(3)*g0_uu(3,3)+
     &                       phi10_x(4)*phi10_x(4)*g0_uu(4,4)+
     &                    2*(phi10_x(1)*phi10_x(2)*g0_uu(1,2)+
     &                       phi10_x(1)*phi10_x(3)*g0_uu(1,3)+
     &                       phi10_x(1)*phi10_x(4)*g0_uu(1,4)+
     &                       phi10_x(2)*phi10_x(3)*g0_uu(2,3)+
     &                       phi10_x(2)*phi10_x(4)*g0_uu(2,4)+
     &                       phi10_x(3)*phi10_x(4)*g0_uu(3,4))

                do a=1,4
                  do b=1,4
                    set_ll(a,b)=phi10_x(a)*phi10_x(b)
     &                         -g0_ll(a,b)*(grad_phi1_sq/2)
                  end do
                end do

                tr_set =set_ll(1,1)*g0_uu(1,1)+
     &                  set_ll(2,2)*g0_uu(2,2)+
     &                  set_ll(3,3)*g0_uu(3,3)+
     &                  set_ll(4,4)*g0_uu(4,4)+
     &               2*(set_ll(1,2)*g0_uu(1,2)+
     &                  set_ll(1,3)*g0_uu(1,3)+
     &                  set_ll(1,4)*g0_uu(1,4)+
     &                  set_ll(2,3)*g0_uu(2,3)+
     &                  set_ll(2,4)*g0_uu(2,4)+
     &                  set_ll(3,4)*g0_uu(3,4))

                do a=1,4
                  do b=1,4
                    do c=1,4
                      do d=1,4
                        weyl(a,b,c,d)=
     &                   2*e0_ll(b,d)*n_l(a)*n_l(c)
     &                  -2*e0_ll(a,d)*n_l(b)*n_l(c)
     &                  -2*e0_ll(b,c)*n_l(a)*n_l(d)
     &                  +2*e0_ll(a,c)*n_l(b)*n_l(d)
     &                  +e0_ll(b,d)*g0_ll(a,c)
     &                  -e0_ll(b,c)*g0_ll(a,d)
     &                  -e0_ll(a,d)*g0_ll(b,c)
     &                  +e0_ll(a,c)*g0_ll(b,d)
                        do e=1,4
                          weyl_x(a,b,c,d,e)=
     &                     2*e0_ll_x(b,d,e)*n_l(a)*n_l(c)
     &                    +2*e0_ll(b,d)*n_l_x(a,e)*n_l(c)
     &                    +2*e0_ll(b,d)*n_l(a)*n_l_x(c,e)
     &                    -2*e0_ll_x(a,d,e)*n_l(b)*n_l(c)
     &                    -2*e0_ll(a,d)*n_l_x(b,e)*n_l(c)
     &                    -2*e0_ll(a,d)*n_l(b)*n_l_x(c,e)
     &                    -2*e0_ll_x(b,c,e)*n_l(a)*n_l(d)
     &                    -2*e0_ll(b,c)*n_l_x(a,e)*n_l(d)
     &                    -2*e0_ll(b,c)*n_l(a)*n_l_x(d,e)
     &                    +2*e0_ll_x(a,c,e)*n_l(b)*n_l(d)
     &                    +2*e0_ll(a,c)*n_l_x(b,e)*n_l(d)
     &                    +2*e0_ll(a,c)*n_l(b)*n_l_x(d,e)
     &                    +e0_ll_x(b,d,e)*g0_ll(a,c)
     &                    +e0_ll(b,d)*g0_ll_x(a,c,e)
     &                    -e0_ll_x(b,c,e)*g0_ll(a,d)
     &                    -e0_ll(b,c)*g0_ll_x(a,d,e)
     &                    -e0_ll_x(a,d,e)*g0_ll(b,c)
     &                    -e0_ll(a,d)*g0_ll_x(b,c,e)
     &                    +e0_ll_x(a,c,e)*g0_ll(b,d)
     &                    +e0_ll(a,c)*g0_ll_x(b,d,e)
                        !TEMPoRARY: eventually add magnetic terms with levi-civita symbols
                        !  do f=1,4
                        !
                        !  end do
                        end do        
                      end do
                    end do
                  end do
                end do

                !TEMPORARY: move this to tensor_init()
                ! check that ricci = -12/L^2
                ! conformal factor from between (14d) and (15a) in ConformalWaveEqnsSummary.pdf
                omega=(1-rho0**2)/2
 
                !--------------------------------------------------------------------------
                ! cfe = 
                !--------------------------------------------------------------------------
                do a=1,4
                  do b=1,4

                    cfe(a,b)=0

                    do c=1,4
                      do d=1,4

                        cfe(a,b)=cfe(a,b)
     &                  +g0_uu(d,c)*e0_ll_xx(a,b,c,d)

                        do e=1,4
                          cfe(a,b)=cfe(a,b)
     &                    -gamma_ull(d,c,e)*g0_uu(e,c)*e0_ll_x(a,b,d)
     &                    -gamma_ull(e,c,b)*g0_uu(d,c)*e0_ll_x(a,e,d)
     &                    -gamma_ull(e,c,a)*g0_uu(d,c)*e0_ll_x(e,b,d)
     &                    -e0_ll(e,b)*g0_uu(d,c)*gamma_ull_x(e,d,a,c)
     &                    -e0_ll(a,e)*g0_uu(d,c)*gamma_ull_x(e,d,b,c)
     &                    -gamma_ull(e,d,b)*g0_uu(d,c)*e0_ll_x(a,e,c)
     &                    -gamma_ull(e,d,a)*g0_uu(d,c)*e0_ll_x(e,b,c) !run10

                          do f=1,4

                            cfe(a,b)=cfe(a,b)
     &                      +gamma_ull(e,d,f)*gamma_ull(f,c,b)
     &                                  *e0_ll(a,e)*g0_uu(d,c) !run4
     &                      +gamma_ull(e,d,f)*gamma_ull(f,c,a)
     &                                  *e0_ll(e,b)*g0_uu(d,c)
     &                      +gamma_ull(e,c,a)*gamma_ull(f,d,b)
     &                                  *e0_ll(e,f)*g0_uu(d,c)
     &                      +gamma_ull(e,d,a)*gamma_ull(f,c,b)
     &                                  *e0_ll(e,f)*g0_uu(d,c)
     &                      +gamma_ull(d,e,f)*gamma_ull(c,d,b)
     &                                  *e0_ll(a,c)*g0_uu(f,e)
     &                      +gamma_ull(d,e,f)*gamma_ull(c,d,a)
     &                                  *e0_ll(c,b)*g0_uu(f,e) !run5
     &                      -(n_l(d)*n_l(e)*weyl(a,f,b,c)*ricci
     &                        *g0_uu(f,d)*g0_uu(c,e)           )/2 !run9

                            do m=1,4
                              do n=1,4

                                cfe(a,b)=cfe(a,b)
     &                          -2*n_l(c)*g0_uu(d,c)*g0_uu(e,f)
     *                           *g0_uu(m,n)*weyl(a,d,b,e)

!                                do p=1,4
!                                  do q=1,4
!                                    cfe(a,b)=cfe(a,b)
!     &                              -gamma_ull(c,d,e)*gamma_ull(e,f,m)
!     &                               *n_l(c)*n_l(n)*weyl(a,p,b,q)
!     &                               *g0_uu(p,n)*g0_uu(q,m)*g0_uu(f,d)
!                                  end do
!                                end do  

                              end do
                            end do  

                          end do
                        end do

                      end do
                    end do

                  end do
                end do

                !--------------------------------------------------------------------------
                ! phi1_res = phi1,ab g^ab 
                !--------------------------------------------------------------------------
                phi1_res= phi10_xx(1,1)*g0_uu(1,1)+
     &                    phi10_xx(2,2)*g0_uu(2,2)+
     &                    phi10_xx(3,3)*g0_uu(3,3)+
     &                    phi10_xx(4,4)*g0_uu(4,4)+
     &                 2*(phi10_xx(1,2)*g0_uu(1,2)+
     &                    phi10_xx(1,3)*g0_uu(1,3)+
     &                    phi10_xx(1,4)*g0_uu(1,4)+
     &                    phi10_xx(2,3)*g0_uu(2,3)+
     &                    phi10_xx(2,4)*g0_uu(2,4)+
     &                    phi10_xx(3,4)*g0_uu(3,4))

                !---------------------------------------------------------------- 
                ! computes diag. Jacobian of g_np1->L.g_np1 transformation
                ! by differentiating L.g wrt. g(a,b)_ij_np1 diag. entries
                ! 
                ! ddgb_J_tx,ddgb_J_ty differ from ddgb_J due to forward/backward stencils
                ! at excision surfaces that affect the cross-derivatives tx,ty 
                ! (these are the only contributions, since the diag Jacobian is diff wrt. g_ij_np1)
                !---------------------------------------------------------------- 
                dgb_J=1/2/dt
                ddgb_J=1/dt/dt

                if (i.eq.1.or.(chr(i-1,j).eq.ex)) then
                   if (i.le.(Nx-3)
     &                 .and.((chr(i+1,j).ne.ex
     &                 .and.chr(i+2,j).ne.ex
     &                 .and.chr(i+3,j).ne.ex))) then
                      ddgb_J_tx=-1/dt/dx
                   else if (i.le.(Nx-2)
     &                      .and.((chr(i+1,j).ne.ex
     &                      .and.chr(i+2,j).ne.ex))) then
                      ddgb_J_tx=-3/4/dt/dx
                   else if (i.le.(Nx-1).and.chr(i+1,j).ne.ex) then
                      ddgb_J_tx=-1/2/dt/dx
                   else
                      write(*,*) 'g_evo_opt: error in chr stencil (A)'
                      write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
                      write(*,*) '    (first error only)'
                      ddgb_J_tx=0
                   end if
                else if (i.eq.Nx.or.(chr(i+1,j).eq.ex)) then
                   if (i.ge.4
     &                 .and.((chr(i-1,j).ne.ex
     &                 .and.chr(i-2,j).ne.ex
     &                 .and.chr(i-3,j).ne.ex))) then
                      ddgb_J_tx=1/dt/dx
                   else if (i.ge.3
     &                      .and.((chr(i-1,j).ne.ex
     &                      .and.chr(i-2,j).ne.ex))) then
                      ddgb_J_tx=3/4/dt/dx
                   else if (i.ge.2.and.chr(i-1,j).ne.ex) then
                      ddgb_J_tx=1/2/dt/dx
                   else
                      write(*,*) 'g_evo_opt: error in chr stencil (B)'
                      write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
                      write(*,*) '    (first error only)'
                      ddgb_J_tx=0
                   end if
                else
                   if ((chr(i+1,j).ne.ex.and.chr(i-1,j).ne.ex)) then
                      ddgb_J_tx=0
                   else
                      write(*,*) 'g_evo_opt: error in chr stencil (C)'
                      write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
                      write(*,*) '    (first error only)'
                      ddgb_J_tx=0
                   end if
                end if

                if ((j.eq.1).or.(chr(i,j-1).eq.ex)) then
                   if (j.le.(Ny-3)
     &                 .and.((chr(i,j+1).ne.ex
     &                 .and.chr(i,j+2).ne.ex
     &                 .and.chr(i,j+3).ne.ex))) then
                      ddgb_J_ty=-1/dt/dy
                   else if (j.le.(Ny-2).and.((chr(i,j+1).ne.ex
     &                      .and.chr(i,j+2).ne.ex))) then
                      ddgb_J_ty=-3/4/dt/dy              
                   else if (j.le.(Ny-1).and.chr(i,j+1).ne.ex) then
                      ddgb_J_ty=-1/2/dt/dy
                   else
                      write(*,*) 'g_evo_opt: error in chr stencil (D)'
                      write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
                      write(*,*) '    (first error only)'
                      ddgb_J_ty=0
                   end if
                else if ((j.eq.Ny).or.(chr(i,j+1).eq.ex)) then
                   if (j.ge.4
     &                 .and.((chr(i,j-1).ne.ex
     &                 .and.chr(i,j-2).ne.ex
     &                 .and.chr(i,j-3).ne.ex))) then
                      ddgb_J_ty=1/dt/dy
                   else if (j.ge.3.and.((chr(i,j-1).ne.ex
     &                      .and.chr(i,j-2).ne.ex))) then
                      ddgb_J_ty=3/4/dt/dy
                   else if (j.ge.2.and.chr(i,j-1).ne.ex) then
                      ddgb_J_ty=1/2/dt/dy
                   else
                      write(*,*) 'g_evo_opt: error in chr stencil (E)'
                      write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
                      write(*,*) '    (first error only)'
                      ddgb_J_ty=0
                   end if
                else
                   if ((chr(i,j+1).ne.ex.and.chr(i,j-1).ne.ex)) then
                      ddgb_J_ty=0
                   else
                      write(*,*) 'g_evo_opt: error in chr stencil (F)'
                      write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
                      write(*,*) '    (first error only)'
                      ddgb_J_ty=0
                   end if
                end if

                !----------------------------------------------------------------
                ! cfe_J =
                ! (only include terms in cfe that contain factors of 
                !  e0_ll_xx(a,b,1,1) and replace those factors by ddgb_J)
                !----------------------------------------------------------------
                do a=1,4
                  do b=1,4

                    cfe_J(a,b)=
     &              +g0_uu(1,1)*ddgb_J
                        
                  end do
                end do

                !----------------------------------------------------------------
                ! computes diag. Jacobian of phi1_np1->L.phi1_np1 transformation
                ! by differentiating L.phi1=box.phi1-dV/dphi1 wrt. phi1_np1
                ! and remember: phi10=phi1*(1-rho0**2)**3 
                ! 
                ! (re-use the dgb_J,ddgb_J,ddgb_J_tx,ddgb_J_ty defined above)
                !----------------------------------------------------------------
                phi1_J=    (
     &                 g0_uu(1,1)*ddgb_J
     &                 +2*x0*g0_uu(1,2)*dgb_J
     &                 +2*g0_uu(1,2)*ddgb_J_tx
     &                 +2*y0*g0_uu(1,3)*dgb_J
     &                 +2*g0_uu(1,3)*ddgb_J_ty
     &                     )

                ! constraint damping terms added to efe,efe_J
                do a=1,4
                  do b=1,4
                    cd_ll(a,b)=-kappa_cd*
     &                 (
     &                  n_l(a)*c_l(b)+n_l(b)*c_l(a)
     &                 -(1+rho_cd)*g0_ll(a,b)*ndotc
     &                 )
                  end do
                end do

                dc_J=1d0/2d0/dt

                cd_J_ll(1,1)=-kappa_cd*
     &              (
     &               n_l(1)*g0_uu(1,1)*dc_J*(1-rho0**2)
     &              -(1+rho_cd)*g0_ll(1,1)*n_u(1)*
     &                    (0.5d0*g0_uu(1,1)*dc_J*(1-rho0**2))
     &              )
                cd_J_ll(1,2)=-kappa_cd*
     &              (
     &               n_l(1)*g0_uu(1,1)*dc_J*(1-rho0**2)
     &              -(1+rho_cd)*g0_ll(1,2)*n_u(2)*
     &                    (g0_uu(1,1)*dc_J*(1-rho0**2))
     &              )
                cd_J_ll(1,3)=-kappa_cd*
     &              (
     &               n_l(1)*g0_uu(1,1)*dc_J*(1-rho0**2)
     &              -(1+rho_cd)*g0_ll(1,3)*n_u(3)*
     &                    (g0_uu(1,1)*dc_J*(1-rho0**2))
     &              )
                cd_J_ll(2,2)=-kappa_cd*
     &              (
     &               2*n_l(2)*g0_uu(1,2)*dc_J*(1-rho0**2)
     &              -(1+rho_cd)*g0_ll(2,2)*
     &                    (-n_u(1)*(0.5d0*g0_uu(2,2)*dc_J*(1-rho0**2))
     &                     +n_u(2)*(g0_uu(1,2)*dc_J*(1-rho0**2)))
     &              )
                cd_J_ll(2,3)=-kappa_cd*
     &              (
     &               n_l(2)*g0_uu(1,2)*dc_J*(1-rho0**2)
     &              +n_l(3)*g0_uu(1,3)*dc_J*(1-rho0**2)
     &              -(1+rho_cd)*g0_ll(2,2)*
     &                    (-n_u(1)*(0.5d0*g0_uu(2,3)*dc_J*(1-rho0**2))
     &                     +n_u(2)*(g0_uu(1,3)*dc_J*(1-rho0**2))
     &                     +n_u(3)*(g0_uu(1,2)*dc_J*(1-rho0**2)))
     &              )
                cd_J_ll(3,3)=-kappa_cd*
     &              (
     &               2*n_l(3)*g0_uu(1,3)*dc_J*(1-rho0**2)
     &              -(1+rho_cd)*g0_ll(3,3)*
     &                    (-n_u(1)*(0.5d0*g0_uu(3,3)*dc_J*(1-rho0**2))
     &                     +n_u(3)*(g0_uu(1,3)*dc_J*(1-rho0**2)))
     &              )
                cd_J_ll(4,4)=-kappa_cd*
     &              (
     &              -(1+rho_cd)*g0_ll(4,4)*n_u(1)*
     &                    (0.5d0*g0_uu(4,4)*dc_J*(1-rho0**2)*y0**2)
     &              )

                if (kappa_cd.ne.0) then
                  efe(1,1)=efe(1,1)+cd_ll(1,1)
                  efe(1,2)=efe(1,2)+cd_ll(1,2)
                  efe(1,3)=efe(1,3)+cd_ll(1,3)
                  efe(2,2)=efe(2,2)+cd_ll(2,2)
                  efe(2,3)=efe(2,3)+cd_ll(2,3)
                  efe(3,3)=efe(3,3)+cd_ll(3,3)
                  efe(4,4)=efe(4,4)+cd_ll(4,4)
                  efe_J(1,1)=efe_J(1,1)+cd_J_ll(1,1)
                  efe_J(1,2)=efe_J(1,2)+cd_J_ll(1,2)
                  efe_J(1,3)=efe_J(1,3)+cd_J_ll(1,3)
                  efe_J(2,2)=efe_J(2,2)+cd_J_ll(2,2)
                  efe_J(2,3)=efe_J(2,3)+cd_J_ll(2,3)
                  efe_J(3,3)=efe_J(3,3)+cd_J_ll(3,3)
                  efe_J(4,4)=efe_J(4,4)+cd_J_ll(4,4)
                end if

                ! updates ebars
                if (is_nan(cfe(2,2)).or.is_nan(cfe_J(2,2)).or.
     &            cfe_J(2,2).eq.0) then
                  dump=.true.
                else
                  eb_xx_np1(i,j)=eb_xx_np1(i,j)-cfe(2,2)/cfe_J(2,2)
                end if

                ! update phi1
                if (is_nan(phi1_res).or.is_nan(phi1_J).or.
     &            phi1_J.eq.0) then
                  dump=.true.
                else
                  phi1_np1(i,j)=phi1_np1(i,j)-phi1_res/phi1_J
                end if

                eb_res(i,j) =
     &            max(abs(cfe(2,2)/cfe_J(2,2)),
     &                abs(cfe(2,3)/cfe_J(2,3)),
     &                abs(cfe(2,4)/cfe_J(2,4)),
     &                abs(cfe(3,3)/cfe_J(3,3)),
     &                abs(cfe(3,4)/cfe_J(3,4)),
     &                abs(cfe(4,4)/cfe_J(4,4)))
                kg_res(i,j)=abs(phi1_res/phi1_J)

                ! save pointwise max of constraint violation
                cl_res(i,j)=
     &            max(abs(c_l(1)),abs(c_l(2)),abs(c_l(3)),abs(c_l(4)))

                if (dump.and.first_nan.or.
     &              (ltrace.and.abs(x(i)).lt.0.1.and.abs(y(j)).lt.0.1)
     &             ) then
                  first_nan=.false.
                  write(*,*)
                  write(*,*) 'g_evo_opt: Nan/zero at i,j,Nx,Ny,dx=',
     &                                              i,j,Nx,Ny,dx
                  write(*,*) 'x,y=',x(i),y(j)
                  write(*,*) 'dt,dx,dy=',dt,dx,dy
                  write(*,*) 'x0,y0,rho0=',x0,y0,rho0

                  write(*,*) ' at tn:'
                  write(*,*) ' g0u_tt/x/y:',g0u_tt0,g0u_tx0,g0u_ty0
                  write(*,*) ' g0u_xx/y:',g0u_xx0,g0u_xy0
                  write(*,*) ' g0u_yy:',g0u_yy0
                  write(*,*) ' phi1:',phi1_n(i,j)
                  write(*,*) ' phi np1,n,nm1:',phi1_np1(i,j),
     &                     phi1_n(i,j),phi1_nm1(i,j)
                  write(*,*) ' res J:'
                  write(*,*) ' tt:',efe(1,1),efe_J(1,1)
                  write(*,*) ' tx:',efe(1,2),efe_J(1,2)
                  write(*,*) ' ty:',efe(1,3),efe_J(1,3)
                  write(*,*) ' xx:',efe(2,2),efe_J(2,2)
                  write(*,*) ' xy:',efe(2,3),efe_J(2,3)
                  write(*,*) ' yy:',efe(3,3),efe_J(3,3)
                  write(*,*) ' psi:',efe(4,4),efe_J(4,4)
                  write(*,*) ' phi1:',phi1_res,phi1_J
                end if

              ! (REGION) next-to-ads-bdy points; set by linear interpolation
              else if (chr2(i,j).eq.ex) then
                call interp_from_ads_bdy(gb_tx_np1,x,y,L,i,j,chr,ex,
     &                 Nx,Ny)
                call interp_from_ads_bdy(gb_ty_np1,x,y,L,i,j,chr,ex,
     &                 Nx,Ny)
                call interp_from_ads_bdy(gb_xx_np1,x,y,L,i,j,chr,ex,
     &                 Nx,Ny)
                call interp_from_ads_bdy(gb_xy_np1,x,y,L,i,j,chr,ex,
     &                 Nx,Ny)
                call interp_from_ads_bdy(gb_yy_np1,x,y,L,i,j,chr,ex,
     &                 Nx,Ny)
                call interp_from_ads_bdy(psi_np1,x,y,L,i,j,chr,ex,
     &                 Nx,Ny)
                call interp_from_ads_bdy(Hb_t_np1,x,y,L,i,j,chr,ex,
     &                 Nx,Ny)
                call interp_from_ads_bdy(Hb_x_np1,x,y,L,i,j,chr,ex,
     &                 Nx,Ny)
                call interp_from_ads_bdy(Hb_y_np1,x,y,L,i,j,chr,ex,
     &                 Nx,Ny)
                call interp_from_ads_bdy(phi1_np1,x,y,L,i,j,
     &                    chr,ex,Nx,Ny)
                gb_tt_np1(i,j)=gb_xx_np1(i,j)+gb_yy_np1(i,j)
     &                        +2*psi_np1(i,j)

              ! (REGION) non-interior points; set to zero in prior to applying bcs 
              else 
                gb_tt_np1(i,j) = 0
                gb_tx_np1(i,j) = 0
                gb_ty_np1(i,j) = 0
                gb_xx_np1(i,j) = 0
                gb_xy_np1(i,j) = 0
                gb_yy_np1(i,j) = 0
                psi_np1(i,j) = 0 
                phi1_np1(i,j) = 0 
                gb_res(i,j) = 0

              endif ! (near start of main loop)

            end do
          end do

        end do

        return
        end

c----------------------------------------------------------------------
        logical function is_nan(x)
        implicit none
        real*8 x

        integer is_a_nan

        call check_nan(x,is_a_nan)

        if (is_a_nan.eq.0) then
           is_nan=.false.
        else
           is_nan=.true.
        end if

        return
        end
c----------------------------------------------------------------------
