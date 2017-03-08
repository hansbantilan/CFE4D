c----------------------------------------------------------------------
c in polar coordinates t,x, for x in [0,1]
c
c evolution routine for gb,phi1 computing the residual at the time just prior to update
c
c L below is the AdS length scale
c----------------------------------------------------------------------
        subroutine g_evo_opt(db_res,gb_res,kg_res,cl_res,
     &                       db_txtx_np1,db_txtx_n,db_txtx_nm1,
     &                       db_txty_np1,db_txty_n,db_txty_nm1,
     &                       db_txtz_np1,db_txtz_n,db_txtz_nm1,
     &                       db_txxy_np1,db_txxy_n,db_txxy_nm1,
     &                       db_txxz_np1,db_txxz_n,db_txxz_nm1,
     &                       db_txyz_np1,db_txyz_n,db_txyz_nm1,
     &                       db_tyty_np1,db_tyty_n,db_tyty_nm1,
     &                       db_tytz_np1,db_tytz_n,db_tytz_nm1,
     &                       db_tyxy_np1,db_tyxy_n,db_tyxy_nm1,
     &                       db_tyxz_np1,db_tyxz_n,db_tyxz_nm1,
     &                       db_tyyz_np1,db_tyyz_n,db_tyyz_nm1,
     &                       db_tztz_np1,db_tztz_n,db_tztz_nm1,
     &                       db_tzxy_np1,db_tzxy_n,db_tzxy_nm1,
     &                       db_tzxz_np1,db_tzxz_n,db_tzxz_nm1,
     &                       db_tzyz_np1,db_tzyz_n,db_tzyz_nm1,
     &                       db_xyxy_np1,db_xyxy_n,db_xyxy_nm1,
     &                       db_xyxz_np1,db_xyxz_n,db_xyxz_nm1,
     &                       db_xyyz_np1,db_xyyz_n,db_xyyz_nm1,
     &                       db_xzxz_np1,db_xzxz_n,db_xzxz_nm1,
     &                       db_xzyz_np1,db_xzyz_n,db_xzyz_nm1,
     &                       db_yzyz_np1,db_yzyz_n,db_yzyz_nm1,
     &                       gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                       gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                       gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                       psi_np1,psi_n,psi_nm1,
     &                       Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                       Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                       phi1_np1,phi1_n,phi1_nm1,
     &                       L,x,dt,chr,ex,
     &                       phys_bdy,ghost_width,Nx,
     &                       background,kappa_cd,rho_cd)
        implicit none
        integer Nx
        integer phys_bdy(2),ghost_width(2)
        integer background
        real*8 kappa_cd,rho_cd
        real*8 db_res(Nx),gb_res(Nx),kg_res(Nx),cl_res(Nx)
        real*8 db_txtx_np1(Nx),db_txtx_n(Nx),db_txtx_nm1(Nx) 
        real*8 db_txty_np1(Nx),db_txty_n(Nx),db_txty_nm1(Nx)
        real*8 db_txtz_np1(Nx),db_txtz_n(Nx),db_txtz_nm1(Nx)
        real*8 db_txxy_np1(Nx),db_txxy_n(Nx),db_txxy_nm1(Nx)
        real*8 db_txxz_np1(Nx),db_txxz_n(Nx),db_txxz_nm1(Nx)
        real*8 db_txyz_np1(Nx),db_txyz_n(Nx),db_txyz_nm1(Nx)
        real*8 db_tyty_np1(Nx),db_tyty_n(Nx),db_tyty_nm1(Nx)
        real*8 db_tytz_np1(Nx),db_tytz_n(Nx),db_tytz_nm1(Nx)
        real*8 db_tyxy_np1(Nx),db_tyxy_n(Nx),db_tyxy_nm1(Nx)
        real*8 db_tyxz_np1(Nx),db_tyxz_n(Nx),db_tyxz_nm1(Nx)
        real*8 db_tyyz_np1(Nx),db_tyyz_n(Nx),db_tyyz_nm1(Nx)
        real*8 db_tztz_np1(Nx),db_tztz_n(Nx),db_tztz_nm1(Nx)
        real*8 db_tzxy_np1(Nx),db_tzxy_n(Nx),db_tzxy_nm1(Nx)
        real*8 db_tzxz_np1(Nx),db_tzxz_n(Nx),db_tzxz_nm1(Nx)
        real*8 db_tzyz_np1(Nx),db_tzyz_n(Nx),db_tzyz_nm1(Nx)
        real*8 db_xyxy_np1(Nx),db_xyxy_n(Nx),db_xyxy_nm1(Nx)
        real*8 db_xyxz_np1(Nx),db_xyxz_n(Nx),db_xyxz_nm1(Nx)
        real*8 db_xyyz_np1(Nx),db_xyyz_n(Nx),db_xyyz_nm1(Nx)
        real*8 db_xzxz_np1(Nx),db_xzxz_n(Nx),db_xzxz_nm1(Nx)
        real*8 db_xzyz_np1(Nx),db_xzyz_n(Nx),db_xzyz_nm1(Nx)
        real*8 db_yzyz_np1(Nx),db_yzyz_n(Nx),db_yzyz_nm1(Nx)
        real*8 gb_tt_np1(Nx),gb_tt_n(Nx),gb_tt_nm1(Nx)
        real*8 gb_tx_np1(Nx),gb_tx_n(Nx),gb_tx_nm1(Nx)
        real*8 gb_xx_np1(Nx),gb_xx_n(Nx),gb_xx_nm1(Nx)
        real*8 psi_np1(Nx),psi_n(Nx),psi_nm1(Nx)
        real*8 Hb_t_np1(Nx),Hb_t_n(Nx),Hb_t_nm1(Nx)
        real*8 Hb_x_np1(Nx),Hb_x_n(Nx),Hb_x_nm1(Nx)
        real*8 phi1_np1(Nx),phi1_n(Nx),phi1_nm1(Nx)
        real*8 L
        real*8 x(Nx),dt,chr(Nx),ex

        integer a,b,c,d,e
        integer m,n
        integer rb,i

        real*8 x0
        real*8 dx

        integer max_ghost_width

        real*8 boxx_u(4),boxx_l(4)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 phi1_t,phi1_x
        real*8 phi1_tt,phi1_tx
        real*8 phi1_xx

        real*8 phi1_res,phi1_J
        real*8 phi10_res,phi10_J

        integer is,ie,js,je,ks,ke,is_a_nan

        logical ltrace,is_nan,dump,first_nan
        parameter (ltrace=.false.)
        data first_nan/.true./

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
        ! given z_xx=1, z_xm=(1-x0^2), z_mn=1
        ! h0_ll(a,b)           = z_ab*gb_ab   
        ! h0_uu(a,b)           = inverse(gads+z_ab*gb)^ab - gads^ab
        ! h0_ll_x(a,b,c)       = (z_ab*gb_ab)_,c
        ! h0_ll_xx(a,b,c,d)    = (z_ab*gb_ab)_,cd
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
        ! given z_t=(1-x0^2)^2, z_x=(1-x0^2), z_y=(1-x0^2)^2
        ! A_l(a)     = z_a*Hb_a   
        ! A_l_x(a,b) = (z_a*Hb_a)_,b
        ! 
        ! phi10_x(a)  = phi1_,a     
        !
        ! grad_phi1_sq = g^cd*phi1_,c*phi1_,d
        !
        ! set_ab = 2*phi1_,a*phi1_,b - g_ab*grad_phi1_sq
        ! tr_set= g^cd*set_cd 
        !
        ! cfe(a,b) = residual ... hardcoded expressions (see below)
        !
        ! t,x,y=1,2,3
        ! 
        ! NOTE: g0_ll_xx,gads_ll_xx,h0_ll_xx,cfe,cfe_J
        !       do *NOT* symmetric components filled in
        !
        !--------------------------------------------------------------
        real*8 cfe(4,4),cfe_J(4,4)
        real*8 term1(4,4),term2(4,4),term3(4,4),term4(4,4)
        real*8 term5(4,4),term6(4,4),term7(4,4),term8(4,4)
        real*8 gammagg(4,4,4),gammahh(4,4,4)
        real*8 gammagh(4,4,4),gammahg(4,4,4) 
        real*8 cuuuu(4,4,4,4),dlll(4,4,4)
 
        real*8 ndotc,n_l(4),n_u(4),c_l(4),c_J_l(4)
        real*8 cd_ll(4,4),cd_J_ll(4,4)

        real*8 tr_set,grad_phi1_sq
        
        real*8 g0u_tt_ads0,g0u_xx_ads0,g0u_xy_ads0,g0u_yy_ads0

        real*8 H0_t_ads0,H0_x_ads0,H0_y_ads0

        real*8 dgb_J,ddgb_J,ddgb_J_tx
        real*8 dphi1_J,ddphi1_J,ddphi1_J_tx
        real*8 dc_J

        real*8 lambda4

        real*8 omega
        real*8 dllll(4,4,4,4)
        real*8 dlulu(4,4,4,4)

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,chi,theta)
        !--------------------------------------------------------------
        real*8 d0_llll(4,4,4,4),g0_ll(4,4),g0_uu(4,4)
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
        ! initialize fixed-size variables 
        !--------------------------------------------------------------
        data a,b,c,d,e/0,0,0,0,0/
        data rb,i/0,0/
        data max_ghost_width/0/
        data is,ie,js,je,ks,ke,is_a_nan/0,0,0,0,0,0,0/

        data lambda4/0.0/

        data g0u_tt_ads0,g0u_xx_ads0/0.0,0.0/
        data g0u_xy_ads0,g0u_yy_ads0/0.0,0.0/

        data H0_t_ads0,H0_x_ads0,H0_y_ads0/0.0,0.0,0.0/

        data dgb_J,ddgb_J,ddgb_J_tx/0.0,0.0,0.0/
        data dphi1_J,ddphi1_J/0.0,0.0/
        data dc_J/0.0/

        data dlll/64*0.0/
        data cuuuu/256*0.0/

        data term1,term2/16*0.0,16*0.0/
        data term3,term4/16*0.0,16*0.0/
        data term5,term6/16*0.0,16*0.0/
        data term7,term8/16*0.0,16*0.0/

        data cfe,cfe_J/16*0.0,16*0.0/
        data cd_ll,cd_J_ll/16*0.0,16*0.0/

        data phi1_t,phi1_x/0.0,0.0/
        data phi1_tt,phi1_tx/0.0,0.0/
        data phi1_xx/0.0/

        data x0/0.0/

        data dx/0.0/

        data phi1_res,phi1_J/0.0,0.0/
        data phi10_res,phi10_J/0.0,0.0/

        data n_l,n_u,c_l,c_J_l/4*0.0,4*0.0,4*0.0,4*0.0/

        data gammagg,gammahh/64*0.0,64*0.0/
        data gammagh,gammahg/64*0.0,64*0.0/

        data d0_llll,g0_ll,g0_uu/256*0,16*0.0,16*0.0/
        data gads_ll,gads_uu/16*0.0,16*0.0/
        data h0_ll,h0_uu/16*0.0,16*0.0/
        data gamma_ull/64*0.0/
        data gamma_ull_x/256*0.0/

        data g0_ll_x,g0_uu_x/64*0.0,64*0.0/
        data gads_ll_x,gads_uu_x/64*0.0,64*0.0/
        data h0_ll_x,h0_uu_x/64*0.0,64*0.0/

        data g0_ll_xx/256*0.0/
        data gads_ll_xx/256*0.0/
        data h0_ll_xx/256*0.0/

        data ricci/0.0/
        data ricci_ll,ricci_lu/16*0.0,16*0.0/
        data einstein_ll,set_ll/16*0.0,16*0.0/
        data riemann_ulll/256*0.0/

        data A_l,Hads_l/4*0.0,4*0.0/
        data A_l_x/16*0.0/

        data phi10_x/4*0.0/
        data phi10_xx/16*0.0/

        !--------------------------------------------------------------
        if (ltrace) write(*,*) 'gb_psi_evo ... N=',Nx

        dx=x(2)-x(1)

        ! AdS4D cosmological constant
        !(lambda4=-(n-1)(n-2)/2/L^2) for n=4 dimensional AdS)
        lambda4=-3/L/L

        ! initialize output variables
        do i=1,Nx
          gb_res(i)=0
          kg_res(i)=0
          cl_res(i)=0
          if (chr(i).eq.ex) then
            phi1_np1(i) = 0
          end if
        end do

        ! set index bounds for main loop
        is=2
        ie=Nx-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)

        !(MAIN LOOP) red-black loop through spacetime pts x(i)
        do rb=0,1

          do i=is+mod(rb,2),ie,2
            x0=x(i)

            dump=.false.

            if (ltrace) write(*,*) 'i:',i

            !(REGION) interior points; evolve 
            if (chr(i).ne.ex) then

              ! computes tensors at point i
              call tensor_init(
     &                gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                psi_np1,psi_n,psi_nm1,
     &                Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                phi1_np1,phi1_n,phi1_nm1,
     &                d0_llll,g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &                gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &                h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                A_l,A_l_x,Hads_l,
     &                gamma_ull,gamma_ull_x,
     &                riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                einstein_ll,set_ll,
     &                phi10_x,phi10_xx,
     &                x,dt,chr,L,ex,Nx,i)

              ! background values
              if (background.ne.1) then
                ricci=-12d0/L/L ! conformal gauge source
                omega=sqrt(-3d0/lambda4)/cosh(2*x0/(1-x0**2))**2
              end if

              ! computes auxiliary objects at point i
              do a=1,4
                do b=1,4
                  do c=1,4
                    do d=1,4
                      dlulu(a,b,c,d)=0
                      do m=1,4
                        dlulu(a,b,c,d)=dlulu(a,b,c,d)
     &                                +dllll(a,m,c,d)*g0_uu(m,b)
                      end do
                      do n=1,4
                        dlulu(a,b,c,d)=dlulu(a,b,c,d)
     &                                +dllll(a,b,c,n)*g0_uu(n,d)
                      end do
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
     &                 -( 0.5d0*( h0_uu(1,1)*h0_ll_x(1,1,a)+
     &                            h0_uu(2,2)*h0_ll_x(2,2,a)+
     &                            h0_uu(3,3)*h0_ll_x(3,3,a)+
     &                            h0_uu(4,4)*h0_ll_x(4,4,a)+
     &                         2*(h0_uu(1,2)*h0_ll_x(1,2,a)+
     &                            h0_uu(1,3)*h0_ll_x(1,3,a)+
     &                            h0_uu(1,4)*h0_ll_x(1,4,a)+
     &                            h0_uu(2,3)*h0_ll_x(2,3,a)+
     &                            h0_uu(2,4)*h0_ll_x(2,4,a)+
     &                            h0_uu(3,4)*h0_ll_x(3,4,a)) )
     &                   -1.0d0*( h0_uu(1,1)*h0_ll_x(1,a,1)+
     &                            h0_uu(2,2)*h0_ll_x(2,a,2)+
     &                            h0_uu(3,3)*h0_ll_x(3,a,3)+
     &                            h0_uu(4,4)*h0_ll_x(4,a,4)+
     &                           (h0_uu(1,2)*h0_ll_x(1,a,2)+
     &                            h0_uu(2,1)*h0_ll_x(2,a,1)+
     &                            h0_uu(1,3)*h0_ll_x(1,a,3)+
     &                            h0_uu(3,1)*h0_ll_x(3,a,1)+
     &                            h0_uu(1,4)*h0_ll_x(1,a,4)+
     &                            h0_uu(4,1)*h0_ll_x(4,a,1)+
     &                            h0_uu(2,3)*h0_ll_x(2,a,3)+
     &                            h0_uu(3,2)*h0_ll_x(3,a,2)+
     &                            h0_uu(2,4)*h0_ll_x(2,a,4)+
     &                            h0_uu(4,2)*h0_ll_x(4,a,2)+
     &                            h0_uu(3,4)*h0_ll_x(3,a,4)+
     &                            h0_uu(4,3)*h0_ll_x(4,a,3)) ) )
     &                 -( 0.5d0*( gads_uu(1,1)*h0_ll_x(1,1,a)+
     &                            gads_uu(2,2)*h0_ll_x(2,2,a)+
     &                            gads_uu(3,3)*h0_ll_x(3,3,a)+
     &                            gads_uu(4,4)*h0_ll_x(4,4,a)+
     &                         2*(gads_uu(1,2)*h0_ll_x(1,2,a)+
     &                            gads_uu(1,3)*h0_ll_x(1,3,a)+
     &                            gads_uu(1,4)*h0_ll_x(1,4,a)+
     &                            gads_uu(2,3)*h0_ll_x(2,3,a)+
     &                            gads_uu(2,4)*h0_ll_x(2,4,a)+
     &                            gads_uu(3,4)*h0_ll_x(3,4,a)) )
     &                   -1.0d0*( gads_uu(1,1)*h0_ll_x(1,a,1)+
     &                            gads_uu(2,2)*h0_ll_x(2,a,2)+
     &                            gads_uu(3,3)*h0_ll_x(3,a,3)+
     &                            gads_uu(4,4)*h0_ll_x(4,a,4)+
     &                           (gads_uu(1,2)*h0_ll_x(1,a,2)+
     &                            gads_uu(2,1)*h0_ll_x(2,a,1)+
     &                            gads_uu(1,3)*h0_ll_x(1,a,3)+
     &                            gads_uu(3,1)*h0_ll_x(3,a,1)+
     &                            gads_uu(1,4)*h0_ll_x(1,a,4)+
     &                            gads_uu(4,1)*h0_ll_x(4,a,1)+
     &                            gads_uu(2,3)*h0_ll_x(2,a,3)+
     &                            gads_uu(3,2)*h0_ll_x(3,a,2)+
     &                            gads_uu(2,4)*h0_ll_x(2,a,4)+
     &                            gads_uu(4,2)*h0_ll_x(4,a,2)+
     &                            gads_uu(3,4)*h0_ll_x(3,a,4)+
     &                            gads_uu(4,3)*h0_ll_x(4,a,3)) ) )
     &                 -( 0.5d0*( h0_uu(1,1)*gads_ll_x(1,1,a)+
     &                            h0_uu(2,2)*gads_ll_x(2,2,a)+
     &                            h0_uu(3,3)*gads_ll_x(3,3,a)+
     &                            h0_uu(4,4)*gads_ll_x(4,4,a)+
     &                         2*(h0_uu(1,2)*gads_ll_x(1,2,a)+
     &                            h0_uu(1,3)*gads_ll_x(1,3,a)+
     &                            h0_uu(1,4)*gads_ll_x(1,4,a)+
     &                            h0_uu(2,3)*gads_ll_x(2,3,a)+
     &                            h0_uu(2,4)*gads_ll_x(2,4,a)+
     &                            h0_uu(3,4)*gads_ll_x(3,4,a)) )
     &                   -1.0d0*( h0_uu(1,1)*gads_ll_x(1,a,1)+
     &                            h0_uu(2,2)*gads_ll_x(2,a,2)+
     &                            h0_uu(3,3)*gads_ll_x(3,a,3)+
     &                            h0_uu(4,4)*gads_ll_x(4,a,4)+
     &                           (h0_uu(1,2)*gads_ll_x(1,a,2)+
     &                            h0_uu(2,1)*gads_ll_x(2,a,1)+
     &                            h0_uu(1,3)*gads_ll_x(1,a,3)+
     &                            h0_uu(3,1)*gads_ll_x(3,a,1)+
     &                            h0_uu(1,4)*gads_ll_x(1,a,4)+
     &                            h0_uu(4,1)*gads_ll_x(4,a,1)+
     &                            h0_uu(2,3)*gads_ll_x(2,a,3)+
     &                            h0_uu(3,2)*gads_ll_x(3,a,2)+
     &                            h0_uu(2,4)*gads_ll_x(2,a,4)+
     &                            h0_uu(4,2)*gads_ll_x(4,a,2)+
     &                            h0_uu(3,4)*gads_ll_x(3,a,4)+
     &                            h0_uu(4,3)*gads_ll_x(4,a,3)) ) )
              end do

              n_l(1)=-1/sqrt(-g0_uu(1,1))
              do a=1,4
                n_u(a)=n_l(1)*g0_uu(a,1)+
     &                 n_l(2)*g0_uu(a,2)+
     &                 n_l(3)*g0_uu(a,3)+
     &                 n_l(4)*g0_uu(a,4)
              end do

              ndotc  =n_u(1)*c_l(1)+
     &                n_u(2)*c_l(2)+
     &                n_u(3)*c_l(3)+
     &                n_u(4)*c_l(4)

              grad_phi1_sq=phi10_x(1)*phi10_x(1)*g0_uu(1,1)+
     &                     phi10_x(2)*phi10_x(2)*g0_uu(2,2)+
     &                     phi10_x(3)*phi10_x(3)*g0_uu(3,3)+
     &                     phi10_x(4)*phi10_x(4)*g0_uu(4,4)+
     &                  2*(phi10_x(1)*phi10_x(2)*g0_uu(1,2)+
     &                     phi10_x(1)*phi10_x(3)*g0_uu(1,3)+
     &                     phi10_x(1)*phi10_x(4)*g0_uu(1,4)+
     &                     phi10_x(2)*phi10_x(3)*g0_uu(2,3)+
     &                     phi10_x(2)*phi10_x(4)*g0_uu(2,4)+
     &                     phi10_x(3)*phi10_x(4)*g0_uu(3,4))
         
              do a=1,4
                do b=1,4
                  set_ll(a,b)=phi10_x(a)*phi10_x(b)
     &                       -g0_ll(a,b)*(grad_phi1_sq/2)
                end do
              end do

              tr_set =set_ll(1,1)*g0_uu(1,1)+
     &                set_ll(2,2)*g0_uu(2,2)+
     &                set_ll(3,3)*g0_uu(3,3)+
     &                set_ll(4,4)*g0_uu(4,4)+
     &             2*(set_ll(1,2)*g0_uu(1,2)+
     &                set_ll(1,3)*g0_uu(1,3)+
     &                set_ll(1,4)*g0_uu(1,4)+
     &                set_ll(2,3)*g0_uu(2,3)+
     &                set_ll(2,4)*g0_uu(2,4)+
     &                set_ll(3,4)*g0_uu(3,4))

              !---------------------------------------------------------------- 
              ! cfe_ab =   term1_ab + term2_ab + term3_ab + term4_ab 
              !          + term5_ab 
              ! for
              ! 
              ! d_term1_abcd = -g^mn d_abcd,mn  
              ! d_term2_abcd = -2 omega d_a^m_b^n d_cdmn
              ! d_term3_abcd = +2 omega d_a^m_d^n d_cnbm
              ! d_term4_abcd = -2 omega d_a^m_c^n d_ambn
              ! d_term5_abcd = +1/2 d_abcd R
              !
              !----------------------------------------------------------------



              !---------------------------------------------------------------- 
              ! Analytically remove the pure AdS terms in the EFEs,
              ! denoted ->, so 
              ! efe_ab =   term1_ab + term2_ab + term3_ab + term4_ab 
              !          + term5_ab + term6_ab + term7_ab + term8_ab 
              !          - 8*PI*(se_ab-1/3*tr(se)*g_ab)
              ! for
              ! 
              ! term1_ab = -1/2 g^cd g_ab,cd  
              !          ->-1/2 (h^cd h_ab,cd + gads^cd h_ab,cd + h^cd
              !          gads_ab,cd) 
              ! term2_ab = -1/2 g^cd,a g_bc,d
              !          ->-1/2 (h^cd,a h_bc,d + gads^cd,a h_bc,d +
              !          h^cd,a gads_bc,d)
              ! term3_ab = -1/2 g^cd,b g_ac,d
              !          ->-1/2 (h^cd,b h_ac,d + gads^cd,b h_ac,d +
              !          h^cd,b gads_ac,d)
              ! term4_ab = -1/2 H_a,b
              !          ->-1/2 Hb_a,b
              ! term5_ab = -1/2 H_b,a
              !          ->-1/2 Hb_b,a
              ! term6_ab = H_c G^c_ab
              !          ->Hb_c Ghh^c_ab + Hb_c Ggh^c_ab + Hb_c
              !          Ghg^c_ab  
              ! term7_ab = -G^c_db G^d_ca
              !          ->-(  Ghh^c_db Ghh^d_ca + Ggg^c_db Ghh^d_ca
              !              + Ggg^c_db Ggh^d_ca + Ggg^c_db Ghg^d_ca
              !              + Ghh^c_db Ggg^d_ca + Ghh^c_db Ggh^d_ca
              !              + Ghh^c_db Ghg^d_ca + Ggh^c_db Ggg^d_ca
              !              + Ggh^c_db Ghh^d_ca + Ggh^c_db Ggh^d_ca
              !              + Ggh^c_db Ghg^d_ca + Ghg^c_db Ggg^d_ca
              !              + Ghg^c_db Ghh^d_ca + Ghg^c_db Ggh^d_ca
              !              + Ghg^c_db Ghg^d_ca  )
              ! term8_ab = - lambda4 g_ab
              !          ->- lambda4 h_ab
              !
              ! where G   = guu(g_ll_x-g_ll_x+g_ll_x)
              ! where Ggg = gadsuu(gads_ll_x-gads_ll_x+gads_ll_x)
              ! where Ghh = huu(h_ll_x-h_ll_x+h_ll_x)
              ! where Ggh = gadsuu(h_ll_x-h_ll_x+h_ll_x)
              ! where Ghg = huu(gads_ll_x-gads_ll_x+gads_ll_x)
              !
              !----------------------------------------------------------------
              do a=1,3
                do b=a,3
                  term1(a,b)=-0.5d0*(                             
     &                          h0_uu(1,1)*h0_ll_xx(a,b,1,1)+
     &                          h0_uu(2,2)*h0_ll_xx(a,b,2,2)+
     &                          h0_uu(3,3)*h0_ll_xx(a,b,3,3)+
     &                          h0_uu(4,4)*h0_ll_xx(a,b,4,4)+
     &                       2*(h0_uu(1,2)*h0_ll_xx(a,b,1,2)+
     &                          h0_uu(1,3)*h0_ll_xx(a,b,1,3)+
     &                          h0_uu(1,4)*h0_ll_xx(a,b,1,4)+
     &                          h0_uu(2,3)*h0_ll_xx(a,b,2,3)+
     &                          h0_uu(2,4)*h0_ll_xx(a,b,2,4)+
     &                          h0_uu(3,4)*h0_ll_xx(a,b,3,4))
     &                       +
     &                          gads_uu(1,1)*h0_ll_xx(a,b,1,1)+
     &                          gads_uu(2,2)*h0_ll_xx(a,b,2,2)+
     &                          gads_uu(3,3)*h0_ll_xx(a,b,3,3)+
     &                          gads_uu(4,4)*h0_ll_xx(a,b,4,4)+
     &                       2*(gads_uu(1,2)*h0_ll_xx(a,b,1,2)+
     &                          gads_uu(1,3)*h0_ll_xx(a,b,1,3)+
     &                          gads_uu(1,4)*h0_ll_xx(a,b,1,4)+
     &                          gads_uu(2,3)*h0_ll_xx(a,b,2,3)+
     &                          gads_uu(2,4)*h0_ll_xx(a,b,2,4)+
     &                          gads_uu(3,4)*h0_ll_xx(a,b,3,4))
     &                       +
     &                          h0_uu(1,1)*gads_ll_xx(a,b,1,1)+
     &                          h0_uu(2,2)*gads_ll_xx(a,b,2,2)+
     &                          h0_uu(3,3)*gads_ll_xx(a,b,3,3)+
     &                          h0_uu(4,4)*gads_ll_xx(a,b,4,4)+
     &                       2*(h0_uu(1,2)*gads_ll_xx(a,b,1,2)+
     &                          h0_uu(1,3)*gads_ll_xx(a,b,1,3)+
     &                          h0_uu(1,4)*gads_ll_xx(a,b,1,4)+
     &                          h0_uu(2,3)*gads_ll_xx(a,b,2,3)+
     &                          h0_uu(2,4)*gads_ll_xx(a,b,2,4)+
     &                          h0_uu(3,4)*gads_ll_xx(a,b,3,4))
     &                              )
     &
                  term2(a,b)=-0.5d0*(                             
     &                          h0_uu_x(1,1,a)* h0_ll_x(b,1,1) +
     &                          h0_uu_x(1,2,a)*(h0_ll_x(b,1,2) +
     &                                          h0_ll_x(b,2,1))+
     &                          h0_uu_x(1,3,a)*(h0_ll_x(b,1,3) +
     &                                          h0_ll_x(b,3,1))+
     &                          h0_uu_x(1,4,a)*(h0_ll_x(b,1,4) +
     &                                          h0_ll_x(b,4,1))+
     &                          h0_uu_x(2,2,a)* h0_ll_x(b,2,2) +
     &                          h0_uu_x(2,3,a)*(h0_ll_x(b,2,3) +
     &                                          h0_ll_x(b,3,2))+
     &                          h0_uu_x(2,4,a)*(h0_ll_x(b,2,4) +
     &                                          h0_ll_x(b,4,2))+
     &                          h0_uu_x(3,3,a)* h0_ll_x(b,3,3) +
     &                          h0_uu_x(3,4,a)*(h0_ll_x(b,3,4) +
     &                                          h0_ll_x(b,4,3))+
     &                          h0_uu_x(4,4,a)* h0_ll_x(b,4,4) 
     &                       +
     &                          gads_uu_x(1,1,a)* h0_ll_x(b,1,1) +
     &                          gads_uu_x(1,2,a)*(h0_ll_x(b,1,2) +
     &                                            h0_ll_x(b,2,1))+
     &                          gads_uu_x(1,3,a)*(h0_ll_x(b,1,3) +
     &                                            h0_ll_x(b,3,1))+
     &                          gads_uu_x(1,4,a)*(h0_ll_x(b,1,4) +
     &                                            h0_ll_x(b,4,1))+
     &                          gads_uu_x(2,2,a)* h0_ll_x(b,2,2) +
     &                          gads_uu_x(2,3,a)*(h0_ll_x(b,2,3) +
     &                                            h0_ll_x(b,3,2))+
     &                          gads_uu_x(2,4,a)*(h0_ll_x(b,2,4) +
     &                                            h0_ll_x(b,4,2))+
     &                          gads_uu_x(3,3,a)* h0_ll_x(b,3,3) +
     &                          gads_uu_x(3,4,a)*(h0_ll_x(b,3,4) +
     &                                            h0_ll_x(b,4,3))+
     &                          gads_uu_x(4,4,a)* h0_ll_x(b,4,4) 
     &                       +
     &                          h0_uu_x(1,1,a)* gads_ll_x(b,1,1) +
     &                          h0_uu_x(1,2,a)*(gads_ll_x(b,1,2) + 
     &                                          gads_ll_x(b,2,1))+ 
     &                          h0_uu_x(1,3,a)*(gads_ll_x(b,1,3) + 
     &                                          gads_ll_x(b,3,1))+ 
     &                          h0_uu_x(1,4,a)*(gads_ll_x(b,1,4) +
     &                                          gads_ll_x(b,4,1))+
     &                          h0_uu_x(2,2,a)* gads_ll_x(b,2,2) +
     &                          h0_uu_x(2,3,a)*(gads_ll_x(b,2,3) +
     &                                          gads_ll_x(b,3,2))+
     &                          h0_uu_x(2,4,a)*(gads_ll_x(b,2,4) +
     &                                          gads_ll_x(b,4,2))+
     &                          h0_uu_x(3,3,a)* gads_ll_x(b,3,3) +
     &                          h0_uu_x(3,4,a)*(gads_ll_x(b,3,4) +
     &                                          gads_ll_x(b,4,3))+
     &                          h0_uu_x(4,4,a)* gads_ll_x(b,4,4) 
     &                            )
     &
                  term3(a,b)=-0.5d0*(                            
     &                          h0_uu_x(1,1,b)* h0_ll_x(a,1,1) +
     &                          h0_uu_x(1,2,b)*(h0_ll_x(a,1,2) +
     &                                          h0_ll_x(a,2,1))+
     &                          h0_uu_x(1,3,b)*(h0_ll_x(a,1,3) +
     &                                          h0_ll_x(a,3,1))+
     &                          h0_uu_x(1,4,b)*(h0_ll_x(a,1,4) +
     &                                          h0_ll_x(a,4,1))+
     &                          h0_uu_x(2,2,b)* h0_ll_x(a,2,2) +
     &                          h0_uu_x(2,3,b)*(h0_ll_x(a,2,3) +
     &                                          h0_ll_x(a,3,2))+
     &                          h0_uu_x(2,4,b)*(h0_ll_x(a,2,4) +
     &                                          h0_ll_x(a,4,2))+
     &                          h0_uu_x(3,3,b)* h0_ll_x(a,3,3) +
     &                          h0_uu_x(3,4,b)*(h0_ll_x(a,3,4) +
     &                                          h0_ll_x(a,4,3))+
     &                          h0_uu_x(4,4,b)* h0_ll_x(a,4,4) 
     &                       +
     &                          gads_uu_x(1,1,b)* h0_ll_x(a,1,1) +
     &                          gads_uu_x(1,2,b)*(h0_ll_x(a,1,2) +
     &                                            h0_ll_x(a,2,1))+
     &                          gads_uu_x(1,3,b)*(h0_ll_x(a,1,3) +
     &                                            h0_ll_x(a,3,1))+
     &                          gads_uu_x(1,4,b)*(h0_ll_x(a,1,4) +
     &                                            h0_ll_x(a,4,1))+
     &                          gads_uu_x(2,2,b)* h0_ll_x(a,2,2) +
     &                          gads_uu_x(2,3,b)*(h0_ll_x(a,2,3) +
     &                                            h0_ll_x(a,3,2))+
     &                          gads_uu_x(2,4,b)*(h0_ll_x(a,2,4) +
     &                                            h0_ll_x(a,4,2))+
     &                          gads_uu_x(3,3,b)* h0_ll_x(a,3,3) +
     &                          gads_uu_x(3,4,b)*(h0_ll_x(a,3,4) +
     &                                            h0_ll_x(a,4,3))+
     &                          gads_uu_x(4,4,b)* h0_ll_x(a,4,4) 
     &                       +
     &                          h0_uu_x(1,1,b)* gads_ll_x(a,1,1) +
     &                          h0_uu_x(1,2,b)*(gads_ll_x(a,1,2) +  
     &                                          gads_ll_x(a,2,1))+  
     &                          h0_uu_x(1,3,b)*(gads_ll_x(a,1,3) +  
     &                                          gads_ll_x(a,3,1))+   
     &                          h0_uu_x(1,4,b)*(gads_ll_x(a,1,4) +
     &                                          gads_ll_x(a,4,1))+
     &                          h0_uu_x(2,2,b)* gads_ll_x(a,2,2) +
     &                          h0_uu_x(2,3,b)*(gads_ll_x(a,2,3) +
     &                                          gads_ll_x(a,3,2))+
     &                          h0_uu_x(2,4,b)*(gads_ll_x(a,2,4) +
     &                                          gads_ll_x(a,4,2))+
     &                          h0_uu_x(3,3,b)* gads_ll_x(a,3,3) +
     &                          h0_uu_x(3,4,b)*(gads_ll_x(a,3,4) +
     &                                          gads_ll_x(a,4,3))+
     &                          h0_uu_x(4,4,b)* gads_ll_x(a,4,4) 
     &                            )
     &
                  term4(a,b)=-0.5d0*A_l_x(a,b)                  
     &
                  term5(a,b)=-0.5d0*A_l_x(b,a)           
     &
                  term6(a,b)=     (           
     &                          Hads_l(1)*gammahh(1,a,b)+      
     &                          Hads_l(2)*gammahh(2,a,b)+
     &                          Hads_l(3)*gammahh(3,a,b)+
     &                          Hads_l(4)*gammahh(4,a,b)
     &                       +
     &                          Hads_l(1)*gammagh(1,a,b)+
     &                          Hads_l(2)*gammagh(2,a,b)+
     &                          Hads_l(3)*gammagh(3,a,b)+
     &                          Hads_l(4)*gammagh(4,a,b)
     &                       +
     &                          Hads_l(1)*gammahg(1,a,b)+
     &                          Hads_l(2)*gammahg(2,a,b)+
     &                          Hads_l(3)*gammahg(3,a,b)+
     &                          Hads_l(4)*gammahg(4,a,b)
     &                       +                  
     &                          A_l(1)*gammagg(1,a,b)  +
     &                          A_l(2)*gammagg(2,a,b)  +
     &                          A_l(3)*gammagg(3,a,b)  +
     &                          A_l(4)*gammagg(4,a,b)  
     &                       +
     &                          A_l(1)*gammahh(1,a,b)  +
     &                          A_l(2)*gammahh(2,a,b)  +
     &                          A_l(3)*gammahh(3,a,b)  +
     &                          A_l(4)*gammahh(4,a,b)  
     &                       +
     &                          A_l(1)*gammagh(1,a,b)  +
     &                          A_l(2)*gammagh(2,a,b)  +
     &                          A_l(3)*gammagh(3,a,b)  +
     &                          A_l(4)*gammagh(4,a,b)  
     &                       +
     &                          A_l(1)*gammahg(1,a,b)  +
     &                          A_l(2)*gammahg(2,a,b)  +
     &                          A_l(3)*gammahg(3,a,b)  +
     &                          A_l(4)*gammahg(4,a,b)  
     &                            ) 
     &                           
                  term7(a,b)=    -(
     &                          gammahh(1,1,b)*gammahh(1,1,a)+
     &                          gammahh(1,2,b)*gammahh(2,1,a)+
     &                          gammahh(1,3,b)*gammahh(3,1,a)+
     &                          gammahh(1,4,b)*gammahh(4,1,a)+
     &                          gammahh(2,1,b)*gammahh(1,2,a)+
     &                          gammahh(2,2,b)*gammahh(2,2,a)+
     &                          gammahh(2,3,b)*gammahh(3,2,a)+
     &                          gammahh(2,4,b)*gammahh(4,2,a)+
     &                          gammahh(3,1,b)*gammahh(1,3,a)+
     &                          gammahh(3,2,b)*gammahh(2,3,a)+
     &                          gammahh(3,3,b)*gammahh(3,3,a)+
     &                          gammahh(3,4,b)*gammahh(4,3,a)+
     &                          gammahh(4,1,b)*gammahh(1,4,a)+
     &                          gammahh(4,2,b)*gammahh(2,4,a)+
     &                          gammahh(4,3,b)*gammahh(3,4,a)+
     &                          gammahh(4,4,b)*gammahh(4,4,a)
     &                       +
     &                          gammagg(1,1,b)*gammahh(1,1,a)+
     &                          gammagg(1,2,b)*gammahh(2,1,a)+
     &                          gammagg(1,3,b)*gammahh(3,1,a)+
     &                          gammagg(1,4,b)*gammahh(4,1,a)+
     &                          gammagg(2,1,b)*gammahh(1,2,a)+
     &                          gammagg(2,2,b)*gammahh(2,2,a)+
     &                          gammagg(2,3,b)*gammahh(3,2,a)+
     &                          gammagg(2,4,b)*gammahh(4,2,a)+
     &                          gammagg(3,1,b)*gammahh(1,3,a)+
     &                          gammagg(3,2,b)*gammahh(2,3,a)+
     &                          gammagg(3,3,b)*gammahh(3,3,a)+
     &                          gammagg(3,4,b)*gammahh(4,3,a)+
     &                          gammagg(4,1,b)*gammahh(1,4,a)+
     &                          gammagg(4,2,b)*gammahh(2,4,a)+
     &                          gammagg(4,3,b)*gammahh(3,4,a)+
     &                          gammagg(4,4,b)*gammahh(4,4,a)
     &                       +
     &                          gammagg(1,1,b)*gammagh(1,1,a)+
     &                          gammagg(1,2,b)*gammagh(2,1,a)+
     &                          gammagg(1,3,b)*gammagh(3,1,a)+
     &                          gammagg(1,4,b)*gammagh(4,1,a)+
     &                          gammagg(2,1,b)*gammagh(1,2,a)+
     &                          gammagg(2,2,b)*gammagh(2,2,a)+
     &                          gammagg(2,3,b)*gammagh(3,2,a)+
     &                          gammagg(2,4,b)*gammagh(4,2,a)+
     &                          gammagg(3,1,b)*gammagh(1,3,a)+
     &                          gammagg(3,2,b)*gammagh(2,3,a)+
     &                          gammagg(3,3,b)*gammagh(3,3,a)+
     &                          gammagg(3,4,b)*gammagh(4,3,a)+
     &                          gammagg(4,1,b)*gammagh(1,4,a)+
     &                          gammagg(4,2,b)*gammagh(2,4,a)+
     &                          gammagg(4,3,b)*gammagh(3,4,a)+
     &                          gammagg(4,4,b)*gammagh(4,4,a)
     &                       +
     &                          gammagg(1,1,b)*gammahg(1,1,a)+
     &                          gammagg(1,2,b)*gammahg(2,1,a)+
     &                          gammagg(1,3,b)*gammahg(3,1,a)+
     &                          gammagg(1,4,b)*gammahg(4,1,a)+
     &                          gammagg(2,1,b)*gammahg(1,2,a)+
     &                          gammagg(2,2,b)*gammahg(2,2,a)+
     &                          gammagg(2,3,b)*gammahg(3,2,a)+
     &                          gammagg(2,4,b)*gammahg(4,2,a)+
     &                          gammagg(3,1,b)*gammahg(1,3,a)+
     &                          gammagg(3,2,b)*gammahg(2,3,a)+
     &                          gammagg(3,3,b)*gammahg(3,3,a)+
     &                          gammagg(3,4,b)*gammahg(4,3,a)+
     &                          gammagg(4,1,b)*gammahg(1,4,a)+
     &                          gammagg(4,2,b)*gammahg(2,4,a)+
     &                          gammagg(4,3,b)*gammahg(3,4,a)+
     &                          gammagg(4,4,b)*gammahg(4,4,a)
     &                       +
     &                          gammahh(1,1,b)*gammagg(1,1,a)+
     &                          gammahh(1,2,b)*gammagg(2,1,a)+
     &                          gammahh(1,3,b)*gammagg(3,1,a)+
     &                          gammahh(1,4,b)*gammagg(4,1,a)+
     &                          gammahh(2,1,b)*gammagg(1,2,a)+
     &                          gammahh(2,2,b)*gammagg(2,2,a)+
     &                          gammahh(2,3,b)*gammagg(3,2,a)+
     &                          gammahh(2,4,b)*gammagg(4,2,a)+
     &                          gammahh(3,1,b)*gammagg(1,3,a)+
     &                          gammahh(3,2,b)*gammagg(2,3,a)+
     &                          gammahh(3,3,b)*gammagg(3,3,a)+
     &                          gammahh(3,4,b)*gammagg(4,3,a)+
     &                          gammahh(4,1,b)*gammagg(1,4,a)+
     &                          gammahh(4,2,b)*gammagg(2,4,a)+
     &                          gammahh(4,3,b)*gammagg(3,4,a)+
     &                          gammahh(4,4,b)*gammagg(4,4,a)
     &                       +
     &                          gammahh(1,1,b)*gammagh(1,1,a)+
     &                          gammahh(1,2,b)*gammagh(2,1,a)+
     &                          gammahh(1,3,b)*gammagh(3,1,a)+
     &                          gammahh(1,4,b)*gammagh(4,1,a)+
     &                          gammahh(2,1,b)*gammagh(1,2,a)+
     &                          gammahh(2,2,b)*gammagh(2,2,a)+
     &                          gammahh(2,3,b)*gammagh(3,2,a)+
     &                          gammahh(2,4,b)*gammagh(4,2,a)+
     &                          gammahh(3,1,b)*gammagh(1,3,a)+
     &                          gammahh(3,2,b)*gammagh(2,3,a)+
     &                          gammahh(3,3,b)*gammagh(3,3,a)+
     &                          gammahh(3,4,b)*gammagh(4,3,a)+
     &                          gammahh(4,1,b)*gammagh(1,4,a)+
     &                          gammahh(4,2,b)*gammagh(2,4,a)+
     &                          gammahh(4,3,b)*gammagh(3,4,a)+
     &                          gammahh(4,4,b)*gammagh(4,4,a)
     &                       +
     &                          gammahh(1,1,b)*gammahg(1,1,a)+
     &                          gammahh(1,2,b)*gammahg(2,1,a)+
     &                          gammahh(1,3,b)*gammahg(3,1,a)+
     &                          gammahh(1,4,b)*gammahg(4,1,a)+
     &                          gammahh(2,1,b)*gammahg(1,2,a)+
     &                          gammahh(2,2,b)*gammahg(2,2,a)+
     &                          gammahh(2,3,b)*gammahg(3,2,a)+
     &                          gammahh(2,4,b)*gammahg(4,2,a)+
     &                          gammahh(3,1,b)*gammahg(1,3,a)+
     &                          gammahh(3,2,b)*gammahg(2,3,a)+
     &                          gammahh(3,3,b)*gammahg(3,3,a)+
     &                          gammahh(3,4,b)*gammahg(4,3,a)+
     &                          gammahh(4,1,b)*gammahg(1,4,a)+
     &                          gammahh(4,2,b)*gammahg(2,4,a)+
     &                          gammahh(4,3,b)*gammahg(3,4,a)+
     &                          gammahh(4,4,b)*gammahg(4,4,a)
     &                       +
     &                          gammagh(1,1,b)*gammagg(1,1,a)+
     &                          gammagh(1,2,b)*gammagg(2,1,a)+
     &                          gammagh(1,3,b)*gammagg(3,1,a)+
     &                          gammagh(1,4,b)*gammagg(4,1,a)+
     &                          gammagh(2,1,b)*gammagg(1,2,a)+
     &                          gammagh(2,2,b)*gammagg(2,2,a)+
     &                          gammagh(2,3,b)*gammagg(3,2,a)+
     &                          gammagh(2,4,b)*gammagg(4,2,a)+
     &                          gammagh(3,1,b)*gammagg(1,3,a)+
     &                          gammagh(3,2,b)*gammagg(2,3,a)+
     &                          gammagh(3,3,b)*gammagg(3,3,a)+
     &                          gammagh(3,4,b)*gammagg(4,3,a)+
     &                          gammagh(4,1,b)*gammagg(1,4,a)+
     &                          gammagh(4,2,b)*gammagg(2,4,a)+
     &                          gammagh(4,3,b)*gammagg(3,4,a)+
     &                          gammagh(4,4,b)*gammagg(4,4,a)
     &                       +
     &                          gammagh(1,1,b)*gammahh(1,1,a)+
     &                          gammagh(1,2,b)*gammahh(2,1,a)+
     &                          gammagh(1,3,b)*gammahh(3,1,a)+
     &                          gammagh(1,4,b)*gammahh(4,1,a)+
     &                          gammagh(2,1,b)*gammahh(1,2,a)+
     &                          gammagh(2,2,b)*gammahh(2,2,a)+
     &                          gammagh(2,3,b)*gammahh(3,2,a)+
     &                          gammagh(2,4,b)*gammahh(4,2,a)+
     &                          gammagh(3,1,b)*gammahh(1,3,a)+
     &                          gammagh(3,2,b)*gammahh(2,3,a)+
     &                          gammagh(3,3,b)*gammahh(3,3,a)+
     &                          gammagh(3,4,b)*gammahh(4,3,a)+
     &                          gammagh(4,1,b)*gammahh(1,4,a)+
     &                          gammagh(4,2,b)*gammahh(2,4,a)+
     &                          gammagh(4,3,b)*gammahh(3,4,a)+
     &                          gammagh(4,4,b)*gammahh(4,4,a)
     &                       +
     &                          gammagh(1,1,b)*gammagh(1,1,a)+
     &                          gammagh(1,2,b)*gammagh(2,1,a)+
     &                          gammagh(1,3,b)*gammagh(3,1,a)+
     &                          gammagh(1,4,b)*gammagh(4,1,a)+
     &                          gammagh(2,1,b)*gammagh(1,2,a)+
     &                          gammagh(2,2,b)*gammagh(2,2,a)+
     &                          gammagh(2,3,b)*gammagh(3,2,a)+
     &                          gammagh(2,4,b)*gammagh(4,2,a)+
     &                          gammagh(3,1,b)*gammagh(1,3,a)+
     &                          gammagh(3,2,b)*gammagh(2,3,a)+
     &                          gammagh(3,3,b)*gammagh(3,3,a)+
     &                          gammagh(3,4,b)*gammagh(4,3,a)+
     &                          gammagh(4,1,b)*gammagh(1,4,a)+
     &                          gammagh(4,2,b)*gammagh(2,4,a)+
     &                          gammagh(4,3,b)*gammagh(3,4,a)+
     &                          gammagh(4,4,b)*gammagh(4,4,a)
     &                       +
     &                          gammagh(1,1,b)*gammahg(1,1,a)+
     &                          gammagh(1,2,b)*gammahg(2,1,a)+
     &                          gammagh(1,3,b)*gammahg(3,1,a)+
     &                          gammagh(1,4,b)*gammahg(4,1,a)+
     &                          gammagh(2,1,b)*gammahg(1,2,a)+
     &                          gammagh(2,2,b)*gammahg(2,2,a)+
     &                          gammagh(2,3,b)*gammahg(3,2,a)+
     &                          gammagh(2,4,b)*gammahg(4,2,a)+
     &                          gammagh(3,1,b)*gammahg(1,3,a)+
     &                          gammagh(3,2,b)*gammahg(2,3,a)+
     &                          gammagh(3,3,b)*gammahg(3,3,a)+
     &                          gammagh(3,4,b)*gammahg(4,3,a)+
     &                          gammagh(4,1,b)*gammahg(1,4,a)+
     &                          gammagh(4,2,b)*gammahg(2,4,a)+
     &                          gammagh(4,3,b)*gammahg(3,4,a)+
     &                          gammagh(4,4,b)*gammahg(4,4,a)
     &                       +
     &                          gammahg(1,1,b)*gammagg(1,1,a)+
     &                          gammahg(1,2,b)*gammagg(2,1,a)+
     &                          gammahg(1,3,b)*gammagg(3,1,a)+
     &                          gammahg(1,4,b)*gammagg(4,1,a)+
     &                          gammahg(2,1,b)*gammagg(1,2,a)+
     &                          gammahg(2,2,b)*gammagg(2,2,a)+
     &                          gammahg(2,3,b)*gammagg(3,2,a)+
     &                          gammahg(2,4,b)*gammagg(4,2,a)+
     &                          gammahg(3,1,b)*gammagg(1,3,a)+
     &                          gammahg(3,2,b)*gammagg(2,3,a)+
     &                          gammahg(3,3,b)*gammagg(3,3,a)+
     &                          gammahg(3,4,b)*gammagg(4,3,a)+
     &                          gammahg(4,1,b)*gammagg(1,4,a)+
     &                          gammahg(4,2,b)*gammagg(2,4,a)+
     &                          gammahg(4,3,b)*gammagg(3,4,a)+
     &                          gammahg(4,4,b)*gammagg(4,4,a)
     &                       +
     &                          gammahg(1,1,b)*gammahh(1,1,a)+
     &                          gammahg(1,2,b)*gammahh(2,1,a)+
     &                          gammahg(1,3,b)*gammahh(3,1,a)+
     &                          gammahg(1,4,b)*gammahh(4,1,a)+
     &                          gammahg(2,1,b)*gammahh(1,2,a)+
     &                          gammahg(2,2,b)*gammahh(2,2,a)+
     &                          gammahg(2,3,b)*gammahh(3,2,a)+
     &                          gammahg(2,4,b)*gammahh(4,2,a)+
     &                          gammahg(3,1,b)*gammahh(1,3,a)+
     &                          gammahg(3,2,b)*gammahh(2,3,a)+
     &                          gammahg(3,3,b)*gammahh(3,3,a)+
     &                          gammahg(3,4,b)*gammahh(4,3,a)+
     &                          gammahg(4,1,b)*gammahh(1,4,a)+
     &                          gammahg(4,2,b)*gammahh(2,4,a)+
     &                          gammahg(4,3,b)*gammahh(3,4,a)+
     &                          gammahg(4,4,b)*gammahh(4,4,a)
     &                       +
     &                          gammahg(1,1,b)*gammagh(1,1,a)+
     &                          gammahg(1,2,b)*gammagh(2,1,a)+
     &                          gammahg(1,3,b)*gammagh(3,1,a)+
     &                          gammahg(1,4,b)*gammagh(4,1,a)+
     &                          gammahg(2,1,b)*gammagh(1,2,a)+
     &                          gammahg(2,2,b)*gammagh(2,2,a)+
     &                          gammahg(2,3,b)*gammagh(3,2,a)+
     &                          gammahg(2,4,b)*gammagh(4,2,a)+
     &                          gammahg(3,1,b)*gammagh(1,3,a)+
     &                          gammahg(3,2,b)*gammagh(2,3,a)+
     &                          gammahg(3,3,b)*gammagh(3,3,a)+
     &                          gammahg(3,4,b)*gammagh(4,3,a)+
     &                          gammahg(4,1,b)*gammagh(1,4,a)+
     &                          gammahg(4,2,b)*gammagh(2,4,a)+
     &                          gammahg(4,3,b)*gammagh(3,4,a)+
     &                          gammahg(4,4,b)*gammagh(4,4,a)
     &                       +
     &                          gammahg(1,1,b)*gammahg(1,1,a)+
     &                          gammahg(1,2,b)*gammahg(2,1,a)+
     &                          gammahg(1,3,b)*gammahg(3,1,a)+
     &                          gammahg(1,4,b)*gammahg(4,1,a)+
     &                          gammahg(2,1,b)*gammahg(1,2,a)+
     &                          gammahg(2,2,b)*gammahg(2,2,a)+
     &                          gammahg(2,3,b)*gammahg(3,2,a)+
     &                          gammahg(2,4,b)*gammahg(4,2,a)+
     &                          gammahg(3,1,b)*gammahg(1,3,a)+
     &                          gammahg(3,2,b)*gammahg(2,3,a)+
     &                          gammahg(3,3,b)*gammahg(3,3,a)+
     &                          gammahg(3,4,b)*gammahg(4,3,a)+
     &                          gammahg(4,1,b)*gammahg(1,4,a)+
     &                          gammahg(4,2,b)*gammahg(2,4,a)+
     &                          gammahg(4,3,b)*gammahg(3,4,a)+
     &                          gammahg(4,4,b)*gammahg(4,4,a)
     &                            )
                  term8(a,b)=-lambda4*h0_ll(a,b)
     &
                  cfe(a,b)=term1(a,b)+term2(a,b)+term3(a,b)+term4(a,b)
     &                    +term5(a,b)+term6(a,b)+term7(a,b)+term8(a,b)
     &                    -8*PI*(set_ll(a,b)-tr_set*g0_ll(a,b)/2)

                end do
              end do

              !--------------------------------------------------------------------------
              ! phi1_res = g^ab phi1,ab + g^ab,a phi1,b + g^cb gamma^a_ab phi1,c  
              !         (= g^ab phi1,ab - g^ab gamma^c_ab phi1,c) 
              !--------------------------------------------------------------------------
              phi1_res= phi10_xx(1,1)*g0_uu(1,1)+
     &                  phi10_xx(2,2)*g0_uu(2,2)+
     &                  phi10_xx(3,3)*g0_uu(3,3)+
     &                  phi10_xx(4,4)*g0_uu(4,4)+
     &               2*(phi10_xx(1,2)*g0_uu(1,2)+
     &                  phi10_xx(1,3)*g0_uu(1,3)+
     &                  phi10_xx(1,4)*g0_uu(1,4)+
     &                  phi10_xx(2,3)*g0_uu(2,3)+
     &                  phi10_xx(2,4)*g0_uu(2,4)+
     &                  phi10_xx(3,4)*g0_uu(3,4))
     &              +
     &                  phi10_x(1)*g0_uu_x(1,1,1)+
     &                  phi10_x(1)*g0_uu_x(2,1,2)+
     &                  phi10_x(1)*g0_uu_x(3,1,3)+
     &                  phi10_x(1)*g0_uu_x(4,1,4)+
     &                  phi10_x(2)*g0_uu_x(1,2,1)+
     &                  phi10_x(2)*g0_uu_x(2,2,2)+
     &                  phi10_x(2)*g0_uu_x(3,2,3)+
     &                  phi10_x(2)*g0_uu_x(4,2,4)+
     &                  phi10_x(3)*g0_uu_x(1,3,1)+
     &                  phi10_x(3)*g0_uu_x(2,3,2)+
     &                  phi10_x(3)*g0_uu_x(3,3,3)+
     &                  phi10_x(3)*g0_uu_x(4,3,4)+
     &                  phi10_x(4)*g0_uu_x(1,4,1)+
     &                  phi10_x(4)*g0_uu_x(2,4,2)+
     &                  phi10_x(4)*g0_uu_x(3,4,3)+
     &                  phi10_x(4)*g0_uu_x(4,4,4)
     &              +
     &                  phi10_x(1)*g0_uu(1,1)*gamma_ull(1,1,1)+
     &                  phi10_x(1)*g0_uu(1,1)*gamma_ull(2,2,1)+ 
     &                  phi10_x(1)*g0_uu(1,1)*gamma_ull(3,3,1)+
     &                  phi10_x(1)*g0_uu(1,1)*gamma_ull(4,4,1)+
     &                  phi10_x(1)*g0_uu(1,2)*gamma_ull(1,1,2)+
     &                  phi10_x(1)*g0_uu(1,2)*gamma_ull(2,2,2)+
     &                  phi10_x(1)*g0_uu(1,2)*gamma_ull(3,3,2)+
     &                  phi10_x(1)*g0_uu(1,2)*gamma_ull(4,4,2)+
     &                  phi10_x(1)*g0_uu(1,3)*gamma_ull(1,1,3)+
     &                  phi10_x(1)*g0_uu(1,3)*gamma_ull(2,2,3)+
     &                  phi10_x(1)*g0_uu(1,3)*gamma_ull(3,3,3)+
     &                  phi10_x(1)*g0_uu(1,3)*gamma_ull(4,4,3)+
     &                  phi10_x(1)*g0_uu(1,4)*gamma_ull(1,1,4)+
     &                  phi10_x(1)*g0_uu(1,4)*gamma_ull(2,2,4)+
     &                  phi10_x(1)*g0_uu(1,4)*gamma_ull(3,3,4)+
     &                  phi10_x(1)*g0_uu(1,4)*gamma_ull(4,4,4)+
     &                  phi10_x(2)*g0_uu(2,1)*gamma_ull(1,1,1)+
     &                  phi10_x(2)*g0_uu(2,1)*gamma_ull(2,2,1)+ 
     &                  phi10_x(2)*g0_uu(2,1)*gamma_ull(3,3,1)+
     &                  phi10_x(2)*g0_uu(2,1)*gamma_ull(4,4,1)+
     &                  phi10_x(2)*g0_uu(2,2)*gamma_ull(1,1,2)+
     &                  phi10_x(2)*g0_uu(2,2)*gamma_ull(2,2,2)+
     &                  phi10_x(2)*g0_uu(2,2)*gamma_ull(3,3,2)+
     &                  phi10_x(2)*g0_uu(2,2)*gamma_ull(4,4,2)+
     &                  phi10_x(2)*g0_uu(2,3)*gamma_ull(1,1,3)+
     &                  phi10_x(2)*g0_uu(2,3)*gamma_ull(2,2,3)+
     &                  phi10_x(2)*g0_uu(2,3)*gamma_ull(3,3,3)+
     &                  phi10_x(2)*g0_uu(2,3)*gamma_ull(4,4,3)+
     &                  phi10_x(2)*g0_uu(2,4)*gamma_ull(1,1,4)+
     &                  phi10_x(2)*g0_uu(2,4)*gamma_ull(2,2,4)+
     &                  phi10_x(2)*g0_uu(2,4)*gamma_ull(3,3,4)+
     &                  phi10_x(2)*g0_uu(2,4)*gamma_ull(4,4,4)+
     &                  phi10_x(3)*g0_uu(3,1)*gamma_ull(1,1,1)+
     &                  phi10_x(3)*g0_uu(3,1)*gamma_ull(2,2,1)+ 
     &                  phi10_x(3)*g0_uu(3,1)*gamma_ull(3,3,1)+
     &                  phi10_x(3)*g0_uu(3,1)*gamma_ull(4,4,1)+
     &                  phi10_x(3)*g0_uu(3,2)*gamma_ull(1,1,2)+
     &                  phi10_x(3)*g0_uu(3,2)*gamma_ull(2,2,2)+
     &                  phi10_x(3)*g0_uu(3,2)*gamma_ull(3,3,2)+
     &                  phi10_x(3)*g0_uu(3,2)*gamma_ull(4,4,2)+
     &                  phi10_x(3)*g0_uu(3,3)*gamma_ull(1,1,3)+
     &                  phi10_x(3)*g0_uu(3,3)*gamma_ull(2,2,3)+
     &                  phi10_x(3)*g0_uu(3,3)*gamma_ull(3,3,3)+
     &                  phi10_x(3)*g0_uu(3,3)*gamma_ull(4,4,3)+
     &                  phi10_x(3)*g0_uu(3,4)*gamma_ull(1,1,4)+
     &                  phi10_x(3)*g0_uu(3,4)*gamma_ull(2,2,4)+
     &                  phi10_x(3)*g0_uu(3,4)*gamma_ull(3,3,4)+
     &                  phi10_x(3)*g0_uu(3,4)*gamma_ull(4,4,4)+
     &                  phi10_x(4)*g0_uu(4,1)*gamma_ull(1,1,1)+
     &                  phi10_x(4)*g0_uu(4,1)*gamma_ull(2,2,1)+ 
     &                  phi10_x(4)*g0_uu(4,1)*gamma_ull(3,3,1)+
     &                  phi10_x(4)*g0_uu(4,1)*gamma_ull(4,4,1)+
     &                  phi10_x(4)*g0_uu(4,2)*gamma_ull(1,1,2)+
     &                  phi10_x(4)*g0_uu(4,2)*gamma_ull(2,2,2)+
     &                  phi10_x(4)*g0_uu(4,2)*gamma_ull(3,3,2)+
     &                  phi10_x(4)*g0_uu(4,2)*gamma_ull(4,4,2)+
     &                  phi10_x(4)*g0_uu(4,3)*gamma_ull(1,1,3)+
     &                  phi10_x(4)*g0_uu(4,3)*gamma_ull(2,2,3)+
     &                  phi10_x(4)*g0_uu(4,3)*gamma_ull(3,3,3)+
     &                  phi10_x(4)*g0_uu(4,3)*gamma_ull(4,4,3)+
     &                  phi10_x(4)*g0_uu(4,4)*gamma_ull(1,1,4)+
     &                  phi10_x(4)*g0_uu(4,4)*gamma_ull(2,2,4)+
     &                  phi10_x(4)*g0_uu(4,4)*gamma_ull(3,3,4)+
     &                  phi10_x(4)*g0_uu(4,4)*gamma_ull(4,4,4)

              !---------------------------------------------------------------- 
              ! computes diag. Jacobian of g_np1->L.g_np1 transformation
              ! by differentiating L.g wrt. g(a,b)_ij_np1 diag. entries
              ! 
              ! ddgb_J_tx differs from ddgb_J due to forward/backward stencils
              ! at excision surfaces that affect the cross-derivative tx
              ! (this is the only contribution, since the diag Jacobian is diff wrt. g_ij_np1)
              !---------------------------------------------------------------- 
              dgb_J=1d0/2d0/dt
              ddgb_J=1d0/dt/dt

              if (i.eq.1.or.(chr(i-1).eq.ex)) then
                 if (i.le.(Nx-3)
     &               .and.((chr(i+1).ne.ex
     &               .and.chr(i+2).ne.ex
     &               .and.chr(i+3).ne.ex))) then
                    ddgb_J_tx=-1d0/dt/dx
                 else if (i.le.(Nx-2)
     &                    .and.((chr(i+1).ne.ex
     &                    .and.chr(i+2).ne.ex))) then
                    ddgb_J_tx=-3d0/4d0/dt/dx
                 else if (i.le.(Nx-1).and.chr(i+1).ne.ex) then
                    ddgb_J_tx=-1d0/2d0/dt/dx
                 else
                    write(*,*) 'g_evo_opt: error in chr stencil (A)'
                    write(*,*) '    i,Nx,dx=',i,Nx,dx
                    write(*,*) '    (first error only)'
                    ddgb_J_tx=0
                 end if
              else if (i.eq.Nx.or.(chr(i+1).eq.ex)) then
                 if (i.ge.4
     &               .and.((chr(i-1).ne.ex
     &               .and.chr(i-2).ne.ex
     &               .and.chr(i-3).ne.ex))) then
                    ddgb_J_tx=1d0/dt/dx
                 else if (i.ge.3
     &                    .and.((chr(i-1).ne.ex
     &                    .and.chr(i-2).ne.ex))) then
                    ddgb_J_tx=3d0/4d0/dt/dx
                 else if (i.ge.2.and.chr(i-1).ne.ex) then
                    ddgb_J_tx=1d0/2d0/dt/dx
                 else
                    write(*,*) 'g_evo_opt: error in chr stencil (B)'
                    write(*,*) '    i,Nx,dx=',i,Nx,dx
                    write(*,*) '    (first error only)'
                    ddgb_J_tx=0
                 end if
              else
                 if ((chr(i+1).ne.ex.and.chr(i-1).ne.ex)) then
                    ddgb_J_tx=0
                 else
                    write(*,*) 'g_evo_opt: error in chr stencil (C)'
                    write(*,*) '    i,Nx,dx=',i,Nx,dx
                    write(*,*) '    (first error only)'
                    ddgb_J_tx=0
                 end if
              end if

              cfe_J(1,1)=    -0.5d0*(
     &                          g0_uu(1,1)*ddgb_J
     &                          -4*x0*g0_uu(1,2)*dgb_J
     &                          +2*g0_uu(1,2)*ddgb_J_tx
     &                              )
     &                    
     &                       -0.5d0*(
     &                          -dgb_J*
     &                          (g0_uu(1,1)*g0_uu(1,1)*g0_ll_x(1,1,1)+
     &                           g0_uu(1,1)*g0_uu(2,1)*g0_ll_x(1,1,2)+
     &                           g0_uu(1,1)*g0_uu(3,1)*g0_ll_x(1,1,3)+
     &                           g0_uu(1,1)*g0_uu(4,1)*g0_ll_x(1,1,4)+
     &                           g0_uu(2,1)*g0_uu(1,1)*g0_ll_x(1,2,1)+
     &                           g0_uu(2,1)*g0_uu(2,1)*g0_ll_x(1,2,2)+
     &                           g0_uu(2,1)*g0_uu(3,1)*g0_ll_x(1,2,3)+
     &                           g0_uu(2,1)*g0_uu(4,1)*g0_ll_x(1,2,4)+
     &                           g0_uu(3,1)*g0_uu(1,1)*g0_ll_x(1,3,1)+
     &                           g0_uu(3,1)*g0_uu(2,1)*g0_ll_x(1,3,2)+
     &                           g0_uu(3,1)*g0_uu(3,1)*g0_ll_x(1,3,3)+
     &                           g0_uu(3,1)*g0_uu(4,1)*g0_ll_x(1,3,4)+
     &                           g0_uu(4,1)*g0_uu(1,1)*g0_ll_x(1,4,1)+
     &                           g0_uu(4,1)*g0_uu(2,1)*g0_ll_x(1,4,2)+
     &                           g0_uu(4,1)*g0_uu(3,1)*g0_ll_x(1,4,3)+
     &                           g0_uu(4,1)*g0_uu(4,1)*g0_ll_x(1,4,4))
     &                          +dgb_J*
     &                          (g0_uu_x(1,1,1))
     &                              )*2
     &
     &                       +      (
     &                          0.5d0*dgb_J*
     &                          ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
     &                           (Hads_l(2)+A_l(2))*g0_uu(2,1)+
     &                           (Hads_l(3)+A_l(3))*g0_uu(3,1)+
     &                           (Hads_l(4)+A_l(4))*g0_uu(4,1))
     &                              )  
     &                    
     &                       -      (
     &                          0.25d0*dgb_J*
     &                          (cuuuu(1,1,1,1)*dlll(1,1,1)+
     &                           cuuuu(1,1,1,2)*dlll(1,1,2)+
     &                           cuuuu(1,1,1,3)*dlll(1,1,3)+
     &                           cuuuu(1,1,1,4)*dlll(1,1,4)+
     &                           cuuuu(2,1,1,1)*dlll(2,1,1)+
     &                           cuuuu(2,1,1,2)*dlll(2,1,2)+
     &                           cuuuu(2,1,1,3)*dlll(2,1,3)+
     &                           cuuuu(2,1,1,4)*dlll(2,1,4)+
     &                           cuuuu(3,1,1,1)*dlll(3,1,1)+
     &                           cuuuu(3,1,1,2)*dlll(3,1,2)+
     &                           cuuuu(3,1,1,3)*dlll(3,1,3)+
     &                           cuuuu(3,1,1,4)*dlll(3,1,4)+
     &                           cuuuu(4,1,1,1)*dlll(4,1,1)+
     &                           cuuuu(4,1,1,2)*dlll(4,1,2)+
     &                           cuuuu(4,1,1,3)*dlll(4,1,3)+
     &                           cuuuu(4,1,1,4)*dlll(4,1,4)
     &                         +
     &                           cuuuu(1,1,1,1)*dlll(1,1,1)+
     &                           cuuuu(1,1,2,1)*dlll(2,1,1)+
     &                           cuuuu(1,1,3,1)*dlll(3,1,1)+
     &                           cuuuu(1,1,4,1)*dlll(4,1,1)+
     &                           cuuuu(1,2,1,1)*dlll(1,1,2)+
     &                           cuuuu(1,2,2,1)*dlll(2,1,2)+
     &                           cuuuu(1,2,3,1)*dlll(3,1,2)+
     &                           cuuuu(1,2,4,1)*dlll(4,1,2)+
     &                           cuuuu(1,3,1,1)*dlll(1,1,3)+
     &                           cuuuu(1,3,2,1)*dlll(2,1,3)+
     &                           cuuuu(1,3,3,1)*dlll(3,1,3)+
     &                           cuuuu(1,3,4,1)*dlll(4,1,3)+
     &                           cuuuu(1,4,1,1)*dlll(1,1,4)+
     &                           cuuuu(1,4,2,1)*dlll(2,1,4)+
     &                           cuuuu(1,4,3,1)*dlll(3,1,4)+
     &                           cuuuu(1,4,4,1)*dlll(4,1,4))
     &                              )

              cfe_J(1,2)=  -0.5d0*(
     &                        g0_uu(1,1)*ddgb_J
     &                        -4*(1)*x0*g0_uu(1,2)*dgb_J
     &                        +2*g0_uu(1,2)*ddgb_J_tx
     &                            )
     &                     
     &                     -0.5d0*(
     &                        -dgb_J*
     &                        (g0_uu(1,1)*g0_uu(1,2)*g0_ll_x(2,1,1)+
     &                         g0_uu(1,1)*g0_uu(2,2)*g0_ll_x(2,1,2)+
     &                         g0_uu(1,1)*g0_uu(3,2)*g0_ll_x(2,1,3)+
     &                         g0_uu(1,1)*g0_uu(4,2)*g0_ll_x(2,1,4)+
     &                         g0_uu(2,1)*g0_uu(1,2)*g0_ll_x(2,2,1)+
     &                         g0_uu(2,1)*g0_uu(2,2)*g0_ll_x(2,2,2)+
     &                         g0_uu(2,1)*g0_uu(3,2)*g0_ll_x(2,2,3)+
     &                         g0_uu(2,1)*g0_uu(4,2)*g0_ll_x(2,2,4)+
     &                         g0_uu(3,1)*g0_uu(1,2)*g0_ll_x(2,3,1)+
     &                         g0_uu(3,1)*g0_uu(2,2)*g0_ll_x(2,3,2)+
     &                         g0_uu(3,1)*g0_uu(3,2)*g0_ll_x(2,3,3)+
     &                         g0_uu(3,1)*g0_uu(4,2)*g0_ll_x(2,3,4)+
     &                         g0_uu(4,1)*g0_uu(1,2)*g0_ll_x(2,4,1)+
     &                         g0_uu(4,1)*g0_uu(2,2)*g0_ll_x(2,4,2)+
     &                         g0_uu(4,1)*g0_uu(3,2)*g0_ll_x(2,4,3)+
     &                         g0_uu(4,1)*g0_uu(4,2)*g0_ll_x(2,4,4))
     &                        +dgb_J*
     &                        (g0_uu_x(1,1,1))
     &                            )
     &                  
     &                     -0.5d0*(
     &                        dgb_J*
     &                        (g0_uu_x(2,1,2))
     &                            )
     &
     &                     -      (
     &                        0.5d0*dgb_J*
     &                        (cuuuu(1,1,1,2)*dlll(1,2,1)+
     &                         cuuuu(1,1,2,2)*dlll(2,2,1)+
     &                         cuuuu(1,1,3,2)*dlll(3,2,1)+
     &                         cuuuu(1,1,4,2)*dlll(4,2,1)+
     &                         cuuuu(1,2,1,2)*dlll(1,2,2)+
     &                         cuuuu(1,2,2,2)*dlll(2,2,2)+
     &                         cuuuu(1,2,3,2)*dlll(3,2,2)+
     &                         cuuuu(1,2,4,2)*dlll(4,2,2)+
     &                         cuuuu(1,3,1,2)*dlll(1,2,3)+
     &                         cuuuu(1,3,2,2)*dlll(2,2,3)+
     &                         cuuuu(1,3,3,2)*dlll(3,2,3)+
     &                         cuuuu(1,3,4,2)*dlll(4,2,3)+
     &                         cuuuu(1,4,1,2)*dlll(1,2,4)+
     &                         cuuuu(1,4,2,2)*dlll(2,2,4)+
     &                         cuuuu(1,4,3,2)*dlll(3,2,4)+
     &                         cuuuu(1,4,4,2)*dlll(4,2,4))
     &                            )

              cfe_J(2,2)=-0.5d0*(
     &                      g0_uu(1,1)*ddgb_J
     &                      -4*x0*g0_uu(1,2)*dgb_J
     &                      +2*g0_uu(1,2)*ddgb_J_tx
     &                          )
     &                
     &                   -0.5d0*(
     &                      +dgb_J*
     &                      (g0_uu_x(2,1,2))
     &                          )
     &                
     &                   -0.5d0*(
     &                      +dgb_J*
     &                      (g0_uu_x(2,1,2))
     &                          )
     &
     &                   +      (
     &                      -0.5d0*dgb_J*
     &                      ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
     &                       (Hads_l(2)+A_l(2))*g0_uu(2,1)+
     &                       (Hads_l(3)+A_l(3))*g0_uu(3,1)+
     &                       (Hads_l(4)+A_l(4))*g0_uu(4,1))
     &                          )  
     &                
     &                   -      (
     &                      0.25d0*dgb_J*
     &                      (cuuuu(1,2,1,1)*dlll(1,2,1)+
     &                       cuuuu(1,2,1,2)*dlll(1,2,2)+
     &                       cuuuu(1,2,1,3)*dlll(1,2,3)+
     &                       cuuuu(1,2,1,4)*dlll(1,2,4)+
     &                       cuuuu(2,2,1,1)*dlll(2,2,1)+
     &                       cuuuu(2,2,1,2)*dlll(2,2,2)+
     &                       cuuuu(2,2,1,3)*dlll(2,2,3)+
     &                       cuuuu(2,2,1,4)*dlll(2,2,4)+
     &                       cuuuu(3,2,1,1)*dlll(3,2,1)+
     &                       cuuuu(3,2,1,2)*dlll(3,2,2)+
     &                       cuuuu(3,2,1,3)*dlll(3,2,3)+
     &                       cuuuu(3,2,1,4)*dlll(3,2,4)+
     &                       cuuuu(4,2,1,1)*dlll(4,2,1)+
     &                       cuuuu(4,2,1,2)*dlll(4,2,2)+
     &                       cuuuu(4,2,1,3)*dlll(4,2,3)+
     &                       cuuuu(4,2,1,4)*dlll(4,2,4)+
     &                       cuuuu(1,1,2,1)*dlll(1,2,1)-
     &                       cuuuu(1,1,2,2)*dlll(1,2,2)-
     &                       cuuuu(1,1,2,3)*dlll(1,2,3)-
     &                       cuuuu(1,1,2,4)*dlll(1,2,4)-
     &                       cuuuu(2,1,2,1)*dlll(2,2,1)-
     &                       cuuuu(2,1,2,2)*dlll(2,2,2)-
     &                       cuuuu(2,1,2,3)*dlll(2,2,3)-
     &                       cuuuu(2,1,2,4)*dlll(2,2,4)-
     &                       cuuuu(3,1,2,1)*dlll(3,2,1)-
     &                       cuuuu(3,1,2,2)*dlll(3,2,2)-
     &                       cuuuu(3,1,2,3)*dlll(3,2,3)-
     &                       cuuuu(3,1,2,4)*dlll(3,2,4)-
     &                       cuuuu(4,1,2,1)*dlll(4,2,1)-
     &                       cuuuu(4,1,2,2)*dlll(4,2,2)-
     &                       cuuuu(4,1,2,3)*dlll(4,2,3)-
     &                       cuuuu(4,1,2,4)*dlll(4,2,4)
     &                     +
     &                       cuuuu(1,1,1,2)*dlll(1,2,1)+
     &                       cuuuu(1,1,2,2)*dlll(2,2,1)+
     &                       cuuuu(1,1,3,2)*dlll(3,2,1)+
     &                       cuuuu(1,1,4,2)*dlll(4,2,1)+
     &                       cuuuu(1,2,1,2)*dlll(1,2,2)+
     &                       cuuuu(1,2,2,2)*dlll(2,2,2)+
     &                       cuuuu(1,2,3,2)*dlll(3,2,2)+
     &                       cuuuu(1,2,4,2)*dlll(4,2,2)+
     &                       cuuuu(1,3,1,2)*dlll(1,2,3)+
     &                       cuuuu(1,3,2,2)*dlll(2,2,3)+
     &                       cuuuu(1,3,3,2)*dlll(3,2,3)+
     &                       cuuuu(1,3,4,2)*dlll(4,2,3)+
     &                       cuuuu(1,4,1,2)*dlll(1,2,4)+
     &                       cuuuu(1,4,2,2)*dlll(2,2,4)+
     &                       cuuuu(1,4,3,2)*dlll(3,2,4)+
     &                       cuuuu(1,4,4,2)*dlll(4,2,4)+
     &                       cuuuu(2,1,1,1)*dlll(1,2,1)-
     &                       cuuuu(2,1,2,1)*dlll(2,2,1)-
     &                       cuuuu(2,1,3,1)*dlll(3,2,1)-
     &                       cuuuu(2,1,4,1)*dlll(4,2,1)-
     &                       cuuuu(2,2,1,1)*dlll(1,2,2)-
     &                       cuuuu(2,2,2,1)*dlll(2,2,2)-
     &                       cuuuu(2,2,3,1)*dlll(3,2,2)-
     &                       cuuuu(2,2,4,1)*dlll(4,2,2)-
     &                       cuuuu(2,3,1,1)*dlll(1,2,3)-
     &                       cuuuu(2,3,2,1)*dlll(2,2,3)-
     &                       cuuuu(2,3,3,1)*dlll(3,2,3)-
     &                       cuuuu(2,3,4,1)*dlll(4,2,3)-
     &                       cuuuu(2,4,1,1)*dlll(1,2,4)-
     &                       cuuuu(2,4,2,1)*dlll(2,2,4)-
     &                       cuuuu(2,4,3,1)*dlll(3,2,4)-
     &                       cuuuu(2,4,4,1)*dlll(4,2,4))
     &                          )

              cfe_J(3,3)=    -0.5d0*(
     &                          x0**2*g0_uu(1,1)*ddgb_J
     &                          +2*(2*x0-4*x0**3)*g0_uu(1,2)*dgb_J
     &                          +2*x0**2*g0_uu(1,2)*ddgb_J_tx
     &                              )
     &                    
     &                       -0.5d0*(
     &                          +x0**2*dgb_J
     &                                      *(g0_uu_x(3,1,3))
     &                              )
     &                    
     &                       -0.5d0*(
     &                          +x0**2*dgb_J
     &                                      *(g0_uu_x(3,1,3))
     &                              )
     &
     &                       +      (
     &                          -0.5d0*x0**2*dgb_J*
     &                          ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
     &                           (Hads_l(2)+A_l(2))*g0_uu(2,1)+
     &                           (Hads_l(3)+A_l(3))*g0_uu(3,1)+
     &                           (Hads_l(4)+A_l(4))*g0_uu(4,1))
     &                              )  
     &                    
     &                       -      (
     &                          0.25d0*x0**2*dgb_J*
     &                          (cuuuu(1,3,1,1)*dlll(1,3,1)+
     &                           cuuuu(1,3,1,2)*dlll(1,3,2)+
     &                           cuuuu(1,3,1,3)*dlll(1,3,3)+
     &                           cuuuu(1,3,1,4)*dlll(1,3,4)+
     &                           cuuuu(2,3,1,1)*dlll(2,3,1)+
     &                           cuuuu(2,3,1,2)*dlll(2,3,2)+
     &                           cuuuu(2,3,1,3)*dlll(2,3,3)+
     &                           cuuuu(2,3,1,4)*dlll(2,3,4)+
     &                           cuuuu(3,3,1,1)*dlll(3,3,1)+
     &                           cuuuu(3,3,1,2)*dlll(3,3,2)+
     &                           cuuuu(3,3,1,3)*dlll(3,3,3)+
     &                           cuuuu(3,3,1,4)*dlll(3,3,4)+
     &                           cuuuu(4,3,1,1)*dlll(4,3,1)+
     &                           cuuuu(4,3,1,2)*dlll(4,3,2)+
     &                           cuuuu(4,3,1,3)*dlll(4,3,3)+
     &                           cuuuu(4,3,1,4)*dlll(4,3,4)+
     &                           cuuuu(1,1,3,1)*dlll(1,3,1)-
     &                           cuuuu(1,1,3,2)*dlll(1,3,2)-
     &                           cuuuu(1,1,3,3)*dlll(1,3,3)-
     &                           cuuuu(1,1,3,4)*dlll(1,3,4)-
     &                           cuuuu(2,1,3,1)*dlll(2,3,1)-
     &                           cuuuu(2,1,3,2)*dlll(2,3,2)-
     &                           cuuuu(2,1,3,3)*dlll(2,3,3)-
     &                           cuuuu(2,1,3,4)*dlll(2,3,4)-
     &                           cuuuu(3,1,3,1)*dlll(3,3,1)-
     &                           cuuuu(3,1,3,2)*dlll(3,3,2)-
     &                           cuuuu(3,1,3,3)*dlll(3,3,3)-
     &                           cuuuu(3,1,3,4)*dlll(3,3,4)-
     &                           cuuuu(4,1,3,1)*dlll(4,3,1)-
     &                           cuuuu(4,1,3,2)*dlll(4,3,2)-
     &                           cuuuu(4,1,3,3)*dlll(4,3,3)-
     &                           cuuuu(4,1,3,4)*dlll(4,3,4)
     &                         +
     &                           cuuuu(1,1,1,3)*dlll(1,3,1)+
     &                           cuuuu(1,1,2,3)*dlll(2,3,1)+
     &                           cuuuu(1,1,3,3)*dlll(3,3,1)+
     &                           cuuuu(1,1,4,3)*dlll(4,3,1)+
     &                           cuuuu(1,2,1,3)*dlll(1,3,2)+
     &                           cuuuu(1,2,2,3)*dlll(2,3,2)+
     &                           cuuuu(1,2,3,3)*dlll(3,3,2)+
     &                           cuuuu(1,2,4,3)*dlll(4,3,2)+
     &                           cuuuu(1,3,1,3)*dlll(1,3,3)+
     &                           cuuuu(1,3,2,3)*dlll(2,3,3)+
     &                           cuuuu(1,3,3,3)*dlll(3,3,3)+
     &                           cuuuu(1,3,4,3)*dlll(4,3,3)+
     &                           cuuuu(1,4,1,3)*dlll(1,3,4)+
     &                           cuuuu(1,4,2,3)*dlll(2,3,4)+
     &                           cuuuu(1,4,3,3)*dlll(3,3,4)+
     &                           cuuuu(1,4,4,3)*dlll(4,3,4)+
     &                           cuuuu(3,1,1,1)*dlll(1,3,1)-
     &                           cuuuu(3,1,2,1)*dlll(2,3,1)-
     &                           cuuuu(3,1,3,1)*dlll(3,3,1)-
     &                           cuuuu(3,1,4,1)*dlll(4,3,1)-
     &                           cuuuu(3,2,1,1)*dlll(1,3,2)-
     &                           cuuuu(3,2,2,1)*dlll(2,3,2)-
     &                           cuuuu(3,2,3,1)*dlll(3,3,2)-
     &                           cuuuu(3,2,4,1)*dlll(4,3,2)-
     &                           cuuuu(3,3,1,1)*dlll(1,3,3)-
     &                           cuuuu(3,3,2,1)*dlll(2,3,3)-
     &                           cuuuu(3,3,3,1)*dlll(3,3,3)-
     &                           cuuuu(3,3,4,1)*dlll(4,3,3)-
     &                           cuuuu(3,4,1,1)*dlll(1,3,4)-
     &                           cuuuu(3,4,2,1)*dlll(2,3,4)-
     &                           cuuuu(3,4,3,1)*dlll(3,3,4)-
     &                           cuuuu(3,4,4,1)*dlll(4,3,4))
     &                              )


              !----------------------------------------------------------------
              ! computes diag. Jacobian of phi1_np1->L.phi1_np1 transformation
              ! by differentiating L.phi1=box.phi1-dV/dphi1 wrt. phi1_np1
              ! and remember: phi10=phi1*(1-x0**2)**2 
              ! 
              ! ddphi1_J_tx differs from ddphi_J due to forward/backward stencils
              ! at excision surfaces that affect the cross-derivative tx 
              ! (these are the only contributions, since the diag Jacobian is diff wrt. phi1_np1)
              !----------------------------------------------------------------
              dphi1_J=1d0/2d0/dt
              ddphi1_J=1d0/dt/dt

              if (i.eq.1.or.(chr(i-1).eq.ex)) then
                 if (i.le.(Nx-3)
     &               .and.((chr(i+1).ne.ex
     &               .and.chr(i+2).ne.ex
     &               .and.chr(i+3).ne.ex))) then
                    ddphi1_J_tx=-1d0/dt/dx
                 else if (i.le.(Nx-2)
     &                    .and.((chr(i+1).ne.ex
     &                    .and.chr(i+2).ne.ex))) then
                    ddphi1_J_tx=-3d0/4d0/dt/dx
                 else if (i.le.(Nx-1).and.chr(i+1).ne.ex) then
                    ddphi1_J_tx=-1d0/2d0/dt/dx
                 else
                    write(*,*) 'g_evo_opt: error in chr stencil (A)'
                    write(*,*) '    i,Nx,dx=',i,Nx,dx
                    write(*,*) '    (first error only)'
                    ddphi1_J_tx=0
                 end if
              else if (i.eq.Nx.or.(chr(i+1).eq.ex)) then
                 if (i.ge.4
     &               .and.((chr(i-1).ne.ex
     &               .and.chr(i-2).ne.ex
     &               .and.chr(i-3).ne.ex))) then
                    ddphi1_J_tx=1d0/dt/dx
                 else if (i.ge.3
     &                    .and.((chr(i-1).ne.ex
     &                    .and.chr(i-2).ne.ex))) then
                    ddphi1_J_tx=3d0/4d0/dt/dx
                 else if (i.ge.2.and.chr(i-1).ne.ex) then
                    ddphi1_J_tx=1d0/2d0/dt/dx
                 else
                    write(*,*) 'g_evo_opt: error in chr stencil (B)'
                    write(*,*) '    i,Nx,dx=',i,Nx,dx
                    write(*,*) '    (first error only)'
                    ddphi1_J_tx=0
                 end if
              else
                 if ((chr(i+1).ne.ex.and.chr(i-1).ne.ex)) then
                    ddphi1_J_tx=0
                 else
                    write(*,*) 'g_evo_opt: error in chr stencil (C)'
                    write(*,*) '    i,Nx,dx=',i,Nx,dx
                    write(*,*) '    (first error only)'
                    ddphi1_J_tx=0
                 end if
              end if

              phi1_J=    (
     &               g0_uu(1,1)*ddphi1_J
     &                         *(1-x0**2)**2
     &               -4*(2)*x0*g0_uu(1,2)*dphi1_J
     &                         *(1-x0**2)**(2)
     &               +2*g0_uu(1,2)*ddphi1_J_tx
     &                         *(1-x0**2)**(2)
     &                   )
     &              +
     &                dphi1_J*(1-x0**2)**2
     &                *(
     &                  g0_uu_x(1,1,1)+
     &                  g0_uu_x(2,1,2)+
     &                  g0_uu_x(3,1,3)+
     &                  g0_uu_x(4,1,4)
     &                 )
     &              +
     &                dphi1_J*(1-x0**2)**2
     &                *(
     &                  g0_uu(1,1)*gamma_ull(1,1,1)+
     &                  g0_uu(1,1)*gamma_ull(2,2,1)+
     &                  g0_uu(1,1)*gamma_ull(3,3,1)+
     &                  g0_uu(1,1)*gamma_ull(4,4,1)+
     &                  g0_uu(1,2)*gamma_ull(1,1,2)+
     &                  g0_uu(1,2)*gamma_ull(2,2,2)+
     &                  g0_uu(1,2)*gamma_ull(3,3,2)+
     &                  g0_uu(1,2)*gamma_ull(4,4,2)+
     &                  g0_uu(1,3)*gamma_ull(1,1,3)+
     &                  g0_uu(1,3)*gamma_ull(2,2,3)+
     &                  g0_uu(1,3)*gamma_ull(3,3,3)+
     &                  g0_uu(1,3)*gamma_ull(4,4,3)+
     &                  g0_uu(1,4)*gamma_ull(1,1,4)+
     &                  g0_uu(1,4)*gamma_ull(2,2,4)+
     &                  g0_uu(1,4)*gamma_ull(3,3,4)+
     &                  g0_uu(1,4)*gamma_ull(4,4,4)
     &                 )

              ! constraint damping terms added to cfe,cfe_J
              do a=1,4
                do b=1,4
                  cd_ll(a,b)=-kappa_cd*
     &                ( n_l(a)*c_l(b)+n_l(b)*c_l(a)
     &                  -(1+rho_cd)*g0_ll(a,b)*ndotc )
                end do
              end do

              dc_J=1d0/2d0/dt

              cd_J_ll(1,1)=-kappa_cd*
     &            ( n_l(1)*g0_uu(1,1)*dc_J
     &              -(1+rho_cd)*g0_ll(1,1)*n_u(1)*0.5d0*
     &                          g0_uu(1,1)*dc_J )
              cd_J_ll(1,2)=-kappa_cd*
     &            ( n_l(1)*g0_uu(1,1)*dc_J
     &              -(1+rho_cd)*g0_ll(1,2)*n_u(2)*g0_uu(1,1)*
     &                          dc_J )
              cd_J_ll(2,2)=-kappa_cd*
     &            ( 2*n_l(2)*g0_uu(1,2)*dc_J
     &              -(1+rho_cd)*g0_ll(2,2)*
     &              (-n_u(1)*0.5d0*g0_uu(2,2)*dc_J
     &               +n_u(2)*g0_uu(1,2)*dc_J) ) 
              cd_J_ll(3,3)=-kappa_cd*
     &            (1+rho_cd)*( n_u(1)*0.5d0*g0_uu(3,3)*
     *                         dc_J*x0**2 )
 
              if (kappa_cd.ne.0) then
                cfe(1,1)=cfe(1,1)+cd_ll(1,1)
                cfe(1,2)=cfe(1,2)+cd_ll(1,2)
                cfe(2,2)=cfe(2,2)+cd_ll(2,2)
                cfe(3,3)=cfe(3,3)+cd_ll(3,3)
                cfe_J(1,1)=cfe_J(1,1)+cd_J_ll(1,1)
                cfe_J(1,2)=cfe_J(1,2)+cd_J_ll(1,2)
                cfe_J(2,2)=cfe_J(2,2)+cd_J_ll(2,2)
                cfe_J(3,3)=cfe_J(3,3)+cd_J_ll(3,3)
              end if

              ! update gbars 
              if (is_nan(cfe(1,1)).or.is_nan(cfe_J(1,1)).or.
     &          cfe_J(1,1).eq.0) then
                dump=.true.
              else
                gb_tt_np1(i)=gb_tt_np1(i)-cfe(1,1)/cfe_J(1,1)
              end if

              if (is_nan(cfe(1,2)).or.is_nan(cfe_J(1,2)).or.
     &          cfe_J(1,2).eq.0) then
                dump=.true.
              else
                gb_tx_np1(i)=gb_tx_np1(i)-cfe(1,2)/cfe_J(1,2)
              end if

              if (is_nan(cfe(2,2)).or.is_nan(cfe_J(2,2)).or.
     &          cfe_J(2,2).eq.0) then
                dump=.true.
              else
                gb_xx_np1(i)=gb_xx_np1(i)-cfe(2,2)/cfe_J(2,2)
              end if

              if (is_nan(cfe(3,3)).or.is_nan(cfe_J(3,3)).or.
     &          cfe_J(3,3).eq.0) then
                dump=.true.
              else
                psi_np1(i)=psi_np1(i)-cfe(3,3)/cfe_J(3,3)
              end if

              ! update phi1 
              if (is_nan(phi1_res).or.is_nan(phi1_J)) then
                dump=.true.
              else
                phi1_np1(i)=phi1_np1(i)-phi1_res/phi1_J 
              end if

              ! save residuals
              gb_res(i) =
     &          max(abs(cfe(1,1)/cfe_J(1,1)),
     &              abs(cfe(1,2)/cfe_J(1,2)),
     &              abs(cfe(2,2)/cfe_J(2,2)),
     &              abs(cfe(3,3)/cfe_J(3,3)))
              kg_res(i)=abs(phi1_res/phi1_J)

              ! save pointwise max of constraint violation
              cl_res(i)=
     &          max(abs(c_l(1)),abs(c_l(2)),abs(c_l(3)),abs(c_l(4)))

              ! check for NaNs
              if (dump.and.first_nan.or.ltrace) then
                first_nan=.false.
                write(*,*)
                write(*,*) 'g_evo_opt: Nan/zero at i,Nx=',
     &                                            i,Nx
                write(*,*) 'x=',x(i)
                write(*,*) 'dt,dx=',dt,dx
                write(*,*) 'x0',x0

                write(*,*) ' at tn:'
                write(*,*) ' gb_tt np1,n,nm1:',gb_tt_np1(i),
     &                 gb_tt_n(i),gb_tt_nm1(i)
                write(*,*) ' gb_tx np1,n,nm1:',gb_tx_np1(i),
     &                 gb_tx_n(i),gb_tx_nm1(i)
                write(*,*) ' gb_xx np1,n,nm1:',gb_xx_np1(i),
     &                 gb_xx_n(i),gb_xx_nm1(i)
                write(*,*) ' psi np1,n,nm1:',psi_np1(i),
     &                 psi_n(i),psi_nm1(i)
                write(*,*) ' gads_tt :',gads_ll(1,1)
                write(*,*) ' gads_tx :',gads_ll(1,2)
                write(*,*) ' gads_xx :',gads_ll(2,2)
                write(*,*) ' gads_psi:',gads_ll(3,3)
                write(*,*) ' h0_tt :',h0_ll(1,1)
                write(*,*) ' h0_tx :',h0_ll(1,2)
                write(*,*) ' h0_xx :',h0_ll(2,2)
                write(*,*) ' h0_psi:',h0_ll(3,3)
                write(*,*) ' g0_tt :',g0_ll(1,1)
                write(*,*) ' g0_tx :',g0_ll(1,2)
                write(*,*) ' g0_xx :',g0_ll(2,2)
                write(*,*) ' g0_psi:',g0_ll(3,3)
                write(*,*) ' g0u_tt :',g0_uu(1,1)
                write(*,*) ' g0u_tx :',g0_uu(1,2)
                write(*,*) ' g0u_xx :',g0_uu(2,2)
                write(*,*) ' g0u_psi:',g0_uu(3,3)
                write(*,*) ' cd_tt:',cd_ll(1,1)
                write(*,*) ' cd_tx:',cd_ll(1,2)
                write(*,*) ' cd_xx:',cd_ll(2,2)
                write(*,*) ' cd_psi:',cd_ll(3,3)
                write(*,*) ' phi np1,n,nm1:',phi1_np1(i),
     &                   phi1_n(i),phi1_nm1(i)
                write(*,*) ' phi1_t/x:',phi1_t,phi1_x
                write(*,*) ' phi1_ti..:',phi1_tt,phi1_tx,phi1_xx

                write(*,*) ' res J:'
                write(*,*) ' tt:',cfe(1,1),cfe_J(1,1)
                write(*,*) ' tx:',cfe(1,2),cfe_J(1,2)
                write(*,*) ' xx:',cfe(2,2),cfe_J(2,2)
                write(*,*) ' psi:',cfe(3,3),cfe_J(3,3)
                write(*,*) ' phi1:',phi1_res,phi1_J
              end if

            ! (REGION) non-interior points; set to zero prior to applying bcs 
            else 
              gb_tt_np1(i) = 0
              gb_tx_np1(i) = 0
              gb_xx_np1(i) = 0
              psi_np1(i)   = 0
              phi1_np1(i)  = 0 

            endif ! (near start of main loop)

          end do

        end do

        ! (REGION) x=0; impose regularity conditions 
        call axi_reg_Hb(Hb_t_np1,Hb_x_np1,chr,ex,L,x,Nx)
        call axi_reg_g(gb_tt_np1,gb_tx_np1,gb_xx_np1,
     &                            psi_np1,chr,ex,L,x,Nx)
        call axi_reg_phi(phi1_np1,chr,ex,L,x,Nx)

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
