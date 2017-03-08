c----------------------------------------------------------------------
c in polar coordinates t,x for x in [0,1]
c
c routine for computing independent residuals 
c----------------------------------------------------------------------
        subroutine ires(efe_all_ires,phj_ires,
     &                  pij_ires,alphaq_ires,thetap_ires,
     &                  db_txtx_np1,db_txtx_n,db_txtx_nm1,
     &                  db_txty_np1,db_txty_n,db_txty_nm1,
     &                  db_txtz_np1,db_txtz_n,db_txtz_nm1,
     &                  db_txxy_np1,db_txxy_n,db_txxy_nm1,
     &                  db_txxz_np1,db_txxz_n,db_txxz_nm1,
     &                  db_txyz_np1,db_txyz_n,db_txyz_nm1,
     &                  db_tyty_np1,db_tyty_n,db_tyty_nm1,
     &                  db_tytz_np1,db_tytz_n,db_tytz_nm1,
     &                  db_tyxy_np1,db_tyxy_n,db_tyxy_nm1,
     &                  db_tyxz_np1,db_tyxz_n,db_tyxz_nm1,
     &                  db_tyyz_np1,db_tyyz_n,db_tyyz_nm1,
     &                  db_tztz_np1,db_tztz_n,db_tztz_nm1,
     &                  db_tzxy_np1,db_tzxy_n,db_tzxy_nm1,
     &                  db_tzxz_np1,db_tzxz_n,db_tzxz_nm1,
     &                  db_tzyz_np1,db_tzyz_n,db_tzyz_nm1,
     &                  db_xyxy_np1,db_xyxy_n,db_xyxy_nm1,
     &                  db_xyxz_np1,db_xyxz_n,db_xyxz_nm1,
     &                  db_xyyz_np1,db_xyyz_n,db_xyyz_nm1,
     &                  db_xzxz_np1,db_xzxz_n,db_xzxz_nm1,
     &                  db_xzyz_np1,db_xzyz_n,db_xzyz_nm1,
     &                  db_yzyz_np1,db_yzyz_n,db_yzyz_nm1,
     &                  gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                  gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                  gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                  psi_np1,psi_n,psi_nm1,
     &                  Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                  Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                  phi1_np1,phi1_n,phi1_nm1,
     &                  x,dt,chr,L,ex,Nx,phys_bdy,ghost_width)
        implicit none
        integer Nx
        integer phys_bdy(2),ghost_width(2)
        real*8 efe_all_ires(Nx),phj_ires(Nx)
        real*8 pij_ires(Nx),alphaq_ires(Nx),thetap_ires(Nx)
        real*8 chr(Nx),ex
        real*8 x(Nx),dt,L
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

        real*8 lambda4

        integer is,ie
        integer a,b,c,d,e,f,g,h

        integer i

        real*8 x0
        real*8 dx

        real*8 efe_ires(4,4)

        real*8 boxx_u(4),boxx_l(4)
        real*8 c_l(4)

        real*8 n_l(4),s_l(4)
        real*8 n_u(4),s_u(4)
        real*8 n_l_x(4,4),s_l_x(4,4)
        real*8 sigma_uu(4,4)
        real*8 thetap_ads
        real*8 thetap_adsbhunitr0

        logical is_nan

        real*8 PI
        parameter (PI=3.141592653589793d0)

        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 phi1_t,phi1_x,phi1_y,phi1_z
        real*8 phi1_tt,phi1_tx,phi1_ty,phi1_tz
        real*8 phi1_xx,phi1_xy,phi1_xz
        real*8 phi1_yy,phi1_yz
        real*8 phi1_zz

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,chi,theta)
        !--------------------------------------------------------------
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
        ! initialize fixed-size variables 
        !--------------------------------------------------------------
        data is,ie/0,0/
        data a,b,c,d,e,f,g,h/0,0,0,0,0,0,0,0/

        data dx/0.0/

        data boxx_u,boxx_l/4*0.0,4*0.0/
        data c_l/4*0.0/

        data n_l,s_l/4*0.0,4*0.0/
        data n_u,s_u/4*0.0,4*0.0/
        data n_l_x,s_l_x/16*0.0,16*0.0/
        data sigma_uu/16*0.0/

        data phi1_t,phi1_x,phi1_y/0.0,0.0,0.0/
        data phi1_tt,phi1_tx,phi1_ty,phi1_tz/0.0,0.0,0.0,0.0/
        data phi1_xx,phi1_xy,phi1_xz/0.0,0.0,0.0/
        data phi1_yy,phi1_yz/0.0,0.0/
        data phi1_zz/0.0/

        data g0_ll,g0_uu/16*0.0,16*0.0/
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

!----------------------------------------------------------------------
        
        dx=(x(2)-x(1))

        ! AdS4D cosmological constant
        !(lambda4=-(n-1)(n-2)/2/L^2) for n=4 dimensional AdS)
        lambda4=-3/L/L

        ! set index bounds for main loop
        is=2
        ie=Nx-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-2
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-2)

        ! (MAIN LOOP) loop through spacetime points x(i)
        do i=is,ie
          x0=x(i)

          if (chr(i).ne.ex) then

            ! computes tensors at point i
            call tensor_init(
     &              db_txtx_np1,db_txtx_n,db_txtx_nm1,
     &              db_txty_np1,db_txty_n,db_txty_nm1,
     &              db_txtz_np1,db_txtz_n,db_txtz_nm1,
     &              db_txxy_np1,db_txxy_n,db_txxy_nm1,
     &              db_txxz_np1,db_txxz_n,db_txxz_nm1,
     &              db_txyz_np1,db_txyz_n,db_txyz_nm1,
     &              db_tyty_np1,db_tyty_n,db_tyty_nm1,
     &              db_tytz_np1,db_tytz_n,db_tytz_nm1,
     &              db_tyxy_np1,db_tyxy_n,db_tyxy_nm1,
     &              db_tyxz_np1,db_tyxz_n,db_tyxz_nm1,
     &              db_tyyz_np1,db_tyyz_n,db_tyyz_nm1,
     &              db_tztz_np1,db_tztz_n,db_tztz_nm1,
     &              db_tzxy_np1,db_tzxy_n,db_tzxy_nm1,
     &              db_tzxz_np1,db_tzxz_n,db_tzxz_nm1,
     &              db_tzyz_np1,db_tzyz_n,db_tzyz_nm1,
     &              db_xyxy_np1,db_xyxy_n,db_xyxy_nm1,
     &              db_xyxz_np1,db_xyxz_n,db_xyxz_nm1,
     &              db_xyyz_np1,db_xyyz_n,db_xyyz_nm1,
     &              db_xzxz_np1,db_xzxz_n,db_xzxz_nm1,
     &              db_xzyz_np1,db_xzyz_n,db_xzyz_nm1,
     &              db_yzyz_np1,db_yzyz_n,db_yzyz_nm1,
     &              gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &              gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &              gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &              psi_np1,psi_n,psi_nm1,
     &              Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &              Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &              phi1_np1,phi1_n,phi1_nm1,
     &              g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &              gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &              h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &              A_l,A_l_x,Hads_l,
     &              gamma_ull,gamma_ull_x,
     &              riemann_ulll,ricci_ll,ricci_lu,ricci,
     &              einstein_ll,set_ll,
     &              phi10_x,phi10_xx,
     &              x,dt,chr,L,ex,Nx,i)

              ! calculates efe_ires functions at point i,j
              !(efe_ires_ab=G_ab+lambda4*g_ab-8*PI*T_ab)
              do a=1,4
                do b=a,4
                  efe_ires(a,b)=einstein_ll(a,b)+lambda4*g0_ll(a,b)
     &                                          -8*PI*set_ll(a,b)
                end do
              end do

              ! calculate efe_all_ires function at point i
              efe_all_ires(i)=
     &        max(abs(efe_ires(1,1)),abs(efe_ires(1,2)),
     &            abs(efe_ires(2,2)),abs(efe_ires(3,3)))

              ! auxiliary objects for phj_ires
              do c=1,4
                boxx_u(c)=-( gamma_ull(c,1,1)*g0_uu(1,1)+
     &                       gamma_ull(c,2,2)*g0_uu(2,2)+
     &                       gamma_ull(c,3,3)*g0_uu(3,3)+
     &                       gamma_ull(c,4,4)*g0_uu(4,4)+
     &                    2*(gamma_ull(c,1,2)*g0_uu(1,2)+
     &                       gamma_ull(c,1,3)*g0_uu(1,3)+
     &                       gamma_ull(c,1,4)*g0_uu(1,4)+
     &                       gamma_ull(c,2,3)*g0_uu(2,3)+
     &                       gamma_ull(c,2,4)*g0_uu(2,4)+
     &                       gamma_ull(c,3,4)*g0_uu(3,4)) )
              end do

              do a=1,4
                boxx_l(a)=boxx_u(1)*g0_ll(a,1)+
     &                    boxx_u(2)*g0_ll(a,2)+
     &                    boxx_u(3)*g0_ll(a,3)+
     &                    boxx_u(4)*g0_ll(a,4)
              end do

              do a=1,4
                c_l(a)=(Hads_l(a)+A_l(a))-boxx_l(a)
              end do

              ! calculate alphaq_ires function at point i
              alphaq_ires(i)=1/sqrt(-g0_uu(1,1))*(1-x0)
              
              ! auxiliary objects for thetap_ires
              n_l(1)=-1/sqrt(-g0_uu(1,1))
              s_u(2)=1/sqrt(g0_ll(2,2))
              do a=1,4
                n_u(a)=n_l(1)*g0_uu(a,1)+
     &                 n_l(2)*g0_uu(a,2)+
     &                 n_l(3)*g0_uu(a,3)+
     &                 n_l(4)*g0_uu(a,4)
              end do
              do a=1,4
                s_l(a)=s_u(1)*g0_ll(a,1)+
     &                 s_u(2)*g0_ll(a,2)+
     &                 s_u(3)*g0_ll(a,3)+
     &                 s_u(4)*g0_ll(a,4)
              end do
              do b=1,4
                n_l_x(1,b)=-1/2.0d0/sqrt(-g0_uu(1,1))**3*g0_uu_x(1,1,b)
                s_l_x(1,b)=1/sqrt(g0_ll(2,2))*g0_ll_x(1,2,b)
     &            -1/2.0d0/sqrt(g0_ll(2,2))**3*g0_ll(1,2)*g0_ll_x(2,2,b)
                s_l_x(2,b)=1/2.0d0/sqrt(g0_ll(2,2))*g0_ll_x(2,2,b)
              end do
              do a=1,4
                do b=1,4
                  sigma_uu(a,b)=g0_uu(a,b)+n_u(a)*n_u(b)-s_u(a)*s_u(b)
                end do
              end do

              ! calculate thetap_ires function at point i
              thetap_ires(i)=0.0d0
              do c=1,4
                do d=1,4
                  thetap_ires(i)=thetap_ires(i)
     &                   +sigma_uu(c,d)*(n_l_x(d,c)+s_l_x(d,c))
                  do e=1,4
                    thetap_ires(i)=thetap_ires(i)
     &                   -sigma_uu(c,d)*gamma_ull(e,c,d)*(n_l(e)+s_l(e))
                  end do
                end do
              end do
              thetap_ads=2/x0*sqrt((1-x0**2)**2+x0**2/L**2)

              thetap_adsbhunitr0=2*sqrt(x0**3+(1-x0**2)**2*(x0**2+x0-1))
     &                           /x0/sqrt(x0)
              thetap_ires(i)=thetap_ires(i)/thetap_ads
              if (is_nan(thetap_ires(i))) thetap_ires(i)=-7.0d0

          end if

        end do

        return
        end

