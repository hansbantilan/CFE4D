c----------------------------------------------------------------------
c this routine initializes the past time level to O(h^3)
c given [gb,phi] and d[[gb,phi]]/dt at t=0
c
c note that at interior boundaries the spatial derivatives
c are calculated using backwards/forwards stencils
c (can't inject from the parent level, as they are out of sync)
c----------------------------------------------------------------------
        subroutine init_nm1(
     &                      db_txtx_np1,db_txtx_n,db_txtx_nm1,
     &                      db_txty_np1,db_txty_n,db_txty_nm1,
     &                      db_txtz_np1,db_txtz_n,db_txtz_nm1,
     &                      db_txxy_np1,db_txxy_n,db_txxy_nm1,
     &                      db_txxz_np1,db_txxz_n,db_txxz_nm1,
     &                      db_txyz_np1,db_txyz_n,db_txyz_nm1,
     &                      db_tyty_np1,db_tyty_n,db_tyty_nm1,
     &                      db_tytz_np1,db_tytz_n,db_tytz_nm1,
     &                      db_tyxy_np1,db_tyxy_n,db_tyxy_nm1,
     &                      db_tyxz_np1,db_tyxz_n,db_tyxz_nm1,
     &                      db_tyyz_np1,db_tyyz_n,db_tyyz_nm1,
     &                      db_tztz_np1,db_tztz_n,db_tztz_nm1,
     &                      db_tzxy_np1,db_tzxy_n,db_tzxy_nm1,
     &                      db_tzxz_np1,db_tzxz_n,db_tzxz_nm1,
     &                      db_tzyz_np1,db_tzyz_n,db_tzyz_nm1,
     &                      db_xyxy_np1,db_xyxy_n,db_xyxy_nm1,
     &                      db_xyxz_np1,db_xyxz_n,db_xyxz_nm1,
     &                      db_xyyz_np1,db_xyyz_n,db_xyyz_nm1,
     &                      db_xzxz_np1,db_xzxz_n,db_xzxz_nm1,
     &                      db_xzyz_np1,db_xzyz_n,db_xzyz_nm1,
     &                      db_yzyz_np1,db_yzyz_n,db_yzyz_nm1,
     &                      gb_tt_np1,gb_tt_n,gb_tt_nm1,gb_tt_t_n,
     &                      gb_tx_np1,gb_tx_n,gb_tx_nm1,gb_tx_t_n,
     &                      gb_xx_np1,gb_xx_n,gb_xx_nm1,gb_xx_t_n,
     &                      psi_np1,psi_n,psi_nm1,psi_t_n,
     &                      Hb_t_np1,Hb_t_n,Hb_t_nm1,Hb_t_t_n,
     &                      Hb_x_np1,Hb_x_n,Hb_x_nm1,Hb_x_t_n,
     &                      phi1_np1,phi1_n,phi1_nm1,phi1_t_n,
     &                      L,phys_bdy,x,dt,chr,ex,Nx)
        implicit none
        integer Nx
        integer phys_bdy(2)
        real*8 dt,ex,L
        real*8 chr(Nx)
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
        real*8 gb_tt_np1(Nx),gb_tt_n(Nx),gb_tt_nm1(Nx),gb_tt_t_n(Nx)
        real*8 gb_tx_np1(Nx),gb_tx_n(Nx),gb_tx_nm1(Nx),gb_tx_t_n(Nx)
        real*8 gb_xx_np1(Nx),gb_xx_n(Nx),gb_xx_nm1(Nx),gb_xx_t_n(Nx)
        real*8 Hb_t_np1(Nx),Hb_t_n(Nx),Hb_t_nm1(Nx),Hb_t_t_n(Nx)
        real*8 Hb_x_np1(Nx),Hb_x_n(Nx),Hb_x_nm1(Nx),Hb_x_t_n(Nx)
        real*8 psi_np1(Nx),psi_n(Nx),psi_nm1(Nx),psi_t_n(Nx)
        real*8 phi1_np1(Nx),phi1_n(Nx),phi1_nm1(Nx),phi1_t_n(Nx)

        real*8 x(Nx)

        real*8 lambda4
        real*8 tr_set

        logical is_nan

        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 gb_tt_t,gb_tt_tt
        real*8 gb_tx_t,gb_tx_tt
        real*8 gb_xx_t,gb_xx_tt
        real*8 psi_t,psi_tt
        real*8 Hb_t_t
        real*8 Hb_x_t
        real*8 phi1_t,phi1_tt

        real*8 h0_ll_tt(4,4),phi10_tt

        real*8 gammagg(4,4,4),gammahh(4,4,4)
        real*8 gammagh(4,4,4),gammahg(4,4,4)

        integer i,is,ie
        integer a,b,c,d

        real*8 x0
        real*8 dx

        real*8 dV1_dphi10 ! eventually have this as an input argument

        logical ltrace
        parameter (ltrace=.false.)

        real*8 PI
        parameter (PI=3.141592653589793d0)

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
        data i,is,ie/0,0,0/

        data x0/0.0/

        data dx/0.0/

        data gammagg,gammahh/64*0.0,64*0.0/
        data gammagh,gammahg/64*0.0,64*0.0/

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
        !---------------------------------------------------------------

        dx=x(2)-x(1)

        ! AdS4D cosmological constant
        !(lambda4=-(n-1)(n-2)/2/L^2) for n=4 dimensional AdS)
        lambda4=-3/L/L

        do i=1,Nx
          phi1_nm1(i)=phi1_n(i)
        end do

        is=1
        ie=Nx
        if (phys_bdy(1).eq.1) is=2
        if (phys_bdy(2).eq.1) ie=Nx-1

        do i=is,ie
          x0=x(i)

          if (chr(i).ne.ex) then

          !-----------------------------------------------------------
          ! some other initializion, which needs to be done before
          ! temporal derivatives are calculated
          !-----------------------------------------------------------

          ! computes tensors at point i 
          call tensor_init(
     &            db_txtx_np1,db_txtx_n,db_txtx_nm1,
     &            db_txty_np1,db_txty_n,db_txty_nm1,
     &            db_txtz_np1,db_txtz_n,db_txtz_nm1,
     &            db_txxy_np1,db_txxy_n,db_txxy_nm1,
     &            db_txxz_np1,db_txxz_n,db_txxz_nm1,
     &            db_txyz_np1,db_txyz_n,db_txyz_nm1,
     &            db_tyty_np1,db_tyty_n,db_tyty_nm1,
     &            db_tytz_np1,db_tytz_n,db_tytz_nm1,
     &            db_tyxy_np1,db_tyxy_n,db_tyxy_nm1,
     &            db_tyxz_np1,db_tyxz_n,db_tyxz_nm1,
     &            db_tyyz_np1,db_tyyz_n,db_tyyz_nm1,
     &            db_tztz_np1,db_tztz_n,db_tztz_nm1,
     &            db_tzxy_np1,db_tzxy_n,db_tzxy_nm1,
     &            db_tzxz_np1,db_tzxz_n,db_tzxz_nm1,
     &            db_tzyz_np1,db_tzyz_n,db_tzyz_nm1,
     &            db_xyxy_np1,db_xyxy_n,db_xyxy_nm1,
     &            db_xyxz_np1,db_xyxz_n,db_xyxz_nm1,
     &            db_xyyz_np1,db_xyyz_n,db_xyyz_nm1,
     &            db_xzxz_np1,db_xzxz_n,db_xzxz_nm1,
     &            db_xzyz_np1,db_xzyz_n,db_xzyz_nm1,
     &            db_yzyz_np1,db_yzyz_n,db_yzyz_nm1,
     &            gb_tt_n,gb_tt_n,gb_tt_n,
     &            gb_tx_n,gb_tx_n,gb_tx_n,
     &            gb_xx_n,gb_xx_n,gb_xx_n,
     &            psi_n,psi_n,psi_n,
     &            Hb_t_n,Hb_t_n,Hb_t_n,
     &            Hb_x_n,Hb_x_n,Hb_x_n,
     &            phi1_n,phi1_n,phi1_n,
     &            g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &            gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &            h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &            A_l,A_l_x,Hads_l,
     &            gamma_ull,gamma_ull_x,
     &            riemann_ulll,ricci_ll,ricci_lu,ricci,
     &            einstein_ll,set_ll,
     &            phi10_x,phi10_xx,
     &            x,dt,chr,L,ex,Nx,i)

          ! auxiliary objects
          do a=1,4
            do b=1,4
              do c=1,4
                gammagg(a,b,c)=0
                gammahh(a,b,c)=0
                gammagh(a,b,c)=0
                gammahg(a,b,c)=0
                do d=1,4
                  gammagg(a,b,c)=gammagg(a,b,c)
     &                           +0.5d0*gads_uu(a,d)
     &                                *(gads_ll_x(c,d,b)
     &                                 -gads_ll_x(b,c,d)
     &                                 +gads_ll_x(d,b,c))
                  gammahh(a,b,c)=gammahh(a,b,c)
     &                           +0.5d0*h0_uu(a,d)
     &                                *(h0_ll_x(c,d,b)
     &                                 -h0_ll_x(b,c,d)
     &                                 +h0_ll_x(d,b,c))
                  gammagh(a,b,c)=gammagh(a,b,c)
     &                           +0.5d0*gads_uu(a,d)
     &                                *(h0_ll_x(c,d,b)
     &                                 -h0_ll_x(b,c,d)
     &                                 +h0_ll_x(d,b,c))
                  gammahg(a,b,c)=gammahg(a,b,c)
     &                           +0.5d0*h0_uu(a,d)
     &                                *(gads_ll_x(c,d,b)
     &                                 -gads_ll_x(b,c,d)
     &                                 +gads_ll_x(d,b,c))
                end do
              end do
            end do
          end do

          tr_set =set_ll(1,1)*g0_uu(1,1)+
     &            set_ll(2,2)*g0_uu(2,2)+
     &            set_ll(3,3)*g0_uu(3,3)+
     &            set_ll(4,4)*g0_uu(4,4)+
     &         2*(set_ll(1,2)*g0_uu(1,2)+
     &            set_ll(1,3)*g0_uu(1,3)+
     &            set_ll(1,4)*g0_uu(1,4)+
     &            set_ll(2,3)*g0_uu(2,3)+
     &            set_ll(2,4)*g0_uu(2,4)+
     &            set_ll(3,4)*g0_uu(3,4))

          ! initial first time derivatives; gb_ii_t_n,Hb_i_t_n,phi1_t_n were set in CFE4D_free_data()

          ! need this in h0_ll_tt,phi10_tt calculations
          phi10_x(1)    =phi1_t_n(i)*(1-x0**2)**2   
          h0_ll_x(1,1,1)=gb_tt_t_n(i)  
          h0_ll_x(1,2,1)=gb_tx_t_n(i)
          h0_ll_x(2,2,1)=gb_xx_t_n(i)  
          h0_ll_x(3,3,1)=psi_t_n(i)  
          h0_ll_x(4,4,1)=psi_t_n(i)  
          A_l_x(1,1)    =Hb_t_t_n(i)*(1-x0**2)
          A_l_x(2,1)    =Hb_x_t_n(i)*(1-x0**2)

          ! need this in gb_ii_nm1/np1,Hb_i_nm1/np1,phi1_nm1/np1 updates
          phi1_t =phi1_t_n(i)                
          gb_tt_t=gb_tt_t_n(i) 
          gb_tx_t=gb_tx_t_n(i) 
          gb_xx_t=gb_xx_t_n(i) 
          psi_t  =psi_t_n(i) 
          Hb_t_t =Hb_t_t_n(i)
          Hb_x_t =Hb_x_t_n(i)

          ! 0 = efe_ab
          do a=1,3
            do b=a,3
              h0_ll_tt(a,b)=2/g0_uu(1,1)*
     &            ( 
     &                       -0.5d0*(                             
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
     &                       -0.5d0*(                             
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
     &                       -0.5d0*(                            
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
     &                       -0.5d0*A_l_x(a,b)                  
     &
     &                       -0.5d0*A_l_x(b,a)           
     &
     &                           +(           
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
     &                           -(
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
     &                       -lambda4*h0_ll(a,b)
     &
     &                       -8*PI*(set_ll(a,b)-tr_set*g0_ll(a,b)/2) 
     &            )
            end do
          end do          

          ! 0 = g^ab phi1,ab - g^ab gamma^c_ab phi1,c 
          phi10_tt=-1/g0_uu(1,1)
     &              *(    
     &                    phi10_xx(2,2)*g0_uu(2,2)+
     &                    phi10_xx(3,3)*g0_uu(3,3)+
     &                    phi10_xx(4,4)*g0_uu(4,4)+
     &                 2*(phi10_xx(1,2)*g0_uu(1,2)+
     &                    phi10_xx(1,3)*g0_uu(1,3)+
     &                    phi10_xx(1,4)*g0_uu(1,4)+
     &                    phi10_xx(2,3)*g0_uu(2,3)+
     &                    phi10_xx(2,4)*g0_uu(2,4)+
     &                    phi10_xx(3,4)*g0_uu(3,4))
     &                -
     &                    phi10_x(1)*( gamma_ull(1,1,1)*g0_uu(1,1)+
     &                                 gamma_ull(1,2,2)*g0_uu(2,2)+
     &                                 gamma_ull(1,3,3)*g0_uu(3,3)+
     &                                 gamma_ull(1,4,4)*g0_uu(4,4)+
     &                              2*(gamma_ull(1,1,2)*g0_uu(1,2)+
     &                                 gamma_ull(1,1,3)*g0_uu(1,3)+
     &                                 gamma_ull(1,1,4)*g0_uu(1,4)+
     &                                 gamma_ull(1,2,3)*g0_uu(2,3)+
     &                                 gamma_ull(1,2,4)*g0_uu(2,4)+
     &                                 gamma_ull(1,3,4)*g0_uu(3,4)) )
     &                -
     &                    phi10_x(2)*( gamma_ull(2,1,1)*g0_uu(1,1)+
     &                                 gamma_ull(2,2,2)*g0_uu(2,2)+
     &                                 gamma_ull(2,3,3)*g0_uu(3,3)+
     &                                 gamma_ull(2,4,4)*g0_uu(4,4)+
     &                              2*(gamma_ull(2,1,2)*g0_uu(1,2)+
     &                                 gamma_ull(2,1,3)*g0_uu(1,3)+
     &                                 gamma_ull(2,1,4)*g0_uu(1,4)+
     &                                 gamma_ull(2,2,3)*g0_uu(2,3)+
     &                                 gamma_ull(2,2,4)*g0_uu(2,4)+
     &                                 gamma_ull(2,3,4)*g0_uu(3,4)) )
     &                -
     &                    phi10_x(3)*( gamma_ull(3,1,1)*g0_uu(1,1)+
     &                                 gamma_ull(3,2,2)*g0_uu(2,2)+
     &                                 gamma_ull(3,3,3)*g0_uu(3,3)+
     &                                 gamma_ull(3,4,4)*g0_uu(4,4)+
     &                              2*(gamma_ull(3,1,2)*g0_uu(1,2)+
     &                                 gamma_ull(3,1,3)*g0_uu(1,3)+
     &                                 gamma_ull(3,1,4)*g0_uu(1,4)+
     &                                 gamma_ull(3,2,3)*g0_uu(2,3)+
     &                                 gamma_ull(3,2,4)*g0_uu(2,4)+
     &                                 gamma_ull(3,3,4)*g0_uu(3,4)) )
     &                -
     &                    phi10_x(4)*( gamma_ull(4,1,1)*g0_uu(1,1)+
     &                                 gamma_ull(4,2,2)*g0_uu(2,2)+
     &                                 gamma_ull(4,3,3)*g0_uu(3,3)+
     &                                 gamma_ull(4,4,4)*g0_uu(4,4)+
     &                              2*(gamma_ull(4,1,2)*g0_uu(1,2)+
     &                                 gamma_ull(4,1,3)*g0_uu(1,3)+
     &                                 gamma_ull(4,1,4)*g0_uu(1,4)+
     &                                 gamma_ull(4,2,3)*g0_uu(2,3)+
     &                                 gamma_ull(4,2,4)*g0_uu(2,4)+
     &                                 gamma_ull(4,3,4)*g0_uu(3,4)) )
     &                -
     &                    dV1_dphi10             
     &                                                                 )

          if (is_nan(h0_ll_tt(1,1)).or.is_nan(h0_ll_tt(1,2)) 
     &    .or.is_nan(h0_ll_tt(2,2)).or.is_nan(h0_ll_tt(3,3))) then
            write(*,*) 'h0_ll_tt(1,1)=',h0_ll_tt(1,1)
            write(*,*) 'h0_ll_tt(1,2)=',h0_ll_tt(1,2)
            write(*,*) 'h0_ll_tt(2,2)=',h0_ll_tt(2,2)
            write(*,*) 'h0_ll_tt(3,3)=',h0_ll_tt(3,3)
            stop
          end if

          ! initial second time derivatives
          gb_tt_tt=h0_ll_tt(1,1) 
          gb_tx_tt=h0_ll_tt(1,2)
          gb_xx_tt=h0_ll_tt(2,2) 
          psi_tt  =h0_ll_tt(3,3)/x0**2 
          phi1_tt =phi10_tt/(1-x0**2)**2

          ! initialize past time level by O(h^3) expansion
           gb_tt_nm1(i)=gb_tt_n(i) - gb_tt_t*dt
     &                      + gb_tt_tt*dt**2/2
           gb_tx_nm1(i)=gb_tx_n(i) - gb_tx_t*dt
     &                      + gb_tx_tt*dt**2/2
           gb_xx_nm1(i)=gb_xx_n(i) - gb_xx_t*dt
     &                      + gb_xx_tt*dt**2/2
           psi_nm1(i)  =psi_n(i) - psi_t*dt
     &                      + psi_tt*dt**2/2
           Hb_t_nm1(i) =Hb_t_n(i) - Hb_t_t*dt
           Hb_x_nm1(i) =Hb_x_n(i) - Hb_x_t*dt
           phi1_nm1(i) =phi1_n(i) - phi1_t*dt
     &                      + phi1_tt*dt**2/2
        
          ! initialize future time level by O(h^3) expansion
           gb_tt_np1(i)=gb_tt_n(i) + gb_tt_t*dt
     &                      + gb_tt_tt*dt**2/2
           gb_tx_np1(i)=gb_tx_n(i) + gb_tx_t*dt
     &                      + gb_tx_tt*dt**2/2
           gb_xx_np1(i)=gb_xx_n(i) + gb_xx_t*dt
     &                      + gb_xx_tt*dt**2/2
           psi_np1(i)  =psi_n(i) + psi_t*dt
     &                      + psi_tt*dt**2/2
           Hb_t_np1(i) =Hb_t_n(i) + Hb_t_t*dt
           Hb_x_np1(i) =Hb_x_n(i) + Hb_x_t*dt     
           phi1_np1(i) =phi1_n(i) + phi1_t*dt
     &                      + phi1_tt*dt**2/2
  
          end if

        end do

        call axi_reg_phi(phi1_nm1,chr,ex,L,x,Nx)
        call axi_reg_phi(phi1_np1,chr,ex,L,x,Nx)
        call axi_reg_g(gb_tt_nm1,gb_tx_nm1,gb_xx_nm1,psi_nm1,
     &                 chr,ex,L,x,Nx)
        call axi_reg_g(gb_tt_np1,gb_tx_np1,gb_xx_np1,psi_np1,
     &                 chr,ex,L,x,Nx)

        return
        end
