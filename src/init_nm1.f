c----------------------------------------------------------------------
c this routine initializes the past time level to O(h^3),
c given [gb,phi] at t=t0, and d[[gb,phi]]/dt
c
c note that at interior boundaries the spatial derivatives
c are calculated using backwards/forwards stencils. (can't 
c inject from the parent level, as they are out of sync)
c----------------------------------------------------------------------
        subroutine init_nm1(gb_tt_np1,gb_tt_n,gb_tt_nm1,gb_tt_t_n,
     &                      gb_tx_np1,gb_tx_n,gb_tx_nm1,gb_tx_t_n,
     &                      gb_ty_np1,gb_ty_n,gb_ty_nm1,gb_ty_t_n,
     &                      gb_xx_np1,gb_xx_n,gb_xx_nm1,gb_xx_t_n,
     &                      gb_xy_np1,gb_xy_n,gb_xy_nm1,gb_xy_t_n,
     &                      gb_yy_np1,gb_yy_n,gb_yy_nm1,gb_yy_t_n,
     &                      psi_np1,psi_n,psi_nm1,psi_t_n,
     &                      Hb_t_np1,Hb_t_n,Hb_t_nm1,Hb_t_t_n,
     &                      Hb_x_np1,Hb_x_n,Hb_x_nm1,Hb_x_t_n,
     &                      Hb_y_np1,Hb_y_n,Hb_y_nm1,Hb_y_t_n,
     &                      phi1_np1,phi1_n,phi1_nm1,phi1_t_n,tfunction,
     &                      L,phys_bdy,x,y,dt,chr,ex,Nx,Ny,regtype)
        implicit none
        integer Nx,Ny
        integer phys_bdy(4)
        integer regtype
        real*8 dt,ex,L
        real*8 chr(Nx,Ny)
        real*8 gb_tt_np1(Nx,Ny),gb_tt_n(Nx,Ny),gb_tt_nm1(Nx,Ny)
        real*8 gb_tt_t_n(Nx,Ny)
        real*8 gb_tx_np1(Nx,Ny),gb_tx_n(Nx,Ny),gb_tx_nm1(Nx,Ny)
        real*8 gb_tx_t_n(Nx,Ny)
        real*8 gb_ty_np1(Nx,Ny),gb_ty_n(Nx,Ny),gb_ty_nm1(Nx,Ny)
        real*8 gb_ty_t_n(Nx,Ny)
        real*8 gb_xx_np1(Nx,Ny),gb_xx_n(Nx,Ny),gb_xx_nm1(Nx,Ny)
        real*8 gb_xx_t_n(Nx,Ny)
        real*8 gb_xy_np1(Nx,Ny),gb_xy_n(Nx,Ny),gb_xy_nm1(Nx,Ny)
        real*8 gb_xy_t_n(Nx,Ny)
        real*8 gb_yy_np1(Nx,Ny),gb_yy_n(Nx,Ny),gb_yy_nm1(Nx,Ny)
        real*8 gb_yy_t_n(Nx,Ny)
        real*8 psi_np1(Nx,Ny),psi_n(Nx,Ny),psi_nm1(Nx,Ny)
        real*8 psi_t_n(Nx,Ny)
        real*8 Hb_t_np1(Nx,Ny),Hb_t_n(Nx,Ny),Hb_t_nm1(Nx,Ny)
        real*8 Hb_t_t_n(Nx,Ny)
        real*8 Hb_x_np1(Nx,Ny),Hb_x_n(Nx,Ny),Hb_x_nm1(Nx,Ny)
        real*8 Hb_x_t_n(Nx,Ny)
        real*8 Hb_y_np1(Nx,Ny),Hb_y_n(Nx,Ny),Hb_y_nm1(Nx,Ny)
        real*8 Hb_y_t_n(Nx,Ny)
        real*8 phi1_np1(Nx,Ny),phi1_n(Nx,Ny),phi1_nm1(Nx,Ny)
        real*8 phi1_t_n(Nx,Ny)

        real*8 tfunction(Nx,Ny)

        real*8 x(Nx),y(Ny)

        real*8 lambda5
        real*8 tr_set

        logical is_nan

        real*8 term1sub(5,5),term2(5,5),term3(5,5),term4(5,5)
        real*8 term5(5,5),term6(5,5),term7(5,5),term8(5,5)

        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 gb_tt_t0,gb_tt_tt0
        real*8 gb_tx_t0,gb_tx_tt0
        real*8 gb_ty_t0,gb_ty_tt0
        real*8 gb_xx_t0,gb_xx_tt0
        real*8 gb_xy_t0,gb_xy_tt0
        real*8 gb_yy_t0,gb_yy_tt0
        real*8 psi_t0,psi_tt0
        real*8 Hb_t_t0
        real*8 Hb_x_t0
        real*8 Hb_y_t0
        real*8 phi1_t0,phi1_tt0

        real*8 h0_ll_tt(5,5),phi10_tt

        real*8 gammagg(5,5,5),gammahh(5,5,5)
        real*8 gammagh(5,5,5),gammahg(5,5,5)

        integer a,b,c,d

        integer i,j,is,ie,js,je
        real*8 dx,dy
        real*8 x0,y0
        real*8 rho0

        real*8 dV1_dphi10 ! EVENTUALLY NEED THIS AS AN INPUT ARGUMENT

        real*8 PI
        parameter (PI=3.141592653589793d0)

        logical ltrace
        parameter (ltrace=.false.)

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,w,y,z)
        !--------------------------------------------------------------
        real*8 g0_ll(5,5),g0_uu(5,5)
        real*8 g0_ll_x(5,5,5),g0_uu_x(5,5,5),g0_ll_xx(5,5,5,5)
        real*8 gads_ll(5,5),gads_uu(5,5)
        real*8 gads_ll_x(5,5,5),gads_uu_x(5,5,5),gads_ll_xx(5,5,5,5)
        real*8 h0_ll(5,5),h0_uu(5,5)
        real*8 h0_ll_x(5,5,5),h0_uu_x(5,5,5),h0_ll_xx(5,5,5,5)
        real*8 gamma_ull(5,5,5),gamma_ull_x(5,5,5,5)
        real*8 riemann_ulll(5,5,5,5)
        real*8 ricci_ll(5,5),ricci_lu(5,5),ricci
        real*8 einstein_ll(5,5),set_ll(5,5)
        real*8 Hads_l(5),A_l(5),A_l_x(5,5)
        real*8 phi10_x(5),phi10_xx(5,5)

        !--------------------------------------------------------------
        ! initialize fixed-size variables
        !--------------------------------------------------------------
        data gammagg,gammahh/125*0.0,125*0.0/
        data gammagh,gammahg/125*0.0,125*0.0/

        data g0_ll,g0_uu/25*0.0,25*0.0/
        data gads_ll,gads_uu/25*0.0,25*0.0/
        data h0_ll,h0_uu/25*0.0,25*0.0/
        data gamma_ull/125*0.0/
        data gamma_ull_x/625*0.0/

        data g0_ll_x,g0_uu_x/125*0.0,125*0.0/
        data gads_ll_x,gads_uu_x/125*0.0,125*0.0/
        data h0_ll_x,h0_uu_x/125*0.0,125*0.0/

        data g0_ll_xx/625*0.0/
        data gads_ll_xx/625*0.0/
        data h0_ll_xx/625*0.0/

        data ricci/0.0/
        data ricci_ll,ricci_lu/25*0.0,25*0.0/
        data einstein_ll,set_ll/25*0.0,25*0.0/
        data riemann_ulll/625*0.0/

        data A_l,Hads_l/5*0.0,5*0.0/
        data A_l_x/25*0.0/

        data phi10_x/5*0.0/
        data phi10_xx/25*0.0/

        data term1sub,term2/25*0.0,25*0.0/
        data term3,term4/25*0.0,25*0.0/
        data term5,term6/25*0.0,25*0.0/
        data term7,term8/25*0.0,25*0.0/
        !--------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        ! CFE4D cosmological constant
        !(lambda5=-(n-1)(n-2)/L^2) for n=5 dimensional AdS)
        lambda5=-6/L/L

        is=1
        js=1
        ie=Nx
        je=Ny
        if (phys_bdy(1).eq.1) is=2
        if (phys_bdy(3).eq.1) js=2
        if (phys_bdy(2).eq.1) ie=Nx-1
        if (phys_bdy(4).eq.1) je=Ny-1

        do i=is,ie
          do j=js,je
            x0=x(i)
            y0=y(j)
            rho0=sqrt(x0**2+y0**2)

            if (chr(i,j).ne.ex) then

              !-----------------------------------------------------------
              ! some other initializion, which needs to be done before
              ! temporal derivatives are calculated
              !-----------------------------------------------------------

              ! computes tensors at point i,j
              call tensor_init(
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
     &                g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &                gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &                h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                A_l,A_l_x,Hads_l,
     &                gamma_ull,gamma_ull_x,
     &                riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                einstein_ll,set_ll,
     &                phi10_x,phi10_xx,
     &                x,y,dt,chr,L,ex,Nx,Ny,i,j)

              ! auxiliary objects
              do a=1,5
                do b=1,5
                  do c=1,5
                    gammagg(a,b,c)=0
                    gammahh(a,b,c)=0
                    gammagh(a,b,c)=0
                    gammahg(a,b,c)=0
                    do d=1,5
                      gammagg(a,b,c)=gammagg(a,b,c)
     &                               +0.5d0*gads_uu(a,d)
     &                                    *(gads_ll_x(c,d,b)
     &                                     -gads_ll_x(b,c,d)
     &                                     +gads_ll_x(d,b,c))
                      gammahh(a,b,c)=gammahh(a,b,c)
     &                               +0.5d0*h0_uu(a,d)
     &                                    *(h0_ll_x(c,d,b)
     &                                     -h0_ll_x(b,c,d)
     &                                     +h0_ll_x(d,b,c))
                      gammagh(a,b,c)=gammagh(a,b,c)
     &                               +0.5d0*gads_uu(a,d)
     &                                    *(h0_ll_x(c,d,b)
     &                                     -h0_ll_x(b,c,d)
     &                                     +h0_ll_x(d,b,c))
                      gammahg(a,b,c)=gammahg(a,b,c)
     &                               +0.5d0*h0_uu(a,d)
     &                                    *(gads_ll_x(c,d,b)
     &                                     -gads_ll_x(b,c,d)
     &                                     +gads_ll_x(d,b,c))
                    end do
                  end do
                end do
              end do

              tr_set =set_ll(1,1)*g0_uu(1,1)+
     &                set_ll(2,2)*g0_uu(2,2)+
     &                set_ll(3,3)*g0_uu(3,3)+
     &                set_ll(4,4)*g0_uu(4,4)+
     &                set_ll(5,5)*g0_uu(5,5)+
     &             2*(set_ll(1,2)*g0_uu(1,2)+
     &                set_ll(1,3)*g0_uu(1,3)+
     &                set_ll(1,4)*g0_uu(1,4)+
     &                set_ll(1,5)*g0_uu(1,5)+
     &                set_ll(2,3)*g0_uu(2,3)+
     &                set_ll(2,4)*g0_uu(2,4)+
     &                set_ll(2,5)*g0_uu(2,5)+
     &                set_ll(3,4)*g0_uu(3,4)+
     &                set_ll(3,5)*g0_uu(3,5)+
     &                set_ll(4,5)*g0_uu(4,5))

              ! initial first time derivatives;
              ! gb_ii_t_n,Hb_i_t_n,phi1_t_n are set in init_ghbdot.f

              ! need this in h0_ll_tt,phi10_tt calculations
              phi10_x(1)    =phi1_t_n(i,j)*(1-rho0**2)**3
              h0_ll_x(1,1,1)=gb_tt_t_n(i,j)*(1-rho0**2)
              h0_ll_x(1,2,1)=gb_tx_t_n(i,j)*(1-rho0**2)
              h0_ll_x(1,3,1)=gb_ty_t_n(i,j)*(1-rho0**2)
              h0_ll_x(2,2,1)=gb_xx_t_n(i,j)*(1-rho0**2)
              h0_ll_x(2,3,1)=gb_xy_t_n(i,j)*(1-rho0**2)
              h0_ll_x(3,3,1)=gb_yy_t_n(i,j)*(1-rho0**2)
              h0_ll_x(4,4,1)=psi_t_n(i,j)*(1-rho0**2)*y0**2
              h0_ll_x(5,5,1)=psi_t_n(i,j)*(1-rho0**2)*y0**2
              A_l_x(1,1)    =Hb_t_t_n(i,j)*(1-rho0**2)**2
              A_l_x(2,1)    =Hb_x_t_n(i,j)*(1-rho0**2)**2
              A_l_x(3,1)    =Hb_y_t_n(i,j)*(1-rho0**2)**2

              ! need this in gb_ii_nm1/np1,Hb_i_nm1/np1,phi1_nm1/np1 updates
              phi1_t0 =phi1_t_n(i,j)*(1-rho0**2)**3 
              gb_tt_t0=gb_tt_t_n(i,j)*(1-rho0**2)
              gb_tx_t0=gb_tx_t_n(i,j)*(1-rho0**2)
              gb_ty_t0=gb_ty_t_n(i,j)*(1-rho0**2)
              gb_xx_t0=gb_xx_t_n(i,j)*(1-rho0**2)
              gb_xy_t0=gb_xy_t_n(i,j)*(1-rho0**2)
              gb_yy_t0=gb_yy_t_n(i,j)*(1-rho0**2)
              psi_t0  =psi_t_n(i,j)*(1-rho0**2)*y0**2
              Hb_t_t0 =Hb_t_t_n(i,j)*(1-rho0**2)**2
              Hb_x_t0 =Hb_x_t_n(i,j)*(1-rho0**2)**2
              Hb_y_t0 =Hb_y_t_n(i,j)*(1-rho0**2)**2
                       
              ! 0 = efe_ab
              do a=1,5
                do b=a,5
                  term1sub(a,b)=-0.5d0*(                             
     &                              h0_uu(2,2)*h0_ll_xx(a,b,2,2)+
     &                              h0_uu(3,3)*h0_ll_xx(a,b,3,3)+
     &                              h0_uu(4,4)*h0_ll_xx(a,b,4,4)+
     &                              h0_uu(5,5)*h0_ll_xx(a,b,5,5)+
     &                           2*(h0_uu(1,2)*h0_ll_xx(a,b,1,2)+
     &                              h0_uu(1,3)*h0_ll_xx(a,b,1,3)+
     &                              h0_uu(1,4)*h0_ll_xx(a,b,1,4)+
     &                              h0_uu(1,5)*h0_ll_xx(a,b,1,5)+
     &                              h0_uu(2,3)*h0_ll_xx(a,b,2,3)+
     &                              h0_uu(2,4)*h0_ll_xx(a,b,2,4)+
     &                              h0_uu(2,5)*h0_ll_xx(a,b,2,5)+
     &                              h0_uu(3,4)*h0_ll_xx(a,b,3,4)+
     &                              h0_uu(3,5)*h0_ll_xx(a,b,3,5)+
     &                              h0_uu(4,5)*h0_ll_xx(a,b,4,5))
     &                           +
     &                              gads_uu(2,2)*h0_ll_xx(a,b,2,2)+
     &                              gads_uu(3,3)*h0_ll_xx(a,b,3,3)+
     &                              gads_uu(4,4)*h0_ll_xx(a,b,4,4)+
     &                              gads_uu(5,5)*h0_ll_xx(a,b,5,5)+
     &                           2*(gads_uu(1,2)*h0_ll_xx(a,b,1,2)+
     &                              gads_uu(1,3)*h0_ll_xx(a,b,1,3)+
     &                              gads_uu(1,4)*h0_ll_xx(a,b,1,4)+
     &                              gads_uu(1,5)*h0_ll_xx(a,b,1,5)+
     &                              gads_uu(2,3)*h0_ll_xx(a,b,2,3)+
     &                              gads_uu(2,4)*h0_ll_xx(a,b,2,4)+
     &                              gads_uu(2,5)*h0_ll_xx(a,b,2,5)+
     &                              gads_uu(3,4)*h0_ll_xx(a,b,3,4)+
     &                              gads_uu(3,5)*h0_ll_xx(a,b,3,5)+
     &                              gads_uu(4,5)*h0_ll_xx(a,b,4,5))
     &                           +
     &                              h0_uu(1,1)*gads_ll_xx(a,b,1,1)+
     &                              h0_uu(2,2)*gads_ll_xx(a,b,2,2)+
     &                              h0_uu(3,3)*gads_ll_xx(a,b,3,3)+
     &                              h0_uu(4,4)*gads_ll_xx(a,b,4,4)+
     &                              h0_uu(5,5)*gads_ll_xx(a,b,5,5)+
     &                           2*(h0_uu(1,2)*gads_ll_xx(a,b,1,2)+
     &                              h0_uu(1,3)*gads_ll_xx(a,b,1,3)+
     &                              h0_uu(1,4)*gads_ll_xx(a,b,1,4)+
     &                              h0_uu(1,5)*gads_ll_xx(a,b,1,5)+
     &                              h0_uu(2,3)*gads_ll_xx(a,b,2,3)+
     &                              h0_uu(2,4)*gads_ll_xx(a,b,2,4)+
     &                              h0_uu(2,5)*gads_ll_xx(a,b,2,5)+
     &                              h0_uu(3,4)*gads_ll_xx(a,b,3,4)+
     &                              h0_uu(3,5)*gads_ll_xx(a,b,3,5)+
     &                              h0_uu(4,5)*gads_ll_xx(a,b,4,5))
     &                                  )
     
                  term2(a,b)=    -0.5d0*(                             
     &                              h0_uu_x(1,1,a)* h0_ll_x(b,1,1) +
     &                              h0_uu_x(1,2,a)*(h0_ll_x(b,1,2) +
     &                                              h0_ll_x(b,2,1))+
     &                              h0_uu_x(1,3,a)*(h0_ll_x(b,1,3) +
     &                                              h0_ll_x(b,3,1))+
     &                              h0_uu_x(1,4,a)*(h0_ll_x(b,1,4) +
     &                                              h0_ll_x(b,4,1))+
     &                              h0_uu_x(1,5,a)*(h0_ll_x(b,1,5) +
     &                                              h0_ll_x(b,5,1))+
     &                              h0_uu_x(2,2,a)* h0_ll_x(b,2,2) +
     &                              h0_uu_x(2,3,a)*(h0_ll_x(b,2,3) +
     &                                              h0_ll_x(b,3,2))+
     &                              h0_uu_x(2,4,a)*(h0_ll_x(b,2,4) +
     &                                              h0_ll_x(b,4,2))+
     &                              h0_uu_x(2,5,a)*(h0_ll_x(b,2,5) +
     &                                              h0_ll_x(b,5,2))+
     &                              h0_uu_x(3,3,a)* h0_ll_x(b,3,3) +
     &                              h0_uu_x(3,4,a)*(h0_ll_x(b,3,4) +
     &                                              h0_ll_x(b,4,3))+
     &                              h0_uu_x(3,5,a)*(h0_ll_x(b,3,5) +
     &                                              h0_ll_x(b,5,3))+
     &                              h0_uu_x(4,4,a)* h0_ll_x(b,4,4) +
     &                              h0_uu_x(4,5,a)*(h0_ll_x(b,4,5) +
     &                                              h0_ll_x(b,5,4))+
     &                              h0_uu_x(5,5,a)* h0_ll_x(b,5,5)
     &                           +
     &                              gads_uu_x(1,1,a)* h0_ll_x(b,1,1) +
     &                              gads_uu_x(1,2,a)*(h0_ll_x(b,1,2) +
     &                                                h0_ll_x(b,2,1))+
     &                              gads_uu_x(1,3,a)*(h0_ll_x(b,1,3) +
     &                                                h0_ll_x(b,3,1))+
     &                              gads_uu_x(1,4,a)*(h0_ll_x(b,1,4) +
     &                                                h0_ll_x(b,4,1))+
     &                              gads_uu_x(1,5,a)*(h0_ll_x(b,1,5) +
     &                                                h0_ll_x(b,5,1))+
     &                              gads_uu_x(2,2,a)* h0_ll_x(b,2,2) +
     &                              gads_uu_x(2,3,a)*(h0_ll_x(b,2,3) +
     &                                                h0_ll_x(b,3,2))+
     &                              gads_uu_x(2,4,a)*(h0_ll_x(b,2,4) +
     &                                                h0_ll_x(b,4,2))+
     &                              gads_uu_x(2,5,a)*(h0_ll_x(b,2,5) +
     &                                                h0_ll_x(b,5,2))+
     &                              gads_uu_x(3,3,a)* h0_ll_x(b,3,3) +
     &                              gads_uu_x(3,4,a)*(h0_ll_x(b,3,4) +
     &                                                h0_ll_x(b,4,3))+
     &                              gads_uu_x(3,5,a)*(h0_ll_x(b,3,5) +
     &                                                h0_ll_x(b,5,3))+
     &                              gads_uu_x(4,4,a)* h0_ll_x(b,4,4) +
     &                              gads_uu_x(4,5,a)*(h0_ll_x(b,4,5) +
     &                                                h0_ll_x(b,5,4))+
     &                              gads_uu_x(5,5,a)* h0_ll_x(b,5,5)
     &                           +
     &                              h0_uu_x(1,1,a)* gads_ll_x(b,1,1) +
     &                              h0_uu_x(1,2,a)*(gads_ll_x(b,1,2) + 
     &                                              gads_ll_x(b,2,1))+ 
     &                              h0_uu_x(1,3,a)*(gads_ll_x(b,1,3) + 
     &                                              gads_ll_x(b,3,1))+ 
     &                              h0_uu_x(1,4,a)*(gads_ll_x(b,1,4) +
     &                                              gads_ll_x(b,4,1))+
     &                              h0_uu_x(1,5,a)*(gads_ll_x(b,1,5) +
     &                                              gads_ll_x(b,5,1))+
     &                              h0_uu_x(2,2,a)* gads_ll_x(b,2,2) +
     &                              h0_uu_x(2,3,a)*(gads_ll_x(b,2,3) +
     &                                              gads_ll_x(b,3,2))+
     &                              h0_uu_x(2,4,a)*(gads_ll_x(b,2,4) +
     &                                              gads_ll_x(b,4,2))+
     &                              h0_uu_x(2,5,a)*(gads_ll_x(b,2,5) +
     &                                              gads_ll_x(b,5,2))+
     &                              h0_uu_x(3,3,a)* gads_ll_x(b,3,3) +
     &                              h0_uu_x(3,4,a)*(gads_ll_x(b,3,4) +
     &                                              gads_ll_x(b,4,3))+
     &                              h0_uu_x(3,5,a)*(gads_ll_x(b,3,5) +
     &                                              gads_ll_x(b,5,3))+
     &                              h0_uu_x(4,4,a)* gads_ll_x(b,4,4) +
     &                              h0_uu_x(4,5,a)*(gads_ll_x(b,4,5) +
     &                                              gads_ll_x(b,5,4))+
     &                              h0_uu_x(5,5,a)* gads_ll_x(b,5,5)
     &                                )
     
                  term3(a,b)=    -0.5d0*(                            
     &                              h0_uu_x(1,1,b)* h0_ll_x(a,1,1) +
     &                              h0_uu_x(1,2,b)*(h0_ll_x(a,1,2) +
     &                                              h0_ll_x(a,2,1))+
     &                              h0_uu_x(1,3,b)*(h0_ll_x(a,1,3) +
     &                                              h0_ll_x(a,3,1))+
     &                              h0_uu_x(1,4,b)*(h0_ll_x(a,1,4) +
     &                                              h0_ll_x(a,4,1))+
     &                              h0_uu_x(1,5,b)*(h0_ll_x(a,1,5) +
     &                                              h0_ll_x(a,5,1))+
     &                              h0_uu_x(2,2,b)* h0_ll_x(a,2,2) +
     &                              h0_uu_x(2,3,b)*(h0_ll_x(a,2,3) +
     &                                              h0_ll_x(a,3,2))+
     &                              h0_uu_x(2,4,b)*(h0_ll_x(a,2,4) +
     &                                              h0_ll_x(a,4,2))+
     &                              h0_uu_x(2,5,b)*(h0_ll_x(a,2,5) +
     &                                              h0_ll_x(a,5,2))+
     &                              h0_uu_x(3,3,b)* h0_ll_x(a,3,3) +
     &                              h0_uu_x(3,4,b)*(h0_ll_x(a,3,4) +
     &                                              h0_ll_x(a,4,3))+
     &                              h0_uu_x(3,5,b)*(h0_ll_x(a,3,5) +
     &                                              h0_ll_x(a,5,3))+
     &                              h0_uu_x(4,4,b)* h0_ll_x(a,4,4) +
     &                              h0_uu_x(4,5,b)*(h0_ll_x(a,4,5) +
     &                                              h0_ll_x(a,5,4))+
     &                              h0_uu_x(5,5,b)* h0_ll_x(a,5,5)
     &                           +
     &                              gads_uu_x(1,1,b)* h0_ll_x(a,1,1) +
     &                              gads_uu_x(1,2,b)*(h0_ll_x(a,1,2) +
     &                                                h0_ll_x(a,2,1))+
     &                              gads_uu_x(1,3,b)*(h0_ll_x(a,1,3) +
     &                                                h0_ll_x(a,3,1))+
     &                              gads_uu_x(1,4,b)*(h0_ll_x(a,1,4) +
     &                                                h0_ll_x(a,4,1))+
     &                              gads_uu_x(1,5,b)*(h0_ll_x(a,1,5) +
     &                                                h0_ll_x(a,5,1))+
     &                              gads_uu_x(2,2,b)* h0_ll_x(a,2,2) +
     &                              gads_uu_x(2,3,b)*(h0_ll_x(a,2,3) +
     &                                                h0_ll_x(a,3,2))+
     &                              gads_uu_x(2,4,b)*(h0_ll_x(a,2,4) +
     &                                                h0_ll_x(a,4,2))+
     &                              gads_uu_x(2,5,b)*(h0_ll_x(a,2,5) +
     &                                                h0_ll_x(a,5,2))+
     &                              gads_uu_x(3,3,b)* h0_ll_x(a,3,3) +
     &                              gads_uu_x(3,4,b)*(h0_ll_x(a,3,4) +
     &                                                h0_ll_x(a,4,3))+
     &                              gads_uu_x(3,5,b)*(h0_ll_x(a,3,5) +
     &                                                h0_ll_x(a,5,3))+
     &                              gads_uu_x(4,4,b)* h0_ll_x(a,4,4) +
     &                              gads_uu_x(4,5,b)*(h0_ll_x(a,4,5) +
     &                                                h0_ll_x(a,5,4))+
     &                              gads_uu_x(5,5,b)* h0_ll_x(a,5,5)
     &                           +
     &                              h0_uu_x(1,1,b)* gads_ll_x(a,1,1) +
     &                              h0_uu_x(1,2,b)*(gads_ll_x(a,1,2) +  
     &                                              gads_ll_x(a,2,1))+  
     &                              h0_uu_x(1,3,b)*(gads_ll_x(a,1,3) +  
     &                                              gads_ll_x(a,3,1))+  
     &                              h0_uu_x(1,4,b)*(gads_ll_x(a,1,4) +
     &                                              gads_ll_x(a,4,1))+
     &                              h0_uu_x(1,5,b)*(gads_ll_x(a,1,5) +
     &                                              gads_ll_x(a,5,1))+
     &                              h0_uu_x(2,2,b)* gads_ll_x(a,2,2) +
     &                              h0_uu_x(2,3,b)*(gads_ll_x(a,2,3) +
     &                                              gads_ll_x(a,3,2))+
     &                              h0_uu_x(2,4,b)*(gads_ll_x(a,2,4) +
     &                                              gads_ll_x(a,4,2))+
     &                              h0_uu_x(2,5,b)*(gads_ll_x(a,2,5) +
     &                                              gads_ll_x(a,5,2))+
     &                              h0_uu_x(3,3,b)* gads_ll_x(a,3,3) +
     &                              h0_uu_x(3,4,b)*(gads_ll_x(a,3,4) +
     &                                              gads_ll_x(a,4,3))+
     &                              h0_uu_x(3,5,b)*(gads_ll_x(a,3,5) +
     &                                              gads_ll_x(a,5,3))+
     &                              h0_uu_x(4,4,b)* gads_ll_x(a,4,4) +
     &                              h0_uu_x(4,5,b)*(gads_ll_x(a,4,5) +
     &                                              gads_ll_x(a,5,4))+
     &                              h0_uu_x(5,5,b)* gads_ll_x(a,5,5)
     &                                )

                  term4(a,b)=    -0.5d0*A_l_x(a,b)                  

                  term5(a,b)=    -0.5d0*A_l_x(b,a)           

                  term6(a,b)=        +(           
     &                              Hads_l(1)*gammahh(1,a,b)+      
     &                              Hads_l(2)*gammahh(2,a,b)+
     &                              Hads_l(3)*gammahh(3,a,b)+
     &                              Hads_l(4)*gammahh(4,a,b)+
     &                              Hads_l(5)*gammahh(5,a,b)
     &                           +
     &                              Hads_l(1)*gammagh(1,a,b)+
     &                              Hads_l(2)*gammagh(2,a,b)+
     &                              Hads_l(3)*gammagh(3,a,b)+
     &                              Hads_l(4)*gammagh(4,a,b)+
     &                              Hads_l(5)*gammagh(5,a,b)
     &                           +
     &                              Hads_l(1)*gammahg(1,a,b)+
     &                              Hads_l(2)*gammahg(2,a,b)+
     &                              Hads_l(3)*gammahg(3,a,b)+
     &                              Hads_l(4)*gammahg(4,a,b)+
     &                              Hads_l(5)*gammahg(5,a,b)
     &                           +                  
     &                              A_l(1)*gammagg(1,a,b)  +
     &                              A_l(2)*gammagg(2,a,b)  +
     &                              A_l(3)*gammagg(3,a,b)  +
     &                              A_l(4)*gammagg(4,a,b)  +
     &                              A_l(5)*gammagg(5,a,b)
     &                           +
     &                              A_l(1)*gammahh(1,a,b)  +
     &                              A_l(2)*gammahh(2,a,b)  +
     &                              A_l(3)*gammahh(3,a,b)  +
     &                              A_l(4)*gammahh(4,a,b)  +
     &                              A_l(5)*gammahh(5,a,b)   
     &                           +
     &                              A_l(1)*gammagh(1,a,b)  +
     &                              A_l(2)*gammagh(2,a,b)  +
     &                              A_l(3)*gammagh(3,a,b)  +
     &                              A_l(4)*gammagh(4,a,b)  +
     &                              A_l(5)*gammagh(5,a,b)
     &                           +
     &                              A_l(1)*gammahg(1,a,b)  +
     &                              A_l(2)*gammahg(2,a,b)  +
     &                              A_l(3)*gammahg(3,a,b)  +
     &                              A_l(4)*gammahg(4,a,b)  +
     &                              A_l(5)*gammahg(5,a,b)
     &                                ) 

                  term7(a,b)=        -(
     &                              gammahh(1,1,b)*gammahh(1,1,a)+
     &                              gammahh(1,2,b)*gammahh(2,1,a)+
     &                              gammahh(1,3,b)*gammahh(3,1,a)+
     &                              gammahh(1,4,b)*gammahh(4,1,a)+
     &                              gammahh(1,5,b)*gammahh(5,1,a)+
     &                              gammahh(2,1,b)*gammahh(1,2,a)+
     &                              gammahh(2,2,b)*gammahh(2,2,a)+
     &                              gammahh(2,3,b)*gammahh(3,2,a)+
     &                              gammahh(2,4,b)*gammahh(4,2,a)+
     &                              gammahh(2,5,b)*gammahh(5,2,a)+
     &                              gammahh(3,1,b)*gammahh(1,3,a)+
     &                              gammahh(3,2,b)*gammahh(2,3,a)+
     &                              gammahh(3,3,b)*gammahh(3,3,a)+
     &                              gammahh(3,4,b)*gammahh(4,3,a)+
     &                              gammahh(3,5,b)*gammahh(5,3,a)+
     &                              gammahh(4,1,b)*gammahh(1,4,a)+
     &                              gammahh(4,2,b)*gammahh(2,4,a)+
     &                              gammahh(4,3,b)*gammahh(3,4,a)+
     &                              gammahh(4,4,b)*gammahh(4,4,a)+
     &                              gammahh(4,5,b)*gammahh(5,4,a)+
     &                              gammahh(5,1,b)*gammahh(1,5,a)+
     &                              gammahh(5,2,b)*gammahh(2,5,a)+
     &                              gammahh(5,3,b)*gammahh(3,5,a)+
     &                              gammahh(5,4,b)*gammahh(4,5,a)+
     &                              gammahh(5,5,b)*gammahh(5,5,a)
     &                           +
     &                              gammagg(1,1,b)*gammahh(1,1,a)+
     &                              gammagg(1,2,b)*gammahh(2,1,a)+
     &                              gammagg(1,3,b)*gammahh(3,1,a)+
     &                              gammagg(1,4,b)*gammahh(4,1,a)+
     &                              gammagg(1,5,b)*gammahh(5,1,a)+
     &                              gammagg(2,1,b)*gammahh(1,2,a)+
     &                              gammagg(2,2,b)*gammahh(2,2,a)+
     &                              gammagg(2,3,b)*gammahh(3,2,a)+
     &                              gammagg(2,4,b)*gammahh(4,2,a)+
     &                              gammagg(2,5,b)*gammahh(5,2,a)+
     &                              gammagg(3,1,b)*gammahh(1,3,a)+
     &                              gammagg(3,2,b)*gammahh(2,3,a)+
     &                              gammagg(3,3,b)*gammahh(3,3,a)+
     &                              gammagg(3,4,b)*gammahh(4,3,a)+
     &                              gammagg(3,5,b)*gammahh(5,3,a)+
     &                              gammagg(4,1,b)*gammahh(1,4,a)+
     &                              gammagg(4,2,b)*gammahh(2,4,a)+
     &                              gammagg(4,3,b)*gammahh(3,4,a)+
     &                              gammagg(4,4,b)*gammahh(4,4,a)+
     &                              gammagg(4,5,b)*gammahh(5,4,a)+
     &                              gammagg(5,1,b)*gammahh(1,5,a)+
     &                              gammagg(5,2,b)*gammahh(2,5,a)+
     &                              gammagg(5,3,b)*gammahh(3,5,a)+
     &                              gammagg(5,4,b)*gammahh(4,5,a)+
     &                              gammagg(5,5,b)*gammahh(5,5,a)
     &                           +
     &                              gammagg(1,1,b)*gammagh(1,1,a)+
     &                              gammagg(1,2,b)*gammagh(2,1,a)+
     &                              gammagg(1,3,b)*gammagh(3,1,a)+
     &                              gammagg(1,4,b)*gammagh(4,1,a)+
     &                              gammagg(1,5,b)*gammagh(5,1,a)+
     &                              gammagg(2,1,b)*gammagh(1,2,a)+
     &                              gammagg(2,2,b)*gammagh(2,2,a)+
     &                              gammagg(2,3,b)*gammagh(3,2,a)+
     &                              gammagg(2,4,b)*gammagh(4,2,a)+
     &                              gammagg(2,5,b)*gammagh(5,2,a)+
     &                              gammagg(3,1,b)*gammagh(1,3,a)+
     &                              gammagg(3,2,b)*gammagh(2,3,a)+
     &                              gammagg(3,3,b)*gammagh(3,3,a)+
     &                              gammagg(3,4,b)*gammagh(4,3,a)+
     &                              gammagg(3,5,b)*gammagh(5,3,a)+
     &                              gammagg(4,1,b)*gammagh(1,4,a)+
     &                              gammagg(4,2,b)*gammagh(2,4,a)+
     &                              gammagg(4,3,b)*gammagh(3,4,a)+
     &                              gammagg(4,4,b)*gammagh(4,4,a)+
     &                              gammagg(4,5,b)*gammagh(5,4,a)+
     &                              gammagg(5,1,b)*gammagh(1,5,a)+
     &                              gammagg(5,2,b)*gammagh(2,5,a)+
     &                              gammagg(5,3,b)*gammagh(3,5,a)+
     &                              gammagg(5,4,b)*gammagh(4,5,a)+
     &                              gammagg(5,5,b)*gammagh(5,5,a)
     &                           +
     &                              gammagg(1,1,b)*gammahg(1,1,a)+
     &                              gammagg(1,2,b)*gammahg(2,1,a)+
     &                              gammagg(1,3,b)*gammahg(3,1,a)+
     &                              gammagg(1,4,b)*gammahg(4,1,a)+
     &                              gammagg(1,5,b)*gammahg(5,1,a)+
     &                              gammagg(2,1,b)*gammahg(1,2,a)+
     &                              gammagg(2,2,b)*gammahg(2,2,a)+
     &                              gammagg(2,3,b)*gammahg(3,2,a)+
     &                              gammagg(2,4,b)*gammahg(4,2,a)+
     &                              gammagg(2,5,b)*gammahg(5,2,a)+
     &                              gammagg(3,1,b)*gammahg(1,3,a)+
     &                              gammagg(3,2,b)*gammahg(2,3,a)+
     &                              gammagg(3,3,b)*gammahg(3,3,a)+
     &                              gammagg(3,4,b)*gammahg(4,3,a)+
     &                              gammagg(3,5,b)*gammahg(5,3,a)+
     &                              gammagg(4,1,b)*gammahg(1,4,a)+
     &                              gammagg(4,2,b)*gammahg(2,4,a)+
     &                              gammagg(4,3,b)*gammahg(3,4,a)+
     &                              gammagg(4,4,b)*gammahg(4,4,a)+
     &                              gammagg(4,5,b)*gammahg(5,4,a)+
     &                              gammagg(5,1,b)*gammahg(1,5,a)+
     &                              gammagg(5,2,b)*gammahg(2,5,a)+
     &                              gammagg(5,3,b)*gammahg(3,5,a)+
     &                              gammagg(5,4,b)*gammahg(4,5,a)+
     &                              gammagg(5,5,b)*gammahg(5,5,a)
     &                           +
     &                              gammahh(1,1,b)*gammagg(1,1,a)+
     &                              gammahh(1,2,b)*gammagg(2,1,a)+
     &                              gammahh(1,3,b)*gammagg(3,1,a)+
     &                              gammahh(1,4,b)*gammagg(4,1,a)+
     &                              gammahh(1,5,b)*gammagg(5,1,a)+
     &                              gammahh(2,1,b)*gammagg(1,2,a)+
     &                              gammahh(2,2,b)*gammagg(2,2,a)+
     &                              gammahh(2,3,b)*gammagg(3,2,a)+
     &                              gammahh(2,4,b)*gammagg(4,2,a)+
     &                              gammahh(2,5,b)*gammagg(5,2,a)+
     &                              gammahh(3,1,b)*gammagg(1,3,a)+
     &                              gammahh(3,2,b)*gammagg(2,3,a)+
     &                              gammahh(3,3,b)*gammagg(3,3,a)+
     &                              gammahh(3,4,b)*gammagg(4,3,a)+
     &                              gammahh(3,5,b)*gammagg(5,3,a)+
     &                              gammahh(4,1,b)*gammagg(1,4,a)+
     &                              gammahh(4,2,b)*gammagg(2,4,a)+
     &                              gammahh(4,3,b)*gammagg(3,4,a)+
     &                              gammahh(4,4,b)*gammagg(4,4,a)+
     &                              gammahh(4,5,b)*gammagg(5,4,a)+
     &                              gammahh(5,1,b)*gammagg(1,5,a)+
     &                              gammahh(5,2,b)*gammagg(2,5,a)+
     &                              gammahh(5,3,b)*gammagg(3,5,a)+
     &                              gammahh(5,4,b)*gammagg(4,5,a)+
     &                              gammahh(5,5,b)*gammagg(5,5,a)
     &                           +
     &                              gammahh(1,1,b)*gammagh(1,1,a)+
     &                              gammahh(1,2,b)*gammagh(2,1,a)+
     &                              gammahh(1,3,b)*gammagh(3,1,a)+
     &                              gammahh(1,4,b)*gammagh(4,1,a)+
     &                              gammahh(1,5,b)*gammagh(5,1,a)+
     &                              gammahh(2,1,b)*gammagh(1,2,a)+
     &                              gammahh(2,2,b)*gammagh(2,2,a)+
     &                              gammahh(2,3,b)*gammagh(3,2,a)+
     &                              gammahh(2,4,b)*gammagh(4,2,a)+
     &                              gammahh(2,5,b)*gammagh(5,2,a)+
     &                              gammahh(3,1,b)*gammagh(1,3,a)+
     &                              gammahh(3,2,b)*gammagh(2,3,a)+
     &                              gammahh(3,3,b)*gammagh(3,3,a)+
     &                              gammahh(3,4,b)*gammagh(4,3,a)+
     &                              gammahh(3,5,b)*gammagh(5,3,a)+
     &                              gammahh(4,1,b)*gammagh(1,4,a)+
     &                              gammahh(4,2,b)*gammagh(2,4,a)+
     &                              gammahh(4,3,b)*gammagh(3,4,a)+
     &                              gammahh(4,4,b)*gammagh(4,4,a)+
     &                              gammahh(4,5,b)*gammagh(5,4,a)+
     &                              gammahh(5,1,b)*gammagh(1,5,a)+
     &                              gammahh(5,2,b)*gammagh(2,5,a)+
     &                              gammahh(5,3,b)*gammagh(3,5,a)+
     &                              gammahh(5,4,b)*gammagh(4,5,a)+
     &                              gammahh(5,5,b)*gammagh(5,5,a)
     &                           +
     &                              gammahh(1,1,b)*gammahg(1,1,a)+
     &                              gammahh(1,2,b)*gammahg(2,1,a)+
     &                              gammahh(1,3,b)*gammahg(3,1,a)+
     &                              gammahh(1,4,b)*gammahg(4,1,a)+
     &                              gammahh(1,5,b)*gammahg(5,1,a)+
     &                              gammahh(2,1,b)*gammahg(1,2,a)+
     &                              gammahh(2,2,b)*gammahg(2,2,a)+
     &                              gammahh(2,3,b)*gammahg(3,2,a)+
     &                              gammahh(2,4,b)*gammahg(4,2,a)+
     &                              gammahh(2,5,b)*gammahg(5,2,a)+
     &                              gammahh(3,1,b)*gammahg(1,3,a)+
     &                              gammahh(3,2,b)*gammahg(2,3,a)+
     &                              gammahh(3,3,b)*gammahg(3,3,a)+
     &                              gammahh(3,4,b)*gammahg(4,3,a)+
     &                              gammahh(3,5,b)*gammahg(5,3,a)+
     &                              gammahh(4,1,b)*gammahg(1,4,a)+
     &                              gammahh(4,2,b)*gammahg(2,4,a)+
     &                              gammahh(4,3,b)*gammahg(3,4,a)+
     &                              gammahh(4,4,b)*gammahg(4,4,a)+
     &                              gammahh(4,5,b)*gammahg(5,4,a)+
     &                              gammahh(5,1,b)*gammahg(1,5,a)+
     &                              gammahh(5,2,b)*gammahg(2,5,a)+
     &                              gammahh(5,3,b)*gammahg(3,5,a)+
     &                              gammahh(5,4,b)*gammahg(4,5,a)+
     &                              gammahh(5,5,b)*gammahg(5,5,a)
     &                           +
     &                              gammagh(1,1,b)*gammagg(1,1,a)+
     &                              gammagh(1,2,b)*gammagg(2,1,a)+
     &                              gammagh(1,3,b)*gammagg(3,1,a)+
     &                              gammagh(1,4,b)*gammagg(4,1,a)+
     &                              gammagh(1,5,b)*gammagg(5,1,a)+
     &                              gammagh(2,1,b)*gammagg(1,2,a)+
     &                              gammagh(2,2,b)*gammagg(2,2,a)+
     &                              gammagh(2,3,b)*gammagg(3,2,a)+
     &                              gammagh(2,4,b)*gammagg(4,2,a)+
     &                              gammagh(2,5,b)*gammagg(5,2,a)+
     &                              gammagh(3,1,b)*gammagg(1,3,a)+
     &                              gammagh(3,2,b)*gammagg(2,3,a)+
     &                              gammagh(3,3,b)*gammagg(3,3,a)+
     &                              gammagh(3,4,b)*gammagg(4,3,a)+
     &                              gammagh(3,5,b)*gammagg(5,3,a)+
     &                              gammagh(4,1,b)*gammagg(1,4,a)+
     &                              gammagh(4,2,b)*gammagg(2,4,a)+
     &                              gammagh(4,3,b)*gammagg(3,4,a)+
     &                              gammagh(4,4,b)*gammagg(4,4,a)+
     &                              gammagh(4,5,b)*gammagg(5,4,a)+
     &                              gammagh(5,1,b)*gammagg(1,5,a)+
     &                              gammagh(5,2,b)*gammagg(2,5,a)+
     &                              gammagh(5,3,b)*gammagg(3,5,a)+
     &                              gammagh(5,4,b)*gammagg(4,5,a)+
     &                              gammagh(5,5,b)*gammagg(5,5,a)
     &                           +
     &                              gammagh(1,1,b)*gammahh(1,1,a)+
     &                              gammagh(1,2,b)*gammahh(2,1,a)+
     &                              gammagh(1,3,b)*gammahh(3,1,a)+
     &                              gammagh(1,4,b)*gammahh(4,1,a)+
     &                              gammagh(1,5,b)*gammahh(5,1,a)+
     &                              gammagh(2,1,b)*gammahh(1,2,a)+
     &                              gammagh(2,2,b)*gammahh(2,2,a)+
     &                              gammagh(2,3,b)*gammahh(3,2,a)+
     &                              gammagh(2,4,b)*gammahh(4,2,a)+
     &                              gammagh(2,5,b)*gammahh(5,2,a)+
     &                              gammagh(3,1,b)*gammahh(1,3,a)+
     &                              gammagh(3,2,b)*gammahh(2,3,a)+
     &                              gammagh(3,3,b)*gammahh(3,3,a)+
     &                              gammagh(3,4,b)*gammahh(4,3,a)+
     &                              gammagh(3,5,b)*gammahh(5,3,a)+
     &                              gammagh(4,1,b)*gammahh(1,4,a)+
     &                              gammagh(4,2,b)*gammahh(2,4,a)+
     &                              gammagh(4,3,b)*gammahh(3,4,a)+
     &                              gammagh(4,4,b)*gammahh(4,4,a)+
     &                              gammagh(4,5,b)*gammahh(5,4,a)+
     &                              gammagh(5,1,b)*gammahh(1,5,a)+
     &                              gammagh(5,2,b)*gammahh(2,5,a)+
     &                              gammagh(5,3,b)*gammahh(3,5,a)+
     &                              gammagh(5,4,b)*gammahh(4,5,a)+
     &                              gammagh(5,5,b)*gammahh(5,5,a)
     &                           +
     &                              gammagh(1,1,b)*gammagh(1,1,a)+
     &                              gammagh(1,2,b)*gammagh(2,1,a)+
     &                              gammagh(1,3,b)*gammagh(3,1,a)+
     &                              gammagh(1,4,b)*gammagh(4,1,a)+
     &                              gammagh(1,5,b)*gammagh(5,1,a)+
     &                              gammagh(2,1,b)*gammagh(1,2,a)+
     &                              gammagh(2,2,b)*gammagh(2,2,a)+
     &                              gammagh(2,3,b)*gammagh(3,2,a)+
     &                              gammagh(2,4,b)*gammagh(4,2,a)+
     &                              gammagh(2,5,b)*gammagh(5,2,a)+
     &                              gammagh(3,1,b)*gammagh(1,3,a)+
     &                              gammagh(3,2,b)*gammagh(2,3,a)+
     &                              gammagh(3,3,b)*gammagh(3,3,a)+
     &                              gammagh(3,4,b)*gammagh(4,3,a)+
     &                              gammagh(3,5,b)*gammagh(5,3,a)+
     &                              gammagh(4,1,b)*gammagh(1,4,a)+
     &                              gammagh(4,2,b)*gammagh(2,4,a)+
     &                              gammagh(4,3,b)*gammagh(3,4,a)+
     &                              gammagh(4,4,b)*gammagh(4,4,a)+
     &                              gammagh(4,5,b)*gammagh(5,4,a)+
     &                              gammagh(5,1,b)*gammagh(1,5,a)+
     &                              gammagh(5,2,b)*gammagh(2,5,a)+
     &                              gammagh(5,3,b)*gammagh(3,5,a)+
     &                              gammagh(5,4,b)*gammagh(4,5,a)+
     &                              gammagh(5,5,b)*gammagh(5,5,a)
     &                           +
     &                              gammagh(1,1,b)*gammahg(1,1,a)+
     &                              gammagh(1,2,b)*gammahg(2,1,a)+
     &                              gammagh(1,3,b)*gammahg(3,1,a)+
     &                              gammagh(1,4,b)*gammahg(4,1,a)+
     &                              gammagh(1,5,b)*gammahg(5,1,a)+
     &                              gammagh(2,1,b)*gammahg(1,2,a)+
     &                              gammagh(2,2,b)*gammahg(2,2,a)+
     &                              gammagh(2,3,b)*gammahg(3,2,a)+
     &                              gammagh(2,4,b)*gammahg(4,2,a)+
     &                              gammagh(2,5,b)*gammahg(5,2,a)+
     &                              gammagh(3,1,b)*gammahg(1,3,a)+
     &                              gammagh(3,2,b)*gammahg(2,3,a)+
     &                              gammagh(3,3,b)*gammahg(3,3,a)+
     &                              gammagh(3,4,b)*gammahg(4,3,a)+
     &                              gammagh(3,5,b)*gammahg(5,3,a)+
     &                              gammagh(4,1,b)*gammahg(1,4,a)+
     &                              gammagh(4,2,b)*gammahg(2,4,a)+
     &                              gammagh(4,3,b)*gammahg(3,4,a)+
     &                              gammagh(4,4,b)*gammahg(4,4,a)+
     &                              gammagh(4,5,b)*gammahg(5,4,a)+
     &                              gammagh(5,1,b)*gammahg(1,5,a)+
     &                              gammagh(5,2,b)*gammahg(2,5,a)+
     &                              gammagh(5,3,b)*gammahg(3,5,a)+
     &                              gammagh(5,4,b)*gammahg(4,5,a)+
     &                              gammagh(5,5,b)*gammahg(5,5,a)
     &                           +
     &                              gammahg(1,1,b)*gammagg(1,1,a)+
     &                              gammahg(1,2,b)*gammagg(2,1,a)+
     &                              gammahg(1,3,b)*gammagg(3,1,a)+
     &                              gammahg(1,4,b)*gammagg(4,1,a)+
     &                              gammahg(1,5,b)*gammagg(5,1,a)+
     &                              gammahg(2,1,b)*gammagg(1,2,a)+
     &                              gammahg(2,2,b)*gammagg(2,2,a)+
     &                              gammahg(2,3,b)*gammagg(3,2,a)+
     &                              gammahg(2,4,b)*gammagg(4,2,a)+
     &                              gammahg(2,5,b)*gammagg(5,2,a)+
     &                              gammahg(3,1,b)*gammagg(1,3,a)+
     &                              gammahg(3,2,b)*gammagg(2,3,a)+
     &                              gammahg(3,3,b)*gammagg(3,3,a)+
     &                              gammahg(3,4,b)*gammagg(4,3,a)+
     &                              gammahg(3,5,b)*gammagg(5,3,a)+
     &                              gammahg(4,1,b)*gammagg(1,4,a)+
     &                              gammahg(4,2,b)*gammagg(2,4,a)+
     &                              gammahg(4,3,b)*gammagg(3,4,a)+
     &                              gammahg(4,4,b)*gammagg(4,4,a)+
     &                              gammahg(4,5,b)*gammagg(5,4,a)+
     &                              gammahg(5,1,b)*gammagg(1,5,a)+
     &                              gammahg(5,2,b)*gammagg(2,5,a)+
     &                              gammahg(5,3,b)*gammagg(3,5,a)+
     &                              gammahg(5,4,b)*gammagg(4,5,a)+
     &                              gammahg(5,5,b)*gammagg(5,5,a)
     &                           +
     &                              gammahg(1,1,b)*gammahh(1,1,a)+
     &                              gammahg(1,2,b)*gammahh(2,1,a)+
     &                              gammahg(1,3,b)*gammahh(3,1,a)+
     &                              gammahg(1,4,b)*gammahh(4,1,a)+
     &                              gammahg(1,5,b)*gammahh(5,1,a)+
     &                              gammahg(2,1,b)*gammahh(1,2,a)+
     &                              gammahg(2,2,b)*gammahh(2,2,a)+
     &                              gammahg(2,3,b)*gammahh(3,2,a)+
     &                              gammahg(2,4,b)*gammahh(4,2,a)+
     &                              gammahg(2,5,b)*gammahh(5,2,a)+
     &                              gammahg(3,1,b)*gammahh(1,3,a)+
     &                              gammahg(3,2,b)*gammahh(2,3,a)+
     &                              gammahg(3,3,b)*gammahh(3,3,a)+
     &                              gammahg(3,4,b)*gammahh(4,3,a)+
     &                              gammahg(3,5,b)*gammahh(5,3,a)+
     &                              gammahg(4,1,b)*gammahh(1,4,a)+
     &                              gammahg(4,2,b)*gammahh(2,4,a)+
     &                              gammahg(4,3,b)*gammahh(3,4,a)+
     &                              gammahg(4,4,b)*gammahh(4,4,a)+
     &                              gammahg(4,5,b)*gammahh(5,4,a)+
     &                              gammahg(5,1,b)*gammahh(1,5,a)+
     &                              gammahg(5,2,b)*gammahh(2,5,a)+
     &                              gammahg(5,3,b)*gammahh(3,5,a)+
     &                              gammahg(5,4,b)*gammahh(4,5,a)+
     &                              gammahg(5,5,b)*gammahh(5,5,a)
     &                           +
     &                              gammahg(1,1,b)*gammagh(1,1,a)+
     &                              gammahg(1,2,b)*gammagh(2,1,a)+
     &                              gammahg(1,3,b)*gammagh(3,1,a)+
     &                              gammahg(1,4,b)*gammagh(4,1,a)+
     &                              gammahg(1,5,b)*gammagh(5,1,a)+
     &                              gammahg(2,1,b)*gammagh(1,2,a)+
     &                              gammahg(2,2,b)*gammagh(2,2,a)+
     &                              gammahg(2,3,b)*gammagh(3,2,a)+
     &                              gammahg(2,4,b)*gammagh(4,2,a)+
     &                              gammahg(2,5,b)*gammagh(5,2,a)+
     &                              gammahg(3,1,b)*gammagh(1,3,a)+
     &                              gammahg(3,2,b)*gammagh(2,3,a)+
     &                              gammahg(3,3,b)*gammagh(3,3,a)+
     &                              gammahg(3,4,b)*gammagh(4,3,a)+
     &                              gammahg(3,5,b)*gammagh(5,3,a)+
     &                              gammahg(4,1,b)*gammagh(1,4,a)+
     &                              gammahg(4,2,b)*gammagh(2,4,a)+
     &                              gammahg(4,3,b)*gammagh(3,4,a)+
     &                              gammahg(4,4,b)*gammagh(4,4,a)+
     &                              gammahg(4,5,b)*gammagh(5,4,a)+
     &                              gammahg(5,1,b)*gammagh(1,5,a)+
     &                              gammahg(5,2,b)*gammagh(2,5,a)+
     &                              gammahg(5,3,b)*gammagh(3,5,a)+
     &                              gammahg(5,4,b)*gammagh(4,5,a)+
     &                              gammahg(5,5,b)*gammagh(5,5,a)
     &                           +
     &                              gammahg(1,1,b)*gammahg(1,1,a)+
     &                              gammahg(1,2,b)*gammahg(2,1,a)+
     &                              gammahg(1,3,b)*gammahg(3,1,a)+
     &                              gammahg(1,4,b)*gammahg(4,1,a)+
     &                              gammahg(1,5,b)*gammahg(5,1,a)+
     &                              gammahg(2,1,b)*gammahg(1,2,a)+
     &                              gammahg(2,2,b)*gammahg(2,2,a)+
     &                              gammahg(2,3,b)*gammahg(3,2,a)+
     &                              gammahg(2,4,b)*gammahg(4,2,a)+
     &                              gammahg(2,5,b)*gammahg(5,2,a)+
     &                              gammahg(3,1,b)*gammahg(1,3,a)+
     &                              gammahg(3,2,b)*gammahg(2,3,a)+
     &                              gammahg(3,3,b)*gammahg(3,3,a)+
     &                              gammahg(3,4,b)*gammahg(4,3,a)+
     &                              gammahg(3,5,b)*gammahg(5,3,a)+
     &                              gammahg(4,1,b)*gammahg(1,4,a)+
     &                              gammahg(4,2,b)*gammahg(2,4,a)+
     &                              gammahg(4,3,b)*gammahg(3,4,a)+
     &                              gammahg(4,4,b)*gammahg(4,4,a)+
     &                              gammahg(4,5,b)*gammahg(5,4,a)+
     &                              gammahg(5,1,b)*gammahg(1,5,a)+
     &                              gammahg(5,2,b)*gammahg(2,5,a)+
     &                              gammahg(5,3,b)*gammahg(3,5,a)+
     &                              gammahg(5,4,b)*gammahg(4,5,a)+
     &                              gammahg(5,5,b)*gammahg(5,5,a)
     &                                )
 
                  term8(a,b)=-2*lambda5*h0_ll(a,b)/3

                  h0_ll_tt(a,b)=2/g0_uu(1,1)*
     &                (
     &                    term1sub(a,b)+term2(a,b)+term3(a,b)+term4(a,b)
     &                      +term5(a,b)+term6(a,b)+term7(a,b)+term8(a,b)
     &                      -8*PI*(set_ll(a,b)-tr_set*g0_ll(a,b)/3)
     &                )
                end do
              end do          

              ! 0 = g^ab phi1,ab - g^ab gamma^c_ab phi1,c 
              phi10_tt=-1/g0_uu(1,1)
     &                  *(
     &                      phi10_xx(2,2)*g0_uu(2,2)+
     &                      phi10_xx(3,3)*g0_uu(3,3)+
     &                      phi10_xx(4,4)*g0_uu(4,4)+
     &                      phi10_xx(5,5)*g0_uu(5,5)+
     &                   2*(phi10_xx(1,2)*g0_uu(1,2)+
     &                      phi10_xx(1,3)*g0_uu(1,3)+
     &                      phi10_xx(1,4)*g0_uu(1,4)+
     &                      phi10_xx(1,5)*g0_uu(1,5)+
     &                      phi10_xx(2,3)*g0_uu(2,3)+
     &                      phi10_xx(2,4)*g0_uu(2,4)+
     &                      phi10_xx(2,5)*g0_uu(2,5)+
     &                      phi10_xx(3,4)*g0_uu(3,4)+
     &                      phi10_xx(3,5)*g0_uu(3,5)+
     &                      phi10_xx(4,5)*g0_uu(4,5))
     &                  -
     &                      phi10_x(1)*( gamma_ull(1,1,1)*g0_uu(1,1)+
     &                                   gamma_ull(1,2,2)*g0_uu(2,2)+
     &                                   gamma_ull(1,3,3)*g0_uu(3,3)+
     &                                   gamma_ull(1,4,4)*g0_uu(4,4)+
     &                                   gamma_ull(1,5,5)*g0_uu(5,5)+
     &                                2*(gamma_ull(1,1,2)*g0_uu(1,2)+
     &                                   gamma_ull(1,1,3)*g0_uu(1,3)+
     &                                   gamma_ull(1,1,4)*g0_uu(1,4)+
     &                                   gamma_ull(1,1,5)*g0_uu(1,5)+
     &                                   gamma_ull(1,2,3)*g0_uu(2,3)+
     &                                   gamma_ull(1,2,4)*g0_uu(2,4)+
     &                                   gamma_ull(1,2,5)*g0_uu(2,5)+
     &                                   gamma_ull(1,3,4)*g0_uu(3,4)+
     &                                   gamma_ull(1,3,5)*g0_uu(3,5)+
     &                                   gamma_ull(1,4,5)*g0_uu(4,5)) )
     &                    -
     &                      phi10_x(2)*( gamma_ull(2,1,1)*g0_uu(1,1)+
     &                                   gamma_ull(2,2,2)*g0_uu(2,2)+
     &                                   gamma_ull(2,3,3)*g0_uu(3,3)+
     &                                   gamma_ull(2,4,4)*g0_uu(4,4)+
     &                                   gamma_ull(2,5,5)*g0_uu(5,5)+
     &                                2*(gamma_ull(2,1,2)*g0_uu(1,2)+
     &                                   gamma_ull(2,1,3)*g0_uu(1,3)+
     &                                   gamma_ull(2,1,4)*g0_uu(1,4)+
     &                                   gamma_ull(2,1,5)*g0_uu(1,5)+
     &                                   gamma_ull(2,2,3)*g0_uu(2,3)+
     &                                   gamma_ull(2,2,4)*g0_uu(2,4)+
     &                                   gamma_ull(2,2,5)*g0_uu(2,5)+
     &                                   gamma_ull(2,3,4)*g0_uu(3,4)+
     &                                   gamma_ull(2,3,5)*g0_uu(3,5)+
     &                                   gamma_ull(2,4,5)*g0_uu(4,5)) )
     &                    -
     &                      phi10_x(3)*( gamma_ull(3,1,1)*g0_uu(1,1)+
     &                                   gamma_ull(3,2,2)*g0_uu(2,2)+
     &                                   gamma_ull(3,3,3)*g0_uu(3,3)+
     &                                   gamma_ull(3,4,4)*g0_uu(4,4)+
     &                                   gamma_ull(3,5,5)*g0_uu(5,5)+
     &                                2*(gamma_ull(3,1,2)*g0_uu(1,2)+
     &                                   gamma_ull(3,1,3)*g0_uu(1,3)+
     &                                   gamma_ull(3,1,4)*g0_uu(1,4)+
     &                                   gamma_ull(3,1,5)*g0_uu(1,5)+
     &                                   gamma_ull(3,2,3)*g0_uu(2,3)+
     &                                   gamma_ull(3,2,4)*g0_uu(2,4)+
     &                                   gamma_ull(3,2,5)*g0_uu(2,5)+
     &                                   gamma_ull(3,3,4)*g0_uu(3,4)+
     &                                   gamma_ull(3,3,5)*g0_uu(3,5)+
     &                                   gamma_ull(3,4,5)*g0_uu(4,5)) )
     &                    -
     &                      phi10_x(4)*( gamma_ull(4,1,1)*g0_uu(1,1)+
     &                                   gamma_ull(4,2,2)*g0_uu(2,2)+
     &                                   gamma_ull(4,3,3)*g0_uu(3,3)+
     &                                   gamma_ull(4,4,4)*g0_uu(4,4)+
     &                                   gamma_ull(4,5,5)*g0_uu(5,5)+
     &                                2*(gamma_ull(4,1,2)*g0_uu(1,2)+
     &                                   gamma_ull(4,1,3)*g0_uu(1,3)+
     &                                   gamma_ull(4,1,4)*g0_uu(1,4)+
     &                                   gamma_ull(4,1,5)*g0_uu(1,5)+
     &                                   gamma_ull(4,2,3)*g0_uu(2,3)+
     &                                   gamma_ull(4,2,4)*g0_uu(2,4)+
     &                                   gamma_ull(4,2,5)*g0_uu(2,5)+
     &                                   gamma_ull(4,3,4)*g0_uu(3,4)+
     &                                   gamma_ull(4,3,5)*g0_uu(3,5)+
     &                                   gamma_ull(4,4,5)*g0_uu(4,5)) )
     &                    -
     &                      phi10_x(5)*( gamma_ull(5,1,1)*g0_uu(1,1)+
     &                                   gamma_ull(5,2,2)*g0_uu(2,2)+
     &                                   gamma_ull(5,3,3)*g0_uu(3,3)+
     &                                   gamma_ull(5,4,4)*g0_uu(4,4)+
     &                                   gamma_ull(5,5,5)*g0_uu(5,5)+
     &                                2*(gamma_ull(5,1,2)*g0_uu(1,2)+
     &                                   gamma_ull(5,1,3)*g0_uu(1,3)+
     &                                   gamma_ull(5,1,4)*g0_uu(1,4)+
     &                                   gamma_ull(5,1,5)*g0_uu(1,5)+
     &                                   gamma_ull(5,2,3)*g0_uu(2,3)+
     &                                   gamma_ull(5,2,4)*g0_uu(2,4)+
     &                                   gamma_ull(5,2,5)*g0_uu(2,5)+
     &                                   gamma_ull(5,3,4)*g0_uu(3,4)+
     &                                   gamma_ull(5,3,5)*g0_uu(3,5)+
     &                                   gamma_ull(5,4,5)*g0_uu(4,5)) )
     &                    -
     &                      dV1_dphi10
     &                                                                )

              if (is_nan(h0_ll_tt(1,1)).or.is_nan(h0_ll_tt(1,2))
     &        .or.is_nan(h0_ll_tt(1,3)).or.is_nan(h0_ll_tt(2,2))
     &        .or.is_nan(h0_ll_tt(2,3)).or.is_nan(h0_ll_tt(3,3))
     &        .or.is_nan(h0_ll_tt(4,4)).or.is_nan(h0_ll_tt(5,5)) ) then
                write(*,*) 'h0_ll_tt(1,1)=',h0_ll_tt(1,1)
                write(*,*) 'h0_ll_tt(1,2)=',h0_ll_tt(1,2)
                write(*,*) 'h0_ll_tt(1,3)=',h0_ll_tt(1,3)
                write(*,*) 'h0_ll_tt(2,2)=',h0_ll_tt(2,2)
                write(*,*) 'h0_ll_tt(2,3)=',h0_ll_tt(2,3)
                write(*,*) 'h0_ll_tt(3,3)=',h0_ll_tt(3,3)
                write(*,*) 'h0_ll_tt(4,4)=',h0_ll_tt(4,4)
                write(*,*) 'h0_ll_tt(5,5)=',h0_ll_tt(5,5)
                stop
              end if

              ! initial second time derivatives
              gb_tt_tt0=h0_ll_tt(1,1)/(1-rho0**2)
              gb_tx_tt0=h0_ll_tt(1,2)/(1-rho0**2)
              gb_ty_tt0=h0_ll_tt(1,3)/(1-rho0**2)
              gb_xx_tt0=h0_ll_tt(2,2)/(1-rho0**2)
              gb_xy_tt0=h0_ll_tt(2,3)/(1-rho0**2)
              gb_yy_tt0=h0_ll_tt(3,3)/(1-rho0**2)
              psi_tt0  =h0_ll_tt(4,4)/(1-rho0**2)/y0**2
              phi1_tt0 =phi10_tt/(1-rho0**2)**3

              ! initialize past time level by O(h^3) expansion
              gb_tt_nm1(i,j)=gb_tt_n(i,j) - gb_tt_t0*dt
     &                         + gb_tt_tt0*dt**2/2
              gb_tx_nm1(i,j)=gb_tx_n(i,j) - gb_tx_t0*dt
     &                         + gb_tx_tt0*dt**2/2
              gb_ty_nm1(i,j)=gb_ty_n(i,j) - gb_ty_t0*dt
     &                         + gb_ty_tt0*dt**2/2
              gb_xx_nm1(i,j)=gb_xx_n(i,j) - gb_xx_t0*dt
     &                         + gb_xx_tt0*dt**2/2
              gb_xy_nm1(i,j)=gb_xy_n(i,j) - gb_xy_t0*dt
     &                         + gb_xy_tt0*dt**2/2
              gb_yy_nm1(i,j)=gb_yy_n(i,j) - gb_yy_t0*dt
     &                         + gb_yy_tt0*dt**2/2
              psi_nm1(i,j)=psi_n(i,j) - psi_t0*dt
     &                         + psi_tt0*dt**2/2
              phi1_nm1(i,j)=phi1_n(i,j) - phi1_t0*dt
     &                         + phi1_tt0*dt**2/2
              
              Hb_t_nm1(i,j)=Hb_t_n(i,j)-Hb_t_t0*dt
              Hb_x_nm1(i,j)=Hb_x_n(i,j)-Hb_x_t0*dt
              Hb_y_nm1(i,j)=Hb_y_n(i,j)-Hb_y_t0*dt

              ! initialize future time level by O(h^3) expansion
              gb_tt_np1(i,j)=gb_tt_n(i,j) + gb_tt_t0*dt
     &                         + gb_tt_tt0*dt**2/2
              gb_tx_np1(i,j)=gb_tx_n(i,j) + gb_tx_t0*dt
     &                         + gb_tx_tt0*dt**2/2
              gb_ty_np1(i,j)=gb_ty_n(i,j) + gb_ty_t0*dt
     &                         + gb_ty_tt0*dt**2/2
              gb_xx_np1(i,j)=gb_xx_n(i,j) + gb_xx_t0*dt
     &                         + gb_xx_tt0*dt**2/2
              gb_xy_np1(i,j)=gb_xy_n(i,j) + gb_xy_t0*dt
     &                         + gb_xy_tt0*dt**2/2
              gb_yy_np1(i,j)=gb_yy_n(i,j) + gb_yy_t0*dt
     &                         + gb_yy_tt0*dt**2/2
              psi_np1(i,j)=psi_n(i,j) + psi_t0*dt
     &                         + psi_tt0*dt**2/2
              phi1_np1(i,j)=phi1_n(i,j) + phi1_t0*dt
     &                         + phi1_tt0*dt**2/2

              Hb_t_np1(i,j)=Hb_t_n(i,j)+Hb_t_t0*dt
              Hb_x_np1(i,j)=Hb_x_n(i,j)+Hb_x_t0*dt
              Hb_y_np1(i,j)=Hb_y_n(i,j)+Hb_y_t0*dt

              ! diagnostic
              tfunction(i,j)=gb_tt_tt0

            end if

          end do
        end do

        return
        end
