c----------------------------------------------------------------------
c routine for computing independent residuals of the CFE4D system
c----------------------------------------------------------------------
        subroutine ires(efe_all_ires,
     &                  efe_tt_ires,efe_tx_ires,efe_ty_ires,
     &                  efe_xx_ires,efe_xy_ires,efe_yy_ires,
     &                  efe_psi_ires,
     &                  gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                  gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                  gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                  gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                  gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                  gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                  psi_np1,psi_n,psi_nm1,
     &                  Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                  Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                  Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                  phi1_np1,phi1_n,phi1_nm1,
     &                  x,y,dt,chr,L,ex,Nx,Ny,phys_bdy,ghost_width)
        implicit none
        integer Nx,Ny
        integer i,j
        integer phys_bdy(4),ghost_width(4)
        real*8 efe_all_ires(Nx,Ny)
        real*8 efe_tt_ires(Nx,Ny),efe_tx_ires(Nx,Ny),efe_ty_ires(Nx,Ny)
        real*8 efe_xx_ires(Nx,Ny),efe_xy_ires(Nx,Ny),efe_yy_ires(Nx,Ny)
        real*8 efe_psi_ires(Nx,Ny)
        real*8 chr(Nx,Ny),ex
        real*8 x(Nx),y(Ny),dt,L
        real*8 lambda5
        real*8 phi1_np1(Nx,Ny),phi1_n(Nx,Ny),phi1_nm1(Nx,Ny)
        real*8 gb_tt_np1(Nx,Ny),gb_tx_np1(Nx,Ny),gb_ty_np1(Nx,Ny)
        real*8 gb_xx_np1(Nx,Ny),gb_xy_np1(Nx,Ny),gb_yy_np1(Nx,Ny)
        real*8 psi_np1(Nx,Ny)
        real*8 gb_tt_n(Nx,Ny),gb_tx_n(Nx,Ny),gb_ty_n(Nx,Ny)
        real*8 gb_xx_n(Nx,Ny),gb_xy_n(Nx,Ny),gb_yy_n(Nx,Ny)
        real*8 psi_n(Nx,Ny)
        real*8 gb_tt_nm1(Nx,Ny),gb_tx_nm1(Nx,Ny),gb_ty_nm1(Nx,Ny)
        real*8 gb_xx_nm1(Nx,Ny),gb_xy_nm1(Nx,Ny),gb_yy_nm1(Nx,Ny)
        real*8 psi_nm1(Nx,Ny)

        real*8 Hb_t_np1(Nx,Ny),Hb_t_n(Nx,Ny),Hb_t_nm1(Nx,Ny)
        real*8 Hb_x_np1(Nx,Ny),Hb_x_n(Nx,Ny),Hb_x_nm1(Nx,Ny)
        real*8 Hb_y_np1(Nx,Ny),Hb_y_n(Nx,Ny),Hb_y_nm1(Nx,Ny)

        integer is,ie,js,je

        integer i1,j1,k1,a,b,c,d,e
        real*8 efe_ires(5,5)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy
        real*8 x0,y0,rho0        

        real*8 boxx_u(5),boxx_l(5) 

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,y,theta,phi)
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
        ! variables for outward null expansion in spherical symmetry
        !--------------------------------------------------------------
        real*8 n_l(5),s_l(5)
        real*8 n_u(5),s_u(5)
        real*8 n_l_x(5,5),n_u_x(5,5),s_l_x(5,5)
        real*8 f0_x(5)
        real*8 f0_xx(5,5)
        real*8 gam_uu(5,5),sig_uu(5,5)
        real*8 gam_uu_x(5,5,5)
        real*8 normsusq

        real*8 g0gamfx(5)
        real*8 nufx,nuxfx(5),gamxfxfx(5)

        real*8 theta(Nx,Ny)

        ! initialize fixed-size variables
        data i1,j1,k1,a,b,c,d,e/0,0,0,0,0,0,0,0/

        data dx,dy/0.0,0.0/
        data x0,y0,rho0/0.0,0.0,0.0/    

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

        data boxx_u,boxx_l/5*0.0,5*0.0/ 

!----------------------------------------------------------------------
        
        dx=(x(2)-x(1))
        dy=(y(2)-y(1))

        ! CFE4D cosmological constant
        !(lambda5=-(n-1)(n-2)/2/L^2) for n=5 dimensional AdS)
        lambda5=-6/L/L

        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)

        ! (MAIN LOOP) loop through spacetime points x(i),y(j)
        do i=is,ie
          do j=js,je

            if (chr(i,j).ne.ex) then

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

              ! calculates efe_ires functions at point i,j
              !(efe_ires_ab=G_ab+lambda5*g_ab-8*PI*T_ab)
              do a=1,5
                do b=a,5
                  efe_ires(a,b)=einstein_ll(a,b)+lambda5*g0_ll(a,b)
     &                                          -8*PI*set_ll(a,b)
                end do
              end do
              efe_tt_ires(i,j)=efe_ires(1,1) 
              efe_tx_ires(i,j)=efe_ires(1,2) 
              efe_ty_ires(i,j)=efe_ires(1,3) 
              efe_xx_ires(i,j)=efe_ires(2,2)
              efe_xy_ires(i,j)=efe_ires(2,3)
              efe_yy_ires(i,j)=efe_ires(3,3)
              efe_psi_ires(i,j)=efe_ires(4,4)

              ! calculate efe_all_ires function at point i,j
              efe_all_ires(i,j)=
     &        max(abs(efe_tt_ires(i,j)),abs(efe_tx_ires(i,j)),
     &            abs(efe_ty_ires(i,j)),abs(efe_xx_ires(i,j)),
     &            abs(efe_xy_ires(i,j)),abs(efe_yy_ires(i,j)),
     &            abs(efe_psi_ires(i,j)))

              x0=x(i)
              y0=y(j)
              rho0=sqrt(x0**2+y0**2)              

              ! calculate boxx^c at point i,j
              ! (boxx^c = -g^ab gamma^c_ab)
              do c=1,5
                boxx_u(c)=-( gamma_ull(c,1,1)*g0_uu(1,1)+
     &                       gamma_ull(c,2,2)*g0_uu(2,2)+
     &                       gamma_ull(c,3,3)*g0_uu(3,3)+
     &                       gamma_ull(c,4,4)*g0_uu(4,4)+
     &                       gamma_ull(c,5,5)*g0_uu(5,5)+
     &                    2*(gamma_ull(c,1,2)*g0_uu(1,2)+
     &                       gamma_ull(c,1,3)*g0_uu(1,3)+
     &                       gamma_ull(c,1,4)*g0_uu(1,4)+
     &                       gamma_ull(c,1,5)*g0_uu(1,5)+
     &                       gamma_ull(c,2,3)*g0_uu(2,3)+
     &                       gamma_ull(c,2,4)*g0_uu(2,4)+
     &                       gamma_ull(c,2,5)*g0_uu(2,5)+
     &                       gamma_ull(c,3,4)*g0_uu(3,4)+
     &                       gamma_ull(c,3,5)*g0_uu(3,5)+
     &                       gamma_ull(c,4,5)*g0_uu(4,5)) )
              end do

              ! calculate boxx_a at point i,j
              ! (boxx_a = g_ab boxx^b)
              do a=1,5
                boxx_l(a)=boxx_u(1)*g0_ll(a,1)+
     &                    boxx_u(2)*g0_ll(a,2)+
     &                    boxx_u(3)*g0_ll(a,3)+
     &                    boxx_u(4)*g0_ll(a,4)+
     &                    boxx_u(5)*g0_ll(a,5)
              end do

              ! define unit time-like vector n, normal to t=const
              ! surfaces
              n_l(1)=-1/sqrt(-g0_uu(1,1))
              do a=1,5
                n_u(a)=n_l(1)*g0_uu(a,1)+
     &                 n_l(2)*g0_uu(a,2)+
     &                 n_l(3)*g0_uu(a,3)+
     &                 n_l(4)*g0_uu(a,4)+
     &                 n_l(5)*g0_uu(a,5)
              end do
              do b=1,5
                n_l_x(1,b)=-1/2.0d0/sqrt(-g0_uu(1,1))**3*g0_uu_x(1,1,b)
              end do
              do a=1,5
                do b=1,5
                  n_u_x(a,b)=n_l_x(1,b)*g0_uu(a,1)+
     &                       n_l_x(2,b)*g0_uu(a,2)+
     &                       n_l_x(3,b)*g0_uu(a,3)+
     &                       n_l_x(4,b)*g0_uu(a,4)+
     &                       n_l_x(5,b)*g0_uu(a,5)+
     &                       n_l(1)*g0_uu_x(a,1,b)+
     &                       n_l(2)*g0_uu_x(a,2,b)+
     &                       n_l(3)*g0_uu_x(a,3,b)+
     &                       n_l(4)*g0_uu_x(a,4,b)+
     &                       n_l(5)*g0_uu_x(a,5,b)
                end do
              end do

              ! define gradients of the flow field f=r-AH_R(chi,phi) 
              ! NOTE: CHECK THESE WITH MATHEMATICA GIVEN
              ! f(x,y)=sqrt(x^2+y^2)
              f0_x(1)=0
              f0_x(2)=x0/rho0
              f0_x(3)=y0/rho0
              f0_x(4)=0
              f0_x(5)=0
              f0_xx(1,1)=0
              f0_xx(1,2)=0
              f0_xx(1,3)=0
              f0_xx(1,4)=0
              f0_xx(1,5)=0
              f0_xx(2,2)=y0**2/rho0**3
              f0_xx(2,3)=-x0*y0/rho0**3
              f0_xx(2,4)=0
              f0_xx(2,5)=0
              f0_xx(3,3)=x0**2/rho0**3
              f0_xx(3,4)=0
              f0_xx(3,5)=0
              f0_xx(4,4)=0
              f0_xx(4,5)=0
              f0_xx(5,5)=0

              do a=1,4
                do b=a+1,5
                  f0_xx(b,a)=f0_xx(a,b)
                end do
              end do

              ! define metric on codimension-1 surfaces
              do a=1,5
                do b=1,5
                  gam_uu(a,b)=g0_uu(a,b)+n_u(a)*n_u(b)
                end do
              end do
              do a=1,5
                do b=1,5
                  do c=1,5
                    gam_uu_x(a,b,c)=g0_uu_x(a,b,c)
     &                             +n_u_x(a,c)*n_u(b)
     &                             +n_u(a)*n_u_x(b,c)
                  end do
                end do
              end do

              ! define unit space-like vector s, orthogonal to n and
              ! projected gradient of the flow field f
              do a=1,5
                s_u(a)=0.0d0
                do b=1,5
                  s_u(a)=s_u(a)+gam_uu(a,b)*f0_x(b)
                end do
              end do
              normsusq=0.0d0
              do a=1,5
                do b=1,5
                  normsusq=normsusq+gam_uu(a,b)*f0_x(a)*f0_x(b)
                end do
              end do
              do a=1,5
                s_u(a)=s_u(a)/sqrt(normsusq)
              end do
              do a=1,5
                s_l(a)=s_u(1)*g0_ll(a,1)+
     &                 s_u(2)*g0_ll(a,2)+
     &                 s_u(3)*g0_ll(a,3)+
     &                 s_u(4)*g0_ll(a,4)+
     &                 s_u(5)*g0_ll(a,5)
              end do

              nufx=0
              do a=1,5
                nufx=nufx
     &              +n_u(a)*f0_x(a)
                nuxfx(a)=0
                gamxfxfx(a)=0
                do c=1,5
                  nuxfx(a)=nuxfx(a)
     &                    +n_u_x(c,a)*f0_x(c)
                  do d=1,5
                    gamxfxfx(a)=gamxfxfx(a)
     &                        +gam_uu_x(c,d,a)*f0_x(c)*f0_x(d)
     &                        +gam_uu(c,d)*f0_xx(c,a)*f0_x(d)
     &                        +gam_uu(c,d)*f0_x(c)*f0_xx(d,a)
                  end do
                end do
              end do
              do a=1,5
                do b=1,5
                  s_l_x(a,b)=
     &                     (f0_xx(a,b)+n_l_x(a,b)*nufx+n_l(a)*nuxfx(b))
     &                     /sqrt(normsusq)
     &                    -(f0_x(a)+n_l(a)*nufx)*gamxfxfx(b)
     &                     /2.0d0/sqrt(normsusq)**3
                end do
              end do

              ! define metric on codimension-2 surfaces
              do a=1,5
                do b=1,5
                  sig_uu(a,b)=g0_uu(a,b)+n_u(a)*n_u(b)-s_u(a)*s_u(b)
                end do
              end do

              ! for theta: outward null expansion
              theta(i,j)=0.0d0
              do c=1,5
                do d=1,5
                  theta(i,j)=theta(i,j)
     &                   +sig_uu(c,d)*(n_l_x(c,d)+s_l_x(c,d))
                  do e=1,5
                    theta(i,j)=theta(i,j)
     &                   -sig_uu(c,d)*gamma_ull(e,c,d)*(n_l(e)+s_l(e))
                  end do
                end do
              end do

              efe_tt_ires(i,j)=!Hads_l(1)+A_l(1)-boxx_l(1)
     &           sqrt((-(-1+rho0**2)**6*gb_xy_n(i,j)**2+
     &                (-4+(-1+rho0**2)**3*gb_xx_n(i,j))*
     &                (-4+(-1+rho0**2)**3*gb_yy_n(i,j)))*
     &                (-4+(-1+rho0**2)**3*psi_n(i,j))**2)!/(-1+rho0**2)**8*y0**4
              efe_tx_ires(i,j)=g0_ll(1,2)!Hads_l(2)+A_l(2)-boxx_l(2)
              efe_ty_ires(i,j)=g0_ll(1,3)!Hads_l(3)+A_l(3)-boxx_l(3)
              efe_xx_ires(i,j)=theta(i,j)
              efe_xy_ires(i,j)=boxx_l(2)
              efe_yy_ires(i,j)=-1/g0_uu(1,1)*(1-x(i)**2-y(j)**2)**2
              efe_psi_ires(i,j)=0
              do a=1,5
                do b=1,5
                  do c=1,5
                    do d=1,5
                      do i1=1,5
                        do j1=1,5
                          do k1=1,5
                            do e=1,5
                              efe_psi_ires(i,j)=efe_psi_ires(i,j)+
     &                                          g0_ll(a,i1)*
     &                                          g0_uu(b,j1)*
     &                                          g0_uu(c,k1)*
     &                                          g0_uu(d,e)*
     &                                          riemann_ulll(a,b,c,d)*
     &                                          riemann_ulll(i1,j1,k1,e)
                            end do
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do

            end if

          end do
        end do

        return
        end
