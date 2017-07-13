c----------------------------------------------------------------------
c this routine calculates Hb, given gb, d(gb)dt.
c----------------------------------------------------------------------
        subroutine init_hb(gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                     gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                     gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                     gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                     gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                     gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                     psi_np1,psi_n,psi_nm1,
     &                     Hb_t_n,Hb_x_n,Hb_y_n,
     &                     L,phys_bdy,x,y,dt,chr,ex,Nx,Ny,regtype)
        implicit none
        integer Nx,Ny
        integer phys_bdy(4)
        integer regtype
        real*8 dt,ex,L
        real*8 chr(Nx,Ny)
        real*8 Hb_t_n(Nx,Ny),Hb_x_n(Nx,Ny),Hb_y_n(Nx,Ny)
        real*8 gb_tt_np1(Nx,Ny),gb_tt_n(Nx,Ny),gb_tt_nm1(Nx,Ny)
        real*8 gb_tx_np1(Nx,Ny),gb_tx_n(Nx,Ny),gb_tx_nm1(Nx,Ny)
        real*8 gb_ty_np1(Nx,Ny),gb_ty_n(Nx,Ny),gb_ty_nm1(Nx,Ny)
        real*8 gb_xx_np1(Nx,Ny),gb_xx_n(Nx,Ny),gb_xx_nm1(Nx,Ny)
        real*8 gb_xy_np1(Nx,Ny),gb_xy_n(Nx,Ny),gb_xy_nm1(Nx,Ny)
        real*8 gb_yy_np1(Nx,Ny),gb_yy_n(Nx,Ny),gb_yy_nm1(Nx,Ny)
        real*8 psi_np1(Nx,Ny),psi_n(Nx,Ny),psi_nm1(Nx,Ny)
        real*8 x(Nx),y(Ny)

        integer i,j,is,ie,js,je,a,b,c,d
        real*8 dx,dy
        real*8 x0,y0
        real*8 rho0

        real*8 PI
        parameter (PI=3.141592653589793d0)

        logical ltrace,extrap_int_boundaries
        parameter (ltrace=.false.,extrap_int_boundaries=.true.)

        real*8 boxx_u(4),boxx_l(4)

        real*8 zeros(Nx,Ny)

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        ! (indices are t,x,w,y,z)
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
        data boxx_u,boxx_l/4*0.0,4*0.0/

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
        dy=y(2)-y(1)

        do i=1,Nx
         do j=1,Ny
           Hb_t_n(i,j)=0
           Hb_x_n(i,j)=0
           Hb_y_n(i,j)=0
           zeros(i,j)=0
         end do
        end do

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
            if (chr(i,j).ne.ex) then
            x0=x(i)
            y0=y(j)
            rho0=sqrt(x0**2+y0**2)

!            ! computes tensors at point i,j
!            call tensor_init(
!     &              gb_tt_np1,gb_tt_n,gb_tt_nm1,
!     &              gb_tx_np1,gb_tx_n,gb_tx_nm1,
!     &              gb_ty_np1,gb_ty_n,gb_ty_nm1,
!     &              gb_xx_np1,gb_xx_n,gb_xx_nm1,
!     &              gb_xy_np1,gb_xy_n,gb_xy_nm1,
!     &              gb_yy_np1,gb_yy_n,gb_yy_nm1,
!     &              psi_np1,psi_n,psi_nm1,
!     &              zeros,zeros,zeros,
!     &              zeros,zeros,zeros,
!     &              zeros,zeros,zeros,
!     &              zeros,zeros,zeros,
!     &              g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
!     &              gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
!     &              h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
!     &              A_l,A_l_x,Hads_l,
!     &              gamma_ull,gamma_ull_x,
!     &              riemann_ulll,ricci_ll,ricci_lu,ricci,
!     &              einstein_ll,set_ll,
!     &              phi10_x,phi10_xx,
!     &              x,y,dt,chr,L,ex,Nx,Ny,i,j)

              ! calculate boxx^c at point i,j
              ! (boxx^c = -g^ab gamma^c_ab)
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

              ! calculate boxx_a at point i,j
              ! (boxx_a = g_ab boxx^b)
              do a=1,4
                boxx_l(a)=boxx_u(1)*g0_ll(a,1)+
     &                    boxx_u(2)*g0_ll(a,2)+
     &                    boxx_u(3)*g0_ll(a,3)+
     &                    boxx_u(4)*g0_ll(a,4)
              end do

              ! here have \box{x}_a = Hads_a + e_a*H0_a, where
              ! e_t=(1-rho0^2)^2, e_x=(1-rho0^2)^2, e_y=(1-rho0^2)^2
              Hb_t_n(i,j)=(boxx_l(1)-Hads_l(1))/(1-rho0**2)**2
              Hb_x_n(i,j)=(boxx_l(2)-Hads_l(2))/(1-rho0**2)**2
              Hb_y_n(i,j)=(boxx_l(3)-Hads_l(3))/(1-rho0**2)**2

            end if
          end do
        end do

        return
        end
