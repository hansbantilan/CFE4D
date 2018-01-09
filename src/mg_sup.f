c----------------------------------------------------------------------
c Relaxation, Lop, Residual routines for MG solver for zeta
c satisfying L.zeta=0, with induced MG right-hand side R giving L.zeta=R
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c if 
c 
c action = 1 : performs 1 RB relaxation sweep of zeta
c              (_lop, _res vars unused) ... sets norm
c action = 2 : calculates MG residual L.zeta-R --> _res vars ... sets norm
c              (_lop unused)
c action = 3 : calculates normal residual L.zeta --> _lop
c              (_res, _rhs unused)

c-----------------------------------------------------------------------      
        subroutine mg_sup(action,zeta,zeta_rhs,zeta_lop,zeta_res,phi1,
     &                    L,cmask,phys_bdy,chr,ex,x,y,norm,Nx,Ny)
        implicit none
        integer Nx,Ny,action
        integer phys_bdy(4)
        real*8 zeta(Nx,Ny),zeta_rhs(Nx,Ny),zeta_lop(Nx,Ny)
        real*8 zeta_res(Nx,Ny)
        real*8 phi1(Nx,Ny)
        real*8 cmask(Nx,Ny),chr(Nx,Ny)
        real*8 x(Nx),y(Ny),norm,ex,L
        real*8 trhoE_grad,trhoE_ptl,alphasq
        real*8 zeta_x(4),ddzeta,ddzeta_Jac,grad_zeta_sq
        real*8 phi1_x(4),ddphi1,ddphi1_Jac,grad_phi1_sq

        real*8 phi10(Nx,Ny)

        real*8 g_ll(4,4)

        integer relax,lop,residual,cdiff_method
        parameter (relax=1,lop=3,residual=2)
        real*8 csr,snr,x0,y0,rho0,x2,y2,x3,y3,x4,y4
        real*8 h,dx,dy

        real*8 zeta0
        real*8 phi10_0

        real*8 Jac,res,rhs,new_rhs
        real*8 lambda4
        real*8 PI
        parameter (PI=0.3141592653589793D1)

        real*8 ddzeta_J_xx,ddzeta_J_yy

        integer i,j,pass,sum,ii,jj
        integer is

        logical first,do_red_black,extrap
        parameter (extrap=.true.) !matches extrap=.true. in df2_int
        data first/.true./
        data do_red_black/.true./
        save first

        ! initialize fixed-size variables
        data i,j,pass,sum,ii,jj/0,0,0,0,0,0/
        data is/0/

        data g_ll/16*0.0/

        data csr,snr,x0,y0,rho0/0.0,0.0,0.0,0.0,0.0/
        data x2,y2,x3,y3,x4,y4/0.0,0.0,0.0,0.0,0.0,0.0/
        data h,dx,dy/0.0,0.0,0.0/

        data zeta0/0.0/
        data phi10_0/0.0/

        data Jac,res,rhs,new_rhs/0.0,0.0,0.0,0.0/
        data lambda4/0.0/

        data ddzeta_J_xx,ddzeta_J_yy/0.0,0.0/

        !--------------------------------------------------------------

        first=.false.

        norm=0
        sum=0

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))

        h=x(2)-x(1)

        if ((abs(1-(y(2)-y(1))/h).gt.1.0d-5)) then
           write(*,*) 'error ... relax() expects dx=dy'
           stop
        end if

        rhs=0

        ! sets CFE4D cosmological constant
        lambda4=-3/L/L

        ! manually reconstruct phi10=phi1*(1-rho^2)^3 
        do i=1,Nx
          do j=1,Ny
            x0=x(i)
            y0=y(j)
            rho0=sqrt(x0**2+y0**2)
            if (phi1(i,j).ne.0) phi10(i,j)=phi1(i,j)*(1-rho0**2)**3
            if (phi1(i,j).eq.0) phi10(i,j)=0
          end do
        end do

        ! (REGION) interior, solve L.zeta=0 Hamiltonian constraint 
        do pass=0,1
          do i=2,Nx-1
            do j=2,Ny-1
              x0=x(i)
              y0=y(j)
              x2=x0*x0
              y2=y0*y0
              x3=x0*x0*x0
              y3=y0*y0*y0
              x4=x0*x0*x0*x0
              y4=y0*y0*y0*y0
              rho0=sqrt(x0**2+y0**2)
              alphasq=((1-rho0)**2+rho0**2/L**2)/((1-rho0)**2)

              if (
     &            ((do_red_black.and.mod(i+j+pass,2).eq.0)
     &             .or.(.not.do_red_black.and.pass.eq.0)) 
     &             .and.(cmask(i,j).eq.1)
     &             .and.(chr(i,j).ne.ex)      
     &           ) then

                ! fill in zeta 
                zeta0=zeta(i,j)

                ! fill in phi10_0
                phi10_0=phi10(i,j)

                ! computes initial energy density at i,j, with initial data
                ! time-symmetric so phi_t=0, and free scalar so V(phi)=0
                call df_int(phi10,phi1_x,ddphi1,ddphi1_Jac,
     &                      grad_phi1_sq,
     &                      x,y,i,j,chr,L,ex,Nx,Ny)
                trhoE_grad=grad_phi1_sq/2
                trhoE_ptl=0

                ! computes normal residual L.zeta
                !(NOTE: the physical energy density rhoE is such that
                ! rhoE_grad=trhoE_grad*zeta^(-2), rhoE_ptl=trhoE_ptl)
                call df_int(zeta,zeta_x,ddzeta,ddzeta_Jac,
     &                      grad_zeta_sq,
     &                      x,y,i,j,chr,L,ex,Nx,Ny)
                res=ddzeta
     &              -lambda4*zeta0/3
     &              +(lambda4+8*PI*trhoE_ptl)*(zeta0**3)/3
     &              +8*PI*trhoE_grad*zeta0/3
            
                ! computes MG residual L.zeta-R
                rhs=res-zeta_rhs(i,j)

                Jac=ddzeta_Jac
     &              -( lambda4/3 )
     &              +(lambda4+8*PI*trhoE_ptl)*(zeta0**2)
     &              +8*PI*trhoE_grad/3

                ! performs action
                if (action.eq.residual) then
                  zeta_res(i,j)=rhs
                else if (action.eq.lop) then
                  zeta_lop(i,j)=res
                else if (action.eq.relax) then
                  zeta(i,j)=zeta(i,j)-rhs/Jac
                end if

                norm=norm+rhs**2
                sum=sum+1
                
              end if
            end do
          end do
        end do

        ! (REGION) y=0 axis, solve d/dy(zeta)=0  
        if (phys_bdy(3).ne.0) then
          do i=2,Nx-1

            ! computes normal residual d/dy(zeta)
            res=zeta(i,1)-(4*zeta(i,2)-zeta(i,3))/3

            ! computes MG residual d/dy(zeta)-R
            rhs=res-zeta_rhs(i,1)

            ! computes diag. Jacobian of zeta->L.zeta transformation
            ! by differentiating L.zeta wrt. z(i,1) diag. entries
            Jac=1

            if (action.eq.residual) then
              zeta_res(i,1)=rhs
            else if (action.eq.lop) then
              zeta_lop(i,1)=res
            else if (action.eq.relax) then
              zeta(i,1)=zeta(i,1)-rhs/Jac                  
            end if

            norm=norm+rhs**2
            sum=sum+1

          end do
        end if

        norm=sqrt(norm/sum)

        return
        end

c-----------------------------------------------------------------------
c The following initializes the ebars
c-----------------------------------------------------------------------
        subroutine init_eb(eb_xx,eb_xy,eb_xz,eb_yy,eb_yz,eb_zz,
     &                     L,cmask,phys_bdy,x,y,chr,ex,Nx,Ny,regtype)
        implicit none
        integer Nx,Ny
        real*8 eb_xx(Nx,Ny),eb_xy(Nx,Ny),eb_xz(Nx,Ny)
        real*8 eb_yy(Nx,Ny),eb_yz(Nx,Ny),eb_zz(Nx,Ny)
        real*8 cmask(Nx,Ny),chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny)
        integer phys_bdy(4),regtype

        real*8 PI
        parameter (PI=0.3141592653589793D1)
        integer i,j
        real*8 x0,y0
        real*8 rho0

        !--------------------------------------------------------------

        ! initialize electric part of the Weyl tensor
        do i=2,Nx-1
          do j=2,Ny-1
 
            x0=x(i)
            y0=y(j)
            rho0=sqrt(x0**2+y0**2)

            if (chr(i,j).ne.ex) then 

              eb_xx(i,j)=0
              eb_xy(i,j)=0
              eb_xz(i,j)=-4*y0
              eb_yy(i,j)=0
              eb_yz(i,j)=4*x0
              eb_zz(i,j)=0 !run30

            endif

          end do
        end do

        return
        end

c-----------------------------------------------------------------------
c The following initializes the bbars
c-----------------------------------------------------------------------
        subroutine init_bb(bb_xx,bb_xy,bb_xz,bb_yy,bb_yz,bb_zz,
     &                     L,cmask,phys_bdy,x,y,chr,ex,Nx,Ny,regtype)
        implicit none
        integer Nx,Ny
        real*8 bb_xx(Nx,Ny),bb_xy(Nx,Ny),bb_xz(Nx,Ny)
        real*8 bb_yy(Nx,Ny),bb_yz(Nx,Ny),bb_zz(Nx,Ny)
        real*8 cmask(Nx,Ny),chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny)
        integer phys_bdy(4),regtype

        real*8 PI
        parameter (PI=0.3141592653589793D1)
        integer i,j
        real*8 x0,y0
        real*8 rho0

        !--------------------------------------------------------------

        ! initialize magnetic part of the Weyl tensor
        do i=2,Nx-1
          do j=2,Ny-1
 
            x0=x(i)
            y0=y(j)
            rho0=sqrt(x0**2+y0**2)

            if (chr(i,j).ne.ex) then 

              bb_xx(i,j)=0
              bb_xy(i,j)=0
              bb_xz(i,j)=0
              bb_yy(i,j)=0
              bb_yz(i,j)=0
              bb_zz(i,j)=0 !run30

            endif

          end do
        end do

        return
        end

c-----------------------------------------------------------------------
c The following initializes the rest of the metric and Hb,
c given zeta
c
c NOTE: if we ever add gb_xy,gb_xz,gb_yz, must define them
c       as AMRD MG_cnst vars
c-----------------------------------------------------------------------
        subroutine init_ghb(zeta,gb_tt,gb_tx,gb_ty,gb_xx,gb_xy,
     &                      gb_yy,psi,Hb_t,Hb_x,Hb_y,L,cmask,phys_bdy,
     &                      x,y,chr,ex,Nx,Ny,regtype,rhoa,rhob)
        implicit none
        integer Nx,Ny
        real*8 zeta(Nx,Ny)
        real*8 gb_xx(Nx,Ny),gb_tt(Nx,Ny)
        real*8 gb_yy(Nx,Ny),gb_tx(Nx,Ny),psi(Nx,Ny)
        real*8 gb_ty(Nx,Ny),gb_xy(Nx,Ny)
        real*8 Hb_t(Nx,Ny)
        real*8 Hb_x(Nx,Ny)
        real*8 Hb_y(Nx,Ny)
        real*8 cmask(Nx,Ny),chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny)
        real*8 rhoa,rhob
        integer phys_bdy(4),regtype

        real*8 zeros(Nx,Ny)

        real*8 zeta0
        
        real*8 PI
        parameter (PI=0.3141592653589793D1)
        integer i,j
        real*8 x0,y0
        real*8 rho0
        real*8 f0

        real*8 tfunction(Nx,Ny)

        real*8 rhop,trans

        real*8 g0_tt_ads0,g0_xx_ads0
        real*8 g0_xy_ads0,g0_yy_ads0,g0_psi_ads0

        !--------------------------------------------------------------

        ! initialize metric given zeta, using full metric expression
        ! g0_ij=g0_ij_ads*zeta^2, and 
        ! g-bar expression g0_ij=g0_ij_ads+gb_ij*(1-rho0^2)
        ! psi expression g0_psi=g0_psi_ads+psi*(1-rho0^2)*y0^2
        do i=2,Nx-1
          do j=2,Ny-1
 
            zeros(i,j)=0

            x0=x(i)
            y0=y(j)
            rho0=sqrt(x0**2+y0**2)

            f0=(1-rho0**2)**2+4*rho0**2/L**2

            ! set gads values
            g0_tt_ads0 =-f0
     &                 /(1-rho0**2)**2
            g0_xx_ads0 =(x0**2*(1+rho0**2)**2/f0+y0**2)
     &                 /(1-rho0**2)**2
     &                 /rho0**2
     &                 *4
            g0_xy_ads0 =((1+rho0**2)**2/f0-1)
     &                 /(1-rho0**2)**2
     &                 /rho0**2
     &                 *x0*y0
     &                 *4
            g0_yy_ads0 =(y0**2*(1+rho0**2)**2/f0+x0**2)
     &                 /(1-rho0**2)**2
     &                 /rho0**2
     &                 *4
            g0_psi_ads0=(y0**2)
     &                 /(1-rho0**2)**2
     &                 *4

            if (chr(i,j).ne.ex) then 

              zeta0=zeta(i,j)

              gb_tx(i,j)=0
              gb_ty(i,j)=0
              gb_xx(i,j)=g0_xx_ads0*(zeta0**2-1)/(1-rho0**2)
              gb_xy(i,j)=g0_xy_ads0*(zeta0**2-1)/(1-rho0**2)
              gb_yy(i,j)=g0_yy_ads0*(zeta0**2-1)/(1-rho0**2)
              psi(i,j)=g0_psi_ads0*(zeta0**2-1)/(1-rho0**2)/(y0**2)

              rhop=(rhob-rho0)/(rhob-rhoa)
              if (rho0.ge.rhob) then
                trans=1
              else if (rho0.ge.rhoa) then
                trans=1-rhop**3*(6*rhop**2-15*rhop+10)
              else
                trans=0
              end if

              !consistent with target gauge -gbtt+gbxx+gbyy+2*gbpsi=0
              gb_tt(i,j)=(gb_xx(i,j)+gb_yy(i,j)+2*psi(i,j))*trans

            endif

          end do
        end do

        return
        end

c----------------------------------------------------------------------
c the following computes ddf, the background Laplacian on f, and 
c grad_f_sq, the background squared gradient of f; 
c evaluated at point i,j and at the initial time.
c----------------------------------------------------------------------
        subroutine df_int(f0,f0_x,ddf,ddf_Jac,grad_f_sq,
     &                    x,y,i,j,chr,L,ex,Nx,Ny)
        implicit none
        integer Nx,Ny
        real*8 chr(Nx,Ny),ex
        real*8 x(Nx),y(Ny),L

        real*8 f0(Nx,Ny),f0_x(4),f0_xx(4,4),ddf,ddf_Jac,grad_f_sq

        real*8 PI
        parameter (PI=0.3141592653589793D1)

        integer i,j,i1,j1,a,b,c,d

        real*8 dx,dy,dt
        real*8 x0,y0,rho0
        real*8 x2,y2,x3,y3,x4,y4

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dt=1.0d0 ! only a placeholder for df2_int (no time derivatives needed in MG ID-solver)

        ! sets x0,y0,rho0 to values at i,j
        x0=x(i)
        y0=y(j)
        rho0=sqrt(x0**2+y0**2)

        x2=x0*x0
        y2=y0*y0
        x3=x0*x0*x0
        y3=y0*y0*y0
        x4=x0*x0*x0*x0
        y4=y0*y0*y0*y0

        ! set first and second derivatives
        !(time symmetric initial data for now, so time derivatives vanish, and
        ! S2 symmetric initial data for now, so phi, theta derivatives vanish) 
        f0_x(1)=0
        f0_x(2)=(f0(i+1,j)-f0(i-1,j))/2/dx!f_x
        f0_x(3)=(f0(i,j+1)-f0(i,j-1))/2/dy!f_y
        f0_x(4)=0

        f0_xx(1,1)=0
        f0_xx(1,2)=0
        f0_xx(1,3)=0
        f0_xx(1,4)=0
        f0_xx(2,2)=(f0(i+1,j)-2*f0(i,j)+f0(i-1,j))/dx/dx
        f0_xx(2,3)=( (f0(i+1,j+1)-f0(i+1,j-1))/2/dy
     &              -(f0(i-1,j+1)-f0(i-1,j-1))/2/dy )/2/dx
        f0_xx(2,4)=0
        f0_xx(3,3)=(f0(i,j+1)-2*f0(i,j)+f0(i,j-1))/dy/dy 
        f0_xx(3,4)=0
        f0_xx(4,4)=0

        do a=1,3
          do b=a+1,4
            f0_xx(b,a)=f0_xx(a,b)
          end do
        end do

        ! calculate ddf, background Laplacian acting on f
        ! and ddf_Jac, Jacobian of f->ddf transformation 
        !(DDf = g^ab D_a D_b f)
        ddf= 
     &        f0_x(3)*
     &        ((-1 + x0**2 - y0**2)*(-1 + x0**2 + y0**2))/(2*y0)
     &       +f0_xx(3,3)*
     &        (-1 + x0**2 + y0**2)**2/4
     &       +f0_x(2)*
     &        (-(x0*(-1 + x0**2 + y0**2)))
     &       +f0_xx(2,3)*
     &        0.0d0
     &       +f0_xx(2,2)*
     &        (-1 + x0**2 + y0**2)**2/4

        ddf_Jac=
     &            (-2/dy/dy)*
     &            (-1 + x0**2 + y0**2)**2/4
     &           +(-2/dx/dx)*
     &            (-1 + x0**2 + y0**2)**2/4

        ! calculate grad_f_sq, squared gradient of f 
        !(Df^2 = g^ab D_a f D_b f)
        grad_f_sq=
     &              f0_x(3)*f0_x(3)*
     &              (-1 + x0**2 + y0**2)**2/4
     &             +2*f0_x(3)*f0_x(2)*
     &              0.0d0
     &             +f0_x(2)*f0_x(2)*
     &              (-1 + x0**2 + y0**2)**2/4

        return
        
        end
