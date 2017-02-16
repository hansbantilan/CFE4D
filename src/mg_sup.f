c----------------------------------------------------------------------
c in polar coordinates t,x for x in [0,1]
c
c Relaxation, Lop, Residual routines for MG solver for zetab
c satisfying L.zetab=0, with induced MG right-hand side R giving L.zetab=R
c
c ( have defined zeta=1+(1-x0**2)**2*zetab, with the power 
c   consistent with
c   g0=gads*zeta^4=gads+gb, and thus with gb=(1-zeta^4)*gads,
c   where zetab~(1-x0**2), gb~(1-x0**2), gads~1/(1-x0**2)^2   )
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c if 
c 
c action = 1 : performs 1 RB relaxation sweep of zetab
c              (_lop, _res vars unused) ... sets norm
c action = 2 : calculates MG residual L.zetab-R --> _res vars ... sets norm
c              (_lop unused)
c action = 3 : calculates normal residual L.zetab --> _lop
c              (_res, _rhs unused)

c-----------------------------------------------------------------------      
        subroutine mg_sup(action,zetab,zetab_rhs,zetab_lop,zetab_res,
     &                    phi1,L,cmask,phys_bdy,chr,ex,x,norm,Nx)
        implicit none
        integer Nx,action
        integer phys_bdy(2)
        real*8 zetab(Nx),zetab_rhs(Nx),zetab_lop(Nx)
        real*8 zetab_res(Nx)
        real*8 phi1(Nx)
        real*8 cmask(Nx),chr(Nx)
        real*8 x(Nx),norm,ex,L

        real*8 trhoE_grad,trhoE_ptl
        real*8 ddzeta,ddzeta_Jac,grad_zetab_sq
        real*8 ddphi1,ddphi1_Jac,grad_phi1_sq

        real*8 phi10(Nx)

        integer relax,lop,residual,cdiff_method
        parameter (relax=1,lop=3,residual=2)

        real*8 x0
        real*8 dx

        real*8 zetab0

        real*8 phi10_0

        real*8 Jac,res,rhs,new_rhs
        real*8 lambda4
        integer i,pass,sum
        integer is

        logical first,do_red_black
        data first/.true./
!        data do_red_black/.false./
        data do_red_black/.true./
        save first

        real*8 PI
        parameter (PI=3.141592653589793d0)

        ! initialize fixed-size variables
        data i,pass,sum/0,0,0/
        data is/0/

        data x0/0.0/
        data dx/0.0/

        data zetab0/0.0/

        data phi10_0/0.0/

        data Jac,res,rhs,new_rhs/0.0,0.0,0.0,0.0/
        data lambda4/0.0/

        !--------------------------------------------------------------

        first=.false.

        norm=0
        sum=0

        dx=(x(2)-x(1))

        rhs=0

        ! sets AdS4D cosmological constant
        lambda4=-3/L/L

        ! manually reconstruct phi10=phi1*(1-x^2)^2 
        do i=1,Nx
          x0=x(i)
          if (phi1(i).ne.0) phi10(i)=phi1(i)*(1-x0**2)**2
          if (phi1(i).eq.0) phi10(i)=0
        end do

        ! (REGION) interior, solve L.zetab=0 Hamiltonian constraint 
        do pass=0,1
          do i=2,Nx-1
            x0=x(i)

            if (
     &          ((do_red_black.and.mod(i+pass,2).eq.0)
     &           .or.(.not.do_red_black.and.pass.eq.0)) 
     &           .and.(cmask(i).eq.1)
     &           .and.(chr(i).ne.ex)      
     &         ) then

              ! fill in zetab 
              zetab0=zetab(i)

              ! fill in phi10_0
              phi10_0=phi10(i)

              ! computes initial energy density at i, with initial data
              ! time-symmetric so phi_t=0, and free scalar so V(phi)=0
              call df_int(phi10,ddphi1,ddphi1_Jac,
     &                    grad_phi1_sq,
     &                    x,i,chr,L,ex,Nx)
              trhoE_grad=grad_phi1_sq/2
              trhoE_ptl=0

              ! computes normal residual L.zetab
              !(NOTE: the physical energy density rhoE is such that
              ! rhoE_grad=trhoE_grad*zeta^(-4)
              ! where zeta=1+(1-x0**2)**2*zetab, and
              ! rhoE_ptl=trhoE_ptl)
              call df_int(zetab,ddzeta,ddzeta_Jac,
     &                    grad_zetab_sq,
     &                    x,i,chr,L,ex,Nx)
              res=ddzeta
     &            -lambda4*(1+(1-x0**2)**2*zetab0)/4
     &            +(lambda4+8*PI*trhoE_ptl)
     &             *((1+(1-x0**2)**2*zetab0)**5)/4
     &            +8*PI*trhoE_grad*(1+(1-x0**2)**2*zetab0)/4

              ! computes MG residual L.zetab-R
              rhs=res-zetab_rhs(i)

              Jac=ddzeta_Jac
     &            -( lambda4*( (1-x0**2)**2 )/4 )
     &            +(lambda4+8*PI*trhoE_ptl)
     &             *( (1+(1-x0**2)**2*zetab0)**4 )*5/4
     &             *(1-x0**2)**2
     &            +8*PI*trhoE_grad*( (1-x0**2)**2 )/4

              ! performs action
              if (action.eq.residual) then
                zetab_res(i)=rhs
              else if (action.eq.lop) then
                zetab_lop(i)=res
              else if (action.eq.relax) then
                zetab(i)=zetab(i)-rhs/Jac
              end if

              norm=norm+rhs**2
              sum=sum+1
              
            end if
          end do
        end do

        ! (REGION) x=0 axis, solve d/dx(zetab)=0  
        if (phys_bdy(1).ne.0) then

          ! computes normal residual d/dy(zetab)
          res=zetab(1)-(4*zetab(2)-zetab(3))/3

          ! computes MG residual d/dy(zetab)-R
          rhs=res-zetab_rhs(1)

          ! computes diag. Jacobian of zetab->L.zetab transformation
          ! by differentiating L.zetab wrt. z(1) diag. entries
          Jac=1

          if (action.eq.residual) then
            zetab_res(1)=rhs
          else if (action.eq.lop) then
            zetab_lop(1)=res
          else if (action.eq.relax) then
            zetab(1)=zetab(1)-rhs/Jac                  
          end if

          norm=norm+rhs**2
          sum=sum+1

        end if

        norm=sqrt(norm/sum)

        return
        end

c-----------------------------------------------------------------------
c in polar coordinates t,x, for x in [0,1]
c
c The following initializes the rest of the metric and Hb,
c given zetab
c
c NOTE: if we ever add gb_xz,gb_yz,... must define them
c       as AMRD MG_cnst vars
c-----------------------------------------------------------------------
        subroutine init_ghb(zetab,phi1,gb_tt,gb_tx,gb_xx,psi,
     &                      L,phys_bdy,x,chr,ex,Nx,rhoa,rhob)
        implicit none
        integer Nx
        integer phys_bdy(2)
        real*8 zetab(Nx)
        real*8 phi1(Nx)
        real*8 gb_tt(Nx),gb_tx(Nx),gb_xx(Nx),psi(Nx)
        real*8 chr(Nx),ex,L
        real*8 x(Nx)
        real*8 rhoa,rhob
 
        real*8 dx

        real*8 zetab0
        
        integer i
        real*8 x0
        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 g0_tt_ads0,g0_xx_ads0,g0_psi_ads0

        real*8 xb,trans

        ! initialize fixed-size variables 
        data zetab0/0.0/

        data i/0/
        data x0/0.0/

        data g0_tt_ads0,g0_xx_ads0,g0_psi_ads0/0.0,0.0,0.0/

        data xb,trans/0.0,0.0/

        !--------------------------------------------------------------

        dx=(x(2)-x(1))

        ! initialize gbars given zetab, using full metric expression
        ! g0_ij=g0_ij_ads*zeta^2, and 
        ! g-bar expression g0_ij=g0_ij_ads+gb_ij
        ! psi expression g0_psi=g0_psi_ads+psi*x0^2
        do i=2,Nx-1
 
          x0=x(i)

          if (chr(i).ne.ex) then 

            zetab0=zetab(i)
            g0_tt_ads0 =-((1-x0**2)**2+x0**2/L**2)
     &                  /(1-x0**2)**2
            g0_xx_ads0 =(1+x0**2)**2/((1-x0**2)**2+x0**2/L**2)
     &                  /(1-x0**2)**2
            g0_psi_ads0=x0**2
     &                  /(1-x0**2)**2

            gb_tt(i)=0
            gb_tx(i)=0
            gb_xx(i)=g0_xx_ads0
     &                 *((1+(1-x0**2)**2*zetab0)**4-1)
            psi(i)=g0_psi_ads0
     &               *((1+(1-x0**2)**2*zetab0)**4-1)
     &               /(x0**2)
            xb=(rhob-x0)/(rhob-rhoa)
            if (x0.ge.rhob) then
              trans=1
            else if (x0.ge.rhoa) then
              trans=1-xb**3*(6*xb**2-15*xb+10)
            else
              trans=0
            end if
            gb_tt(i)=(gb_xx(i)+2*psi(i))*trans

          endif

        end do

        ! x=0 regularization of gbars
        call axi_reg_g(gb_tt,gb_tx,gb_xx,psi,chr,ex,L,x,Nx)

        return
        end

c----------------------------------------------------------------------
c in polar coordinates t,x, for x in [0,1]
c
c the following computes ddf, the background Laplacian on f, and 
c grad_f_sq, the background squared gradient of f; 
c evaluated at point i and at the initial time.
c----------------------------------------------------------------------
        subroutine df_int(f,ddf,ddf_Jac,grad_f_sq,
     &                    x,i,chr,L,ex,Nx)
        implicit none
        integer Nx
        integer i
        real*8 chr(Nx),ex
        real*8 x(Nx),L
        real*8 f(Nx),fdot(Nx)
        real*8 f0,f0_x(4),f0_xx(4,4),ddf,ddf_Jac,grad_f_sq

        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer a,b,c,d

        real*8 dx
        real*8 x0

        real*8 gads_uu_22,boxx_u_2

        ! initialize fixed-size variables 
        data a,b,c,d/0,0,0,0/

        data dx/0.0/
        data x0/0.0/
 
        data f0/0.0/
        data f0_x/4*0.0/
        data f0_xx/16*0.0/

        !--------------------------------------------------------------

        dx=(x(2)-x(1))

        ! sets x0 to values at i
        x0=x(i)

        ! gives values to pure AdS objects
        gads_uu_22=(1-x0**2)**2*((1-x0**2)**2+x0**2/L**2)/(1+x0**2)**2
        boxx_u_2=(1-x0**2)*
     &           (2*L**2*(1-x0**2)**3+x0**2*(3+x0**4))/
     &           (L**2*x0*(1+x0**2)**3)

        ! set first and second derivatives
        !(only the spatial part of the metric here, so don't need time derivatives)
        f0=f(i)
        f0_x(2)=(f(i+1)-f(i-1))/2/dx
        f0_xx(2,2)=(f(i+1)-2*f(i)+f(i-1))/dx/dx

        ! calculate ddzeta, background Laplacian acting on zeta
        ! and ddzeta_Jac, Jacobian of zeta->ddzeta transformation 
        ! (DDf = g^ab D_a D_b f)
        ! all in terms of zetab
        ! assuming zeta=1+(1-x0**2)**2*zetab
        ddf=
     &      gads_uu_22*
     &      ( 
     &        (-4*(1-x0**2)+8*x0**2)*f0
     &        -8*x0*(1-x0**2)*f0_x(2)
     &        +(1-x0**2)**2*f0_xx(2,2) 
     &      )
     &     +
     &      boxx_u_2*
     &      ( -4*x0*(1-x0**2)*f0+(1-x0**2)**2*f0_x(2) )
 
        ddf_Jac=
     &      gads_uu_22*
     &      (      
     &        (-4*(1-x0**2)+8*x0**2)
     &        +(1-x0**2)**2*(-2/dx/dx)   
     &      )
     &     +
     &      boxx_u_2*
     &      ( -4*x0*(1-x0**2) )

        ! calculate grad_f_sq, squared gradient of f
        !(Df^2 = g^ab D_a f D_b f)
        grad_f_sq=gads_uu_22*f0_x(2)*f0_x(2)

        return
        
        end
