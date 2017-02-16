c----------------------------------------------------------------------
c in polar coordinates t,x for x in [0,1]
c
c applies regularity conditions for the scalar field 
c----------------------------------------------------------------------

        subroutine axi_reg_phi(phi,chr,ex,L,x,Nx)
        implicit none
        integer Nx
        real*8 phi(Nx)
        real*8 chr(Nx),ex,L
        real*8 x(Nx)
        integer i
        real*8 dx

        real*8 PI
        parameter (PI=3.141592653589793d0)

        dx=x(2)-x(1)

        if (x(1).gt.dx/2) return

        if (chr(1).ne.ex .and. chr(2).ne.ex .and .chr(3).ne.ex) then
          phi(1)=(4*phi(2)-phi(3))/3
        else if (chr(1).eq.ex) then
          phi(1)=0
        else
          write(*,*) 'WARNING axi_reg_phi, dx=',dx
        end if

        return
        end

c----------------------------------------------------------------------
c in polar coordinates t,x for x in [0,1]
c
c applies regularity conditions for the metric
c----------------------------------------------------------------------

        subroutine axi_reg_g(gb_tt,gb_tx,gb_xx,psi,chr,ex,L,x,Nx)
        implicit none
        integer Nx
        integer regtype
        real*8 gb_tt(Nx),gb_tx(Nx),gb_xx(Nx),psi(Nx)
        real*8 chr(Nx),ex,L
        real*8 x(Nx)
        integer i
        real*8 dx

        real*8 PI
        parameter (PI=3.141592653589793d0)

        dx=x(2)-x(1)

        if (x(1).gt.dx/2) return

        if (chr(1).ne.ex .and. chr(2).ne.ex .and .chr(3).ne.ex) then

!          !version of axis regularity where gb_xx set at x=dx!
!          gb_tt(1)=(4*gb_tt(2)-gb_tt(3))/3                            !gbtt,x=0 at x=0 
!          psi(1)=(4*psi(2)-psi(3))/3                                  !psi,x=0 at x=0
!          gb_xx(1)=psi(1)                                             !gbxx=psi at x=0
!          gb_xx(2)=(3*gb_xx(1)+gb_xx(3))/4                            !gbxx,x=0 at x=dx
!          gb_tx(1)=0                                                  !gbtx=0 at x=0

          gb_tt(1)=(4*gb_tt(2)-gb_tt(3))/3                            !gbtt,x=0 at x=0 
          psi(1)=(4*psi(2)-psi(3))/3                                  !psi,x=0 at x=0
          gb_xx(1)=psi(1)                                             !gbxx=psi at x=0
          gb_tx(1)=0                                                  !gbtx=0 at x=0
        else if (chr(1).eq.ex) then
          gb_tt(1)=0
          gb_tx(1)=0
          gb_xx(1)=0
          psi(1)=0
        else
          write(*,*) 'WARNING axi_reg_g, dx=',dx
        end if

        return
        end

c----------------------------------------------------------------------
c in polar coordinates t,x for x in [0,1]
c
c applies regularity conditions for the source functions
c----------------------------------------------------------------------

        subroutine axi_reg_Hb(Hb_t,Hb_x,chr,ex,L,x,Nx)
        implicit none
        integer Nx
        integer regtype
        real*8 Hb_t(Nx),Hb_x(Nx)
        real*8 chr(Nx),ex,L
        real*8 x(Nx)
        integer i
        real*8 dx

        real*8 PI
        parameter (PI=3.141592653589793d0)

        dx=x(2)-x(1)

        if (x(1).gt.dx/2) return

        if (chr(1).ne.ex .and. chr(2).ne.ex .and .chr(3).ne.ex) then
          Hb_t(1)=(4*Hb_t(2)-Hb_t(3))/3
          Hb_x(1)=0
        else if (chr(1).eq.ex) then
          Hb_t(1)=0
          Hb_x(1)=0
        else
          write(*,*) 'WARNING axi_reg_Hb, dx=',dx
        end if

        return
        end
