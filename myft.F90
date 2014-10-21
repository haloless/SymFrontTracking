

#ifndef SDIM
#define SDIM 2
#endif


      program myft
      !
      use const_module
      use geom_module
      use ns_module
      use front_module
      use prob_module
      use simple_module
      use sip_module
      !
      implicit none
      double precision, allocatable :: u(:,:), v(:,:)
      double precision, allocatable :: p(:,:), pp(:,:)
      double precision, allocatable :: ut(:,:), vt(:,:)
      double precision, allocatable :: usave(:,:), vsave(:,:)
      double precision, allocatable :: tmp1(:,:), tmp2(:,:)
      double precision, allocatable :: uu(:,:), vv(:,:) ! node velocity
      !
      double precision, allocatable :: fx(:,:), fy(:,:)
      !
      double precision, allocatable :: x(:), y(:), xh(:), yh(:)
      double precision, allocatable :: r(:), rh(:) ! cell/face radius
      !
      integer, allocatable :: mask(:,:)
      double precision, allocatable :: color(:,:)
      double precision, allocatable :: dens(:,:), visc(:,:)
      !
      integer :: bctype(SDIM,2)
      !
      double precision, allocatable :: xf(:), yf(:)
      double precision, allocatable :: uf(:), vf(:)
      double precision, allocatable :: tx(:), ty(:)
      
      integer :: Nf
      
      !
      integer :: istep, iplot
      integer :: icycle
      logical :: isconv
      double precision :: source
      double precision :: time
      
      
      integer :: i,j,l
      
      !
      integer :: untin
      namelist /fortin/ &
        coordsys, Lx,Ly, nx,ny, &
        bc_xlo, bc_xhi, bc_ylo, bc_yhi, &
        gx,gy, rho1,rho2, mu1,mu2, sigma, rro, &
        dt, nstep, maxiter, maxerr, beta, &
        Nf, &
        niter, sormax, slarge, urf, sor, nswp, cds_mix_factor, &
        probtype, &
        xc, yc, rad, &
        eta, omega, theta, &
        plot_int
      
      !
      integer :: untchk
      double precision :: chkvol, chkvol0
      
      ! read namelists
      untin = 9
      open(unit=untin, file='prob.fortin', form='formatted', status='old')
      read(unit=untin, nml=fortin)
      close(unit=untin)
      
      
      ! check
      if (coordsys .eq. COORDSYS_CART) then
      elseif (coordsys .eq. COORDSYS_RZ) then
      else
        print *, 'Invalid coordsys=',coordsys
        stop
      endif
      
      bctype(1,1) = bc_xlo
      bctype(1,2) = bc_xhi
      bctype(2,1) = bc_ylo
      bctype(2,2) = bc_yhi
      
      untchk = 19
      open(unit=untchk, file='chk.csv', status='replace')
      write(untchk,*) 'time,vol'

      ! setup grid
      dx = Lx / nx
      dy = Ly / ny
      
      if (bc_xhi == BCTYPE_PER) then
        iumax = nx+1
        jumax = ny+1
      else
        iumax = nx
        jumax = ny+1
      endif
      if (bc_yhi == BCTYPE_PER) then
        ivmax = nx+1
        jvmax = ny+1
      else
        ivmax = nx+1
        jvmax = ny
      endif
      
      print *, 'nx=',nx, 'ny=',ny
      print *, 'dx=',dx, 'dy=',dy
      print *, 'iumax=',iumax, 'jumax=',jumax
      print *, 'ivmax=',ivmax, 'jvmax=',jvmax
      
      ! allocate buffer
      ! staggered velocity has an extra slot for periodic BC
      allocate(u(nx+2,ny+2), v(nx+2,ny+2))      
      allocate(ut(nx+2,ny+2), vt(nx+2,ny+2))
      allocate(usave(nx+2,ny+2), vsave(nx+2,ny+2))
      allocate(p(nx+2,ny+2), pp(nx+2,ny+2))
      allocate(tmp1(nx+2,ny+2), tmp2(nx+2,ny+2))
      allocate(fx(nx+2,ny+2), fy(nx+2,ny+2))
      !
      allocate(uu(nx+1,ny+1), vv(nx+1,ny+1))
      !
      allocate(x(nx+2), y(ny+2), xh(nx+1), yh(ny+1))
      allocate(r(nx+2), rh(nx+1))
      !
      allocate(mask(nx+2,ny+2))
      allocate(color(nx+2,ny+2))
      allocate(dens(nx+2,ny+2), visc(nx+2,ny+2))
      !
      allocate(xf(maxnf+2), yf(maxnf+2))
      allocate(uf(maxnf+2), vf(maxnf+2))
      allocate(tx(maxnf+2), ty(maxnf+2))
      ! SIMPLE
      allocate(axlo(nx+2,ny+2), axhi(nx+2,ny+2))
      allocate(aylo(nx+2,ny+2), ayhi(nx+2,ny+2))
      allocate(acen(nx+2,ny+2), scen(nx+2,ny+2))
      allocate(udiag(nx+2,ny+2), vdiag(nx+2,ny+2))
      ! SIP
      allocate(lw(nx+2,ny+2), ue(nx+2,ny+2))
      allocate(ls(nx+2,ny+2), un(nx+2,ny+2))
      allocate(lpr(nx+2,ny+2), res(nx+2,ny+2))
      
      ! grid position
      do i = 1, nx+2
        x(i) = dx * (dble(i)-1.5d0)
      enddo
      do j = 1, ny+2
        y(j) = dy * (dble(j)-1.5d0)
      enddo
      do i = 1, nx+1
        xh(i) = dx * (dble(i)-1.0d0)
      enddo
      do j = 1, ny+1
        yh(j) = dy * (dble(j)-1.0d0)
      enddo
      ! radius
      if (coordsys == COORDSYS_CART) then
        r(:) = 1.0d0
        rh(:) = 1.0d0
      elseif (coordsys == COORDSYS_RZ) then
        r(:) = x(:)
        rh(:) = xh(:)
      endif
      
      
      ! initialize fluid state
      u(:,:) = 0.0d0
      v(:,:) = 0.0d0
      p(:,:) = 0.0d0
      call bcfill_umac(u, nx,ny)
      call bcfill_vmac(v, nx,ny)

      
      ! initialize density and viscosity
      mask(:,:) = 1
      color(:,:) = 0.0d0
      dens(:,:) = rho1
      visc(:,:) = mu1
      if (probtype == PROBTYPE_DROPLET) then
        do i = 2, nx+1
        do j = 2, ny+1
          if ((x(i)-xc)**2 + (y(j)-yc)**2 < rad**2) then
            color(i,j) = 1.0d0
          endif
        enddo
        enddo
      elseif (probtype == PROBTYPE_RT_INSTAB) then
        do i = 2, nx+1
        do j = 2, ny+1
          if (y(j) > 0.5d0*Ly + eta*Lx*cos(omega*x(i)/Lx+theta)) then
            color(i,j) = 1.0d0
          endif
        enddo
        enddo
      else
        print *, 'Unknown probtype=',probtype
        stop
      endif
      call bcfill_cell_homo(color, nx,ny)
      !color(1,:) = color(2,:)
      !color(nx+2,:) = color(nx+1,:)
      !color(:,1) = color(:,2)
      !color(:,ny+2) = color(:,ny+1)
      call ft_color_prop(color, dens,visc, nx,ny,dx,dy)
      
      
      ! setup front elements
      uf(:) = 0.0d0
      vf(:) = 0.0d0
      if (probtype == PROBTYPE_DROPLET) then
        do l = 2, Nf+1
          xf(l) = xc + rad*sin(PI*(l-2)/(Nf-1))
          yf(l) = yc - rad*cos(PI*(l-2)/(Nf-1))
          !xf(l) = xc - rad*sin(2.0d0*PI*(l-1)/Nf)
          !yf(l) = yc + rad*cos(2.0d0*PI*(l-1)/Nf)
        enddo
        ! two wall nodes
        xf(2) = 0.0d0
        xf(Nf+1) = 0.0d0
      else if (probtype == PROBTYPE_RT_INSTAB) then
        do l = 2, Nf+1
          xf(l) = Lx/(Nf-1) * (l-2)
          yf(l) = 0.5d0*Ly + eta*Lx*cos(omega*xf(l)/Lx+theta)
          !yf(l) = Ly*0.5d0 + Lx*0.1d0*cos(xf(l))
          !yf(l) = Ly*0.5d0 + Lx*0.1d0*cos(2.0d0*PI*xf(l)/(2.0d0*Lx))
        enddo
        ! two wall nodes
        xf(2) = 0.0d0
        xf(Nf+1) = Lx
      else
        print *, 'Unknown probtype=',probtype
        stop
      endif
      ! add two ghost nodes
      call ft_bndry_pos(Nf, xf,yf)
      
      
      ! smooth properties on initialization
      if (.true.) then
        call ft_color_grad(Nf, xf,yf, mask,color, fx,fy, nx,ny, dx,dy)
        call ft_color_smooth(fx,fy, mask,color, dens,visc, r,rh,nx,ny,dx,dy)
      endif
      
      
      ! initial output
      time = 0.0d0
      istep = 0
      iplot = 0
      call vel_node(u,v, uu,vv, nx,ny, dx,dy)
      call output_csv(iplot, time, &
        color,dens,visc, p, uu,vv, &
        x,y, xh,yh, nx,ny)
      call output_front(iplot, xf, yf, uf, vf, tx, ty, Nf)
      iplot = iplot + 1
      
      call chk_color_volume(chkvol0, color, r,rh,nx,ny,dx,dy)
      write(untchk,*) time,',',chkvol0

      
      ! main loop
      do istep = 1, nstep
        time = time + dt
        print *, 'step=',istep, 'time=',time
        
        ! save velocity
        usave(:,:) = u(:,:)
        vsave(:,:) = v(:,:)
        
		!
        call surf_tension(Nf, xf,yf, tx,ty, fx,fy, nx,ny,dx,dy,coordsys)
		
        if (.true.) then
          call chk_avg_dens(rro, dens, r,rh,nx,ny,dx,dy)
          print *, 'rro=',rro
        endif
        
        isconv = .false.
        ! SIMPLE outer loop
        do icycle=1, niter
          print *, 'step=',istep, 'cycle=',icycle
          
          call bcfill_umac(u, nx,ny)
          call bcfill_vmac(v, nx,ny)
          
          call calc_umac(ut,vt,u,v,usave,vsave, p, dens,visc, fx,fy, &
            gx,gy, r,rh, nx,ny,dx,dy,dt)
          call calc_vmac(ut,vt,u,v,usave,vsave, p, dens,visc, fx,fy, &
            gx,gy, r,rh, nx,ny,dx,dy,dt)
          !
          call bcfill_umac(ut, nx,ny)
          call bcfill_vmac(vt, nx,ny)
          !
          call calc_pres2(ut,vt, p,pp, dens, r,rh,nx,ny,dx,dy,dt)
          !
          u(:,:) = ut(:,:)
          v(:,:) = vt(:,:)
          
          print *, 'RESOR=',resor
          source = max(resor(IUMAC), resor(IVMAC), resor(IPRES))
          if (source > slarge) then
            print *,'outer iteration diverging... '
            stop
          endif
          if (source < sormax) then
            isconv = .true.
            exit
          endif
        enddo ! end outer loop
        
        if (isconv) then
          print *, 'SIMPLE converged'
        else
          print *, 'SIMPLE max iteration reached'
        endif
        
        call bcfill_umac(u, nx,ny)
        call bcfill_vmac(v, nx,ny)
        call bcfill_cell_homo(p, nx,ny)
        
        !
        call ft_adv(Nf, xf,yf, uf,vf, u,v, nx,ny, dx,dy,dt)
        call ft_recons(Nf, xf,yf, dx,dy)
        !
        call ft_color_grad(Nf, xf,yf, mask,color, fx,fy, nx,ny,dx,dy)
        call ft_color_smooth(fx,fy, mask,color, dens,visc, r,rh, nx,ny,dx,dy)
        
        if (.false.) then
          call chk_color_volume(chkvol, color, r,rh,nx,ny,dx,dy)
          call ft_volconst(Nf,xf,yf,uf,vf, nx,ny,dx,dy,dt, chkvol,chkvol0)
          call ft_color_grad(Nf, xf,yf, mask,color, fx,fy, nx,ny,dx,dy)
          call ft_color_smooth(fx,fy, mask,color, dens,visc, r,rh, nx,ny,dx,dy)
        endif

        
        
        
        if (mod(istep,plot_int)==0 .or. istep==nstep) then
          call vel_node(u,v, uu,vv, nx,ny, dx,dy)
          call output_csv(iplot, time, &
            color,dens,visc, p, uu,vv, &
            x,y, xh,yh, nx,ny)
          call output_front(iplot, xf, yf, uf, vf, tx, ty, Nf);
          iplot = iplot + 1
          
          !
          call chk_color_volume(chkvol, color, r,rh,nx,ny,dx,dy)
          write(untchk,*) time,',',chkvol
          flush(untchk)
        endif
        
      enddo ! end main loop
      
      close(untchk)
      
      end program myft
      
      
      
      
      
      
      
      
      
      
      
      
      subroutine vel_node( &
        u, v, uu, vv, &
        nx, ny, dx, dy)
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: dx, dy
      double precision, intent(in) :: u(nx+2,ny+2), v(nx+2,ny+2)
      double precision, intent(out) :: uu(nx+1,ny+1), vv(nx+1,ny+1)
      !
      integer :: i, j
      
      do i = 1, nx+1
      do j = 1, ny+1
        uu(i,j) = 0.5d0 * (u(i,j)+u(i,j+1))
        vv(i,j) = 0.5d0 * (v(i,j)+v(i+1,j))
      enddo
      enddo
      
      return
      end subroutine vel_node
      
      subroutine chk_color_volume(&
        chkvol, color, &
        r, rh, nx,ny,dx,dy)
      !use ns_module, only: coordsys
      implicit none
      double precision, intent(out) :: chkvol
      integer, intent(in) :: nx,ny
      double precision, intent(in) :: dx,dy
      double precision, intent(in) :: r(nx+2), rh(nx+1)
      double precision, intent(in) :: color(nx+2,ny+2)
      !
      integer :: i,j
      
      chkvol = 0.0d0
      do i = 2, nx+1
      do j = 2, ny+1
        chkvol = chkvol + color(i,j)*r(i)*dx*dy
      enddo
      enddo
      
      return
      end subroutine chk_color_volume
      
      subroutine chk_avg_dens( &
        chkrho, dens, &
        r, rh, nx,ny,dx,dy)
      implicit none
      double precision, intent(out) :: chkrho
      integer, intent(in) :: nx,ny
      double precision, intent(in) :: dx,dy
      double precision, intent(in) :: r(nx+2), rh(nx+1)
      double precision, intent(in) :: dens(nx+2,ny+2)
      !
      integer :: i,j
      double precision :: vsum
      
      chkrho = 0.0d0
      vsum = 0.0d0
      do i = 2, nx+1
      do j = 2, ny+1
        chkrho = chkrho + dens(i,j)*r(i)*dx*dy
        vsum = vsum + r(i)*dx*dy
      enddo
      enddo
      chkrho = chkrho / vsum
      
      return
      end subroutine chk_avg_dens
      
      
      
      subroutine calc_umac( & 
        unew, vnew, &
        uold, vold, &
        usave, vsave, &
        p, &
        dens, visc, &
        fx, fy, &
        gravx, gravy, &
        r, rh, &
        nx, ny, dx, dy, dt)
      use geom_module, only: coordsys, iumax, jvmax
      use simple_module
      !
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: dx, dy, dt
      double precision, intent(in) :: r(nx+2), rh(nx+1)
      double precision, intent(out) :: unew(nx+2,ny+2), vnew(nx+2,ny+2)
      double precision, intent(in) :: uold(nx+2,ny+2), vold(nx+2,ny+2)
      double precision, intent(in) :: usave(nx+2,ny+2), vsave(nx+2,ny+2)
      double precision, intent(in) :: p(nx+2,ny+2)
      double precision, intent(in) :: fx(nx+2,ny+2), fy(nx+2,ny+2)
      double precision, intent(in) :: dens(nx+2,ny+2), visc(nx+2,ny+2)
      double precision, intent(in) :: gravx, gravy
      !
      integer, parameter :: ivar = IUMAC
      integer :: i,j
      double precision :: rhop,mup,mue,muw,mun,mus
      double precision :: se,sw,sn,ss,sp, vol
      double precision :: ue,uw,un,us, vn,vs
      double precision :: coef
      double precision :: suds, scds ! source term of UDS/CDS advection
      
      axlo(:,:) = 0.0d0
      axhi(:,:) = 0.0d0
      aylo(:,:) = 0.0d0
      ayhi(:,:) = 0.0d0
      acen(:,:) = 0.0d0
      scen(:,:) = 0.0d0
      
      
      do i = 2, iumax
      do j = 2, ny+1
        !
        rhop = 0.5d0 * (dens(i,j)+dens(i+1,j))
        mup = 0.5d0 * (visc(i,j)+visc(i+1,j))
        mue = visc(i+1,j)
        muw = visc(i,j)
        mun = 0.25d0 * (visc(i,j)+visc(i+1,j)+visc(i+1,j+1)+visc(i,j+1))
        mus = 0.25d0 * (visc(i,j)+visc(i+1,j)+visc(i+1,j-1)+visc(i,j-1))
        !
        se = dy * r(i+1)
        sw = dy * r(i)
        sn = dx * rh(i)
        ss = dx * rh(i)
        sp = dy * rh(i)
        vol = dx * dy * rh(i)
        
        ! advective velocity
        ue = 0.5d0 * (uold(i+1,j) + uold(i,j))
        uw = 0.5d0 * (uold(i,j) + uold(i-1,j))
        un = 0.5d0 * (uold(i,j+1) + uold(i,j))
        us = 0.5d0 * (uold(i,j) + uold(i,j-1))
        vn = 0.5d0 * (vold(i,j) + vold(i+1,j))
        vs = 0.5d0 * (vold(i,j-1) + vold(i+1,j-1))
        
        ! convective
        suds = 0.0d0
        scds = 0.0d0
        !
        coef = se * ue
        if (coef > 0.0d0) then
          acen(i,j) = acen(i,j) + coef
          suds = suds + coef*uold(i,j)
        else
          axhi(i,j) = axhi(i,j) + coef
          suds = suds + coef*uold(i+1,j)
        endif
        scds = scds + coef*ue
        !
        coef = sw * uw
        if (coef < 0.0d0) then
          acen(i,j) = acen(i,j) - coef
          suds = suds - coef*uold(i,j)
        else
          axlo(i,j) = axlo(i,j) - coef
          suds = suds - coef*uold(i-1,j)
        endif
        scds = scds - coef*uw
        !
        coef = sn * vn
        if (coef > 0.0d0) then
          acen(i,j) = acen(i,j) + coef
          suds = suds + coef*uold(i,j)
        else
          ayhi(i,j) = ayhi(i,j) + coef
          suds = suds + coef*uold(i,j+1)
        endif
        scds = scds + coef*un
        !
        coef = ss * vs
        if (coef < 0.0d0) then
          acen(i,j) = acen(i,j) - coef
          suds = suds - coef*uold(i,j)
        else
          aylo(i,j) = aylo(i,j) - coef
          suds = suds - coef*uold(i,j-1)
        endif
        scds = scds - coef*us
        
        ! diffusive
        coef = 1.0d0/rhop * se*mue / dx
        axhi(i,j) = axhi(i,j) - coef
        acen(i,j) = acen(i,j) + coef
        !
        coef = 1.0d0/rhop * se*muw / dx
        axlo(i,j) = axlo(i,j) - coef
        acen(i,j) = acen(i,j) + coef
        !
        coef = 1.0d0/rhop * sn*mun / dy
        ayhi(i,j) = ayhi(i,j) - coef
        acen(i,j) = acen(i,j) + coef
        !
        coef = 1.0d0/rhop * ss*mus / dy
        aylo(i,j) = aylo(i,j) - coef
        acen(i,j) = acen(i,j) + coef
        !
        if (coordsys == 1) then ! hoop stress part
          acen(i,j) = acen(i,j) + vol * 2.0d0 * mup/rhop / (rh(i)**2)
        endif
        
        ! source
        ! time-stepping part
        if (.true.) then
          acen(i,j) = acen(i,j) + vol/dt
          scen(i,j) = scen(i,j) + usave(i,j)*vol/dt
        endif
        ! pressure source
        scen(i,j) = scen(i,j) - sp/rhop * (p(i+1,j)-p(i,j))
        ! defect-correction of convection
        if (.true.) then
          scen(i,j) = scen(i,j) - cds_mix_factor*(scds-suds)
        endif
        ! cross term of viscous stress
        if (.true.) then
          scen(i,j) = scen(i,j) + 1.0d0/rhop * &
            ( se*mue*(uold(i+1,j)-uold(i,j))/dx &
            - sw*muw*(uold(i,j)-uold(i-1,j))/dx &
            + sn*mun*(vold(i+1,j)-vold(i,j))/dx &
            - ss*mus*(vold(i+1,j-1)-vold(i,j-1))/dx)
        endif
        ! jump due to gravity and surface tension
        if (.true.) then
          scen(i,j) = scen(i,j) + vol * (fx(i,j)/rhop - gravx)
        endif
      enddo
      enddo
      
      do i = 2, iumax
      do j = 2, ny+1
        ! under-relaxation
        acen(i,j) = acen(i,j) / urf(ivar)
        scen(i,j) = scen(i,j) + (1.0d0-urf(ivar))*acen(i,j)*uold(i,j)
        ! save for pressure correction
        udiag(i,j) = acen(i,j)
      enddo
      enddo
      
      ! solve
      unew(:,:) = uold(:,:)
      call sip_sol(unew,iumax,ny+1,nx,ny, &
        nswp(ivar),sor(ivar),resor(ivar))
      
      return
      end subroutine calc_umac
      
      subroutine calc_vmac( & 
        unew, vnew, &
        uold, vold, &
        usave, vsave, &
        p, &
        dens, visc, &
        fx, fy, &
        gravx, gravy, &
        r, rh, &
        nx, ny, dx, dy, dt)
      use geom_module, only: coordsys, iumax, jvmax
      use simple_module
      !
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: dx, dy, dt
      double precision, intent(in) :: r(nx+2), rh(nx+1)
      double precision, intent(out) :: unew(nx+2,ny+2), vnew(nx+2,ny+2)
      double precision, intent(in) :: uold(nx+2,ny+2), vold(nx+2,ny+2)
      double precision, intent(in) :: usave(nx+2,ny+2), vsave(nx+2,ny+2)
      double precision, intent(in) :: p(nx+2,ny+2)
      double precision, intent(in) :: dens(nx+2,ny+2), visc(nx+2,ny+2)
      double precision, intent(in) :: fx(nx+2,ny+2), fy(nx+2,ny+2)
      double precision, intent(in) :: gravx, gravy
      !
      integer, parameter :: ivar = IVMAC
      integer :: i,j
      double precision :: rhop,mup, mue,muw,mun,mus
      double precision :: se,sw,sn,ss,sp, vol
      double precision :: ve,vw,vn,vs, ue,uw
      double precision :: coef
      double precision :: suds, scds ! source term of UDS/CDS advection
      
      axlo(:,:) = 0.0d0
      axhi(:,:) = 0.0d0
      aylo(:,:) = 0.0d0
      ayhi(:,:) = 0.0d0
      acen(:,:) = 0.0d0
      scen(:,:) = 0.0d0
      
      do i = 2, nx+1
      do j = 2, jvmax
        !
        rhop = 0.5d0 * (dens(i,j)+dens(i,j+1))
        mup = 0.5d0 * (visc(i,j)+visc(i,j+1))
        mue = 0.25d0 * (visc(i,j)+visc(i+1,j)+visc(i+1,j+1)+visc(i,j+1))
        muw = 0.25d0 * (visc(i,j)+visc(i,j+1)+visc(i-1,j+1)+visc(i-1,j))
        mun = visc(i,j+1)
        mus = visc(i,j)
        !
        se = dy * rh(i)
        sw = dy * rh(i-1)
        sn = dx * r(i)
        ss = dx * r(i)
        sp = dx * r(i)
        vol = dx * dy * r(i)
        
        ! advective velocity
        ve = 0.5d0 * (vold(i+1,j) + vold(i,j))
        vw = 0.5d0 * (vold(i,j) + vold(i-1,j))
        vn = 0.5d0 * (vold(i,j+1) + vold(i,j))
        vs = 0.5d0 * (vold(i,j) + vold(i,j-1))
        ue = 0.5d0 * (uold(i,j) + uold(i,j+1))
        uw = 0.5d0 * (uold(i-1,j) + uold(i-1,j+1))
        ! convective
        suds = 0.0d0
        scds = 0.0d0
        !
        coef = se * ue
        if (coef > 0.0d0) then
          acen(i,j) = acen(i,j) + coef
          suds = suds + coef*vold(i,j)
        else
          axhi(i,j) = axhi(i,j) + coef
          suds = suds + coef*vold(i+1,j)
        endif
        scds = scds + coef*ve
        !
        coef = sw * uw
        if (coef < 0.0d0) then
          acen(i,j) = acen(i,j) - coef
          suds = suds - coef*vold(i,j)
        else
          axlo(i,j) = axlo(i,j) - coef
          suds = suds - coef*vold(i-1,j)
        endif
        scds = scds - coef*vw
        !
        coef = sn * vn
        if (coef > 0.0d0) then
          acen(i,j) = acen(i,j) + coef
          suds = suds + coef*vold(i,j)
        else
          ayhi(i,j) = ayhi(i,j) + coef
          suds = suds + coef*vold(i,j+1)
        endif
        scds = scds + coef*vn
        !
        coef = ss * vs
        if (coef < 0.0d0) then
          acen(i,j) = acen(i,j) - coef
          suds = suds - coef*vold(i,j)
        else
          aylo(i,j) = aylo(i,j) - coef
          suds = suds - coef*vold(i,j-1)
        endif
        scds = scds - coef*vs
        
        ! diffusive
        coef = 1.0d0/rhop * se*mue / dx
        axhi(i,j) = axhi(i,j) - coef
        acen(i,j) = acen(i,j) + coef
        !
        coef = 1.0d0/rhop * sw*muw / dx
        axlo(i,j) = axlo(i,j) - coef
        acen(i,j) = acen(i,j) + coef
        !
        coef = 1.0d0/rhop * sn*mun / dy
        ayhi(i,j) = ayhi(i,j) - coef
        acen(i,j) = acen(i,j) + coef
        !
        coef = 1.0d0/rhop * ss*mus / dy
        aylo(i,j) = aylo(i,j) - coef
        acen(i,j) = acen(i,j) + coef
        
        ! source
        ! time-stepping part
        if (.true.) then
          acen(i,j) = acen(i,j) + vol/dt
          scen(i,j) = scen(i,j) + vsave(i,j)*vol/dt
        endif
        ! pressure 
        scen(i,j) = scen(i,j) - sp/rhop * (p(i,j+1)-p(i,j))
        ! defect-correction
        if (.true.) then
          scen(i,j) = scen(i,j) - cds_mix_factor*(scds-suds)
        endif
        ! cross term of viscous stress
        if (.true.) then
          scen(i,j) = scen(i,j) + 1.0d0/rhop * &
            ( se*mue*(uold(i,j+1)-uold(i,j))/dy &
            - sw*muw*(uold(i-1,j+1)-uold(i-1,j))/dy &
            + sn*mun*(vold(i,j+1)-vold(i,j))/dy &
            - ss*mus*(vold(i,j)-vold(i,j-1))/dy)
        endif
        ! jump due to gravity and surface tension
        if (.true.) then
          scen(i,j) = scen(i,j) + vol*(fy(i,j)/rhop - gravy)
        endif
      enddo
      enddo
      
      do i = 2, nx+1
      do j = 2, jvmax
        ! under-relaxation
        acen(i,j) = acen(i,j) / urf(ivar)
        scen(i,j) = scen(i,j) + (1.0d0-urf(ivar))*acen(i,j)*vold(i,j)
        ! save for pressure correction
        vdiag(i,j) = acen(i,j)
      enddo
      enddo
      
      ! solve
      vnew(:,:) = vold(:,:)
      call sip_sol(vnew,nx+1,jvmax,nx,ny, &
        nswp(ivar),sor(ivar),resor(ivar))
      
      return
      end subroutine calc_vmac

      subroutine calc_pres2( &
        u, v, p, pp, &
        dens, &
        r, rh, nx, ny, dx, dy, dt)
      use const_module
      use geom_module, only: bc_xlo,bc_ylo,bc_xhi,bc_yhi, iumax,jvmax
      use simple_module
      !
      implicit none
      integer, intent(in) :: nx,ny
      double precision, intent(in) :: dx,dy,dt
      double precision, intent(inout) :: u(nx+2,ny+2), v(nx+2,ny+2)
      double precision, intent(inout) :: p(nx+2,ny+2)
      double precision, intent(out) :: pp(nx+2,ny+2) ! pressure correction
      double precision, intent(in) :: dens(nx+2,ny+2)
      double precision, intent(in) :: r(nx+2), rh(nx+1)
      !
      integer, parameter :: ivar = IPRES
      integer :: i,j
      double precision :: se,sw,sn,ss
      double precision :: rhoe,rhow,rhon,rhos
      double precision :: coef, rhssum
      double precision :: sp, rhop
      double precision :: ppref
      !
      integer :: periodic(2)
      double precision :: sol(2:nx+1,2:ny+1), rhs(2:nx+1,2:ny+1)
      double precision :: mat(0:4,2:nx+1,2:ny+1)
      
      rhs(:,:) = 0.0d0
      mat(:,:,:) = 0.0d0
      
      periodic(:) = 0
      if ((bc_xlo.eq.BCTYPE_PER) .and. (bc_xhi.eq.BCTYPE_PER)) then
        periodic(1) = nx
      endif
      if ((bc_ylo.eq.BCTYPE_PER) .and. (bc_yhi.eq.BCTYPE_PER)) then
        periodic(2) = ny
      endif

      
      do i = 2, nx+1
      do j = 2, ny+1
        se = dy * rh(i)
        sw = dy * rh(i-1)
        sn = dx * r(i)
        ss = dx * r(i)
        !
        rhoe = 0.5d0 * (dens(i,j)+dens(i+1,j))
        rhow = 0.5d0 * (dens(i,j)+dens(i-1,j))
        rhon = 0.5d0 * (dens(i,j)+dens(i,j+1))
        rhos = 0.5d0 * (dens(i,j)+dens(i,j-1))
        
        !
        if ((bc_xlo.eq.BCTYPE_PER) .or. (i.ne.2)) then
          coef = sw*sw / (rhow*udiag(i-1,j))
          mat(0,i,j) = -coef
          mat(4,i,j) = mat(4,i,j) + coef
        endif
        if ((bc_xhi.eq.BCTYPE_PER) .or. (i.ne.nx+1)) then
          coef = se*se / (rhoe*udiag(i,j))
          mat(2,i,j) = -coef
          mat(4,i,j) = mat(4,i,j) + coef
        endif
        if ((bc_ylo.eq.BCTYPE_PER) .or. (j.ne.2)) then
          coef = ss*ss / (rhos*vdiag(i,j-1))
          mat(1,i,j) = -coef
          mat(4,i,j) = mat(4,i,j) + coef
        endif
        if ((bc_yhi.eq.BCTYPE_PER) .or. (j.ne.ny+1)) then
          coef = sn*sn / (rhon*vdiag(i,j))
          mat(3,i,j) = -coef
          mat(4,i,j) = mat(4,i,j) + coef
        endif
        
        ! source
        rhs(i,j) = -(se*u(i,j) - sw*u(i-1,j) + sn*v(i,j) - ss*v(i,j-1))
      enddo
      enddo
      
      if (.true.) then
        ! check RHS sum = 0
        rhssum = 0.0d0
        do i = 2, nx+1
        do j = 2, ny+1
          rhssum = rhssum + rhs(i,j)
        enddo
        enddo
        if (abs(rhssum) > 1.0e-12) then
          print *,'sum(RHS)=',rhssum
        endif
      endif
      
      if (.true.) then
        mat(4,ipref,jpref) = 1.0d20
      endif
      
      ! guess for pressure correction
      sol(:,:) = 0.0d0
      call solve_poisson2(mat,rhs,sol, nx,ny,periodic)
      
      ! get solution
      do i = 2, nx+1
      do j = 2, ny+1
        pp(i,j) = sol(i,j)
      enddo
      enddo
      ! fill pressure boundary
      call bcfill_cell_homo(pp, nx,ny)
      
      ! correct u-component, no relaxation!
      do i = 2, iumax
      do j = 2, ny+1
        sp = dy * rh(i)
        rhop = 0.5d0 * (dens(i,j)+dens(i+1,j))
        u(i,j) = u(i,j) - sp/(rhop*udiag(i,j)) * (pp(i+1,j)-pp(i,j))
      enddo
      enddo
      ! correct v-component, no relaxation!
      do i = 2, nx+1
      do j = 2, jvmax
        sp = dx * r(i)
        rhop = 0.5d0 * (dens(i,j)+dens(i,j+1))
        v(i,j) = v(i,j) - sp/(rhop*vdiag(i,j)) * (pp(i,j+1)-pp(i,j))
      enddo
      enddo
      ! correct pressure, use relaxation
      ppref = pp(ipref,jpref)
      do i = 2, nx+1
      do j = 2, ny+1
        p(i,j) = p(i,j) + urf(ivar)*(pp(i,j)-ppref)
      enddo
      enddo
      
      return
      end subroutine calc_pres2
      
      !
      subroutine sip_sol( &
        var, &
        iend,jend, nx,ny, &
        nswp, sor, resor)
      use simple_module, only: axlo,axhi,aylo,ayhi,acen,scen
      use sip_module
      implicit none
      integer, intent(in) :: nx,ny
      integer, intent(in) :: iend,jend
      double precision, intent(inout) :: var(nx+2,ny+2)
      integer, intent(in) :: nswp
      double precision, intent(in) :: sor
      double precision, intent(out) :: resor
      !
      integer :: i,j
      double precision :: p1,p2
      double precision, parameter :: epsil = 1.0d-20
      integer :: l
      double precision :: resl, rsm
      logical :: conv
      
      
      ! initialize 
      un(:,:) = 0.0d0
      ue(:,:) = 0.0d0
      res(:,:) = 0.0d0
      
      ! generate coefficients of L and U
      do i = 2, iend
      do j = 2, jend
        lw(i,j) = axlo(i,j) / (1.0d0 + alpha*un(i-1,j))
        ls(i,j) = aylo(i,j) / (1.0d0 + alpha*ue(i,j-1))
        
        p1 = alpha * lw(i,j) * un(i-1,j)
        p2 = alpha * ls(i,j) * ue(i,j-1)
        
        lpr(i,j) = acen(i,j)+p1+p2 - lw(i,j)*ue(i-1,j) - ls(i,j)*un(i,j-1)
        lpr(i,j) = 1.0d0 / (lpr(i,j)+epsil)
        
        un(i,j) = (ayhi(i,j)-p1) * lpr(i,j)
        ue(i,j) = (axhi(i,j)-p2) * lpr(i,j)
      enddo
      enddo
      
      conv = .false.
      
      ! iteration
      do l = 1, nswp
        resl = 0.0d0
        
        ! residual and forward substitution
        do i = 2, iend
        do j = 2, jend
          res(i,j) = scen(i,j) &
          - ayhi(i,j)*var(i,j+1) - aylo(i,j)*var(i,j-1) &
          - axhi(i,j)*var(i+1,j) - axlo(i,j)*var(i-1,j) &
          - acen(i,j)*var(i,j)
          
          resl = resl + abs(res(i,j))
          
          res(i,j) = res(i,j) - ls(i,j)*res(i,j-1) - lw(i,j)*res(i-1,j)
          res(i,j) = res(i,j) * lpr(i,j)
        enddo
        enddo
        
        if (l.eq.1) then
          resor = resl
        endif
        rsm = resl / (resor+epsil)
        
        ! backward substitution and correction
        do i = iend, 2, -1
        do j = jend, 2, -1
          res(i,j) = res(i,j) - un(i,j)*res(i,j+1) - ue(i,j)*res(i+1,j)
          var(i,j) = var(i,j) + res(i,j)
        enddo
        enddo
        
        ! check convergence
        if (rsm < sor) then
          conv = .true.
          exit
        endif
      enddo ! inner iteration
      
      if (conv) then
        print *, 'SIP conv: iter=',l,'/',nswp, '; eps_res=',rsm
      else
        print *, 'SIP fail: iter=',l,'/',nswp, '; eps_res=',rsm
      endif
      
      return
      end subroutine sip_sol
      
      
      
      
      
      






