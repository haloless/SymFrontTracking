


      subroutine vel_conv( &
        ustar, vstar, &
        u, v, &
        dens, visc, &
        r, rh, &
        nx, ny, dx, dy, dt)
      use geom_module, only: iumax, jvmax
      implicit none
      !
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: dx, dy, dt
      double precision, intent(in) :: r(nx+2), rh(nx+1)
      double precision, intent(out) :: ustar(nx+2,ny+2), vstar(nx+2,ny+2)
      double precision, intent(in) :: u(nx+2,ny+2), v(nx+2,ny+2)
      double precision, intent(in) :: dens(nx+2,ny+2), visc(nx+2,ny+2)
      !
      integer :: i, j
      double precision :: ue,uw,un,us, ve,vw,vn,vs
      double precision :: re,rw, rp
      
      ! U advection
      !do i = 2, nx
      do i = 2, iumax
      do j = 2, ny+1
        ue = 0.5d0 * (u(i+1,j) + u(i,j))
        uw = 0.5d0 * (u(i,j) + u(i-1,j))
        un = 0.5d0 * (u(i,j+1) + u(i,j))
        us = 0.5d0 * (u(i,j) + u(i,j-1))
        vn = 0.5d0 * (v(i,j) + v(i+1,j))
        vs = 0.5d0 * (v(i,j-1) + v(i+1,j-1))
        
        re = r(i+1)
        rw = r(i)
        rp = rh(i)
        
        ustar(i,j) = u(i,j) - dt * ( &
          1.0d0/rp*(re*ue*ue-rw*uw*uw)/dx + (un*vn-us*vs)/dy)
      enddo
      enddo
      
      ! V advection
      do i = 2, nx+1
      !do j = 2, ny
      do j = 2, jvmax
        ve = 0.5d0 * (v(i+1,j) + v(i,j))
        vw = 0.5d0 * (v(i,j) + v(i-1,j))
        vn = 0.5d0 * (v(i,j+1) + v(i,j))
        vs = 0.5d0 * (v(i,j) + v(i,j-1))
        ue = 0.5d0 * (u(i,j) + u(i,j+1))
        uw = 0.5d0 * (u(i-1,j) + u(i-1,j+1))
        
        re = rh(i)
        rw = rh(i-1)
        rp = r(i)
        
        vstar(i,j) = v(i,j) - dt * ( &
          1.0d0/rp*(re*ue*ve-rw*uw*vw)/dx + (vn*vn-vs*vs)/dy)
      enddo
      enddo
      
      return
      end subroutine vel_conv
	  
	  
	  
      subroutine vel_diff( &
        ustar, vstar, &
        u, v, &
        dens, visc, &
        r, rh, &
        nx, ny, dx, dy, dt)
      !use ns_module, only: coordsys
      use geom_module, only: coordsys, iumax, jvmax
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: dx, dy
      double precision, intent(in) :: dt
      double precision, intent(out) :: ustar(nx+2,ny+2), vstar(nx+2,ny+2)
      double precision, intent(in) :: u(nx+2,ny+2), v(nx+2,ny+2)
      double precision, intent(in) :: dens(nx+2,ny+2), visc(nx+2,ny+2)
      double precision, intent(in) :: r(nx+2), rh(nx+1)
      !
      integer :: i, j
      double precision :: due, duw, dun, dus
      double precision :: dve, dvw, dvn, dvs
      double precision :: me, mw, mn, ms
      double precision :: rhop, mp
      double precision :: re, rw, rp
	  
      ! U viscosity
      !do i = 2, nx
      do i = 2, iumax
      do j = 2, ny+1
        me = visc(i+1,j)
        mw = visc(i,j)
        due = 2.0d0 * me * (u(i+1,j)-u(i,j)) / dx
        duw = 2.0d0 * mw * (u(i,j)-u(i-1,j)) / dx
        
        mn = 0.25d0 * (visc(i,j)+visc(i+1,j)+visc(i+1,j+1)+visc(i,j+1))
        ms = 0.25d0 * (visc(i,j)+visc(i+1,j)+visc(i+1,j-1)+visc(i,j-1))
        dun = mn * ((u(i,j+1)-u(i,j))/dy + (v(i+1,j)-v(i,j))/dx)
        dus = ms * ((u(i,j)-u(i,j-1))/dy + (v(i+1,j-1)-v(i,j-1))/dx)
        
        rhop = 0.5d0 * (dens(i+1,j)+dens(i,j))
        
        re = r(i+1)
        rw = r(i)
        rp = rh(i)
        
        ustar(i,j) = ustar(i,j) + dt/rhop * &
          (1.0d0/rp*(re*due-rw*duw)/dx + (dun-dus)/dy)
        
        if (coordsys == 1) then ! hoop stress
          mp = 0.5d0 * (visc(i+1,j)+visc(i,j))
          !ustar(i,j) = ustar(i,j) - dt/rhop*2.0d0*mp/(rp**2)
          ustar(i,j) = ustar(i,j) / (1.0d0+dt/rhop*2.0d0*mp/(rp**2))
        endif
      enddo
      enddo
      
      ! V viscosity
      do i = 2, nx+1
      !do j = 2, ny
      do j = 2, jvmax
        me = 0.25d0 * (visc(i,j)+visc(i+1,j)+visc(i+1,j+1)+visc(i,j+1))
        mw = 0.25d0 * (visc(i,j)+visc(i,j+1)+visc(i-1,j+1)+visc(i-1,j))
        dve = me * ((u(i,j+1)-u(i,j))/dy + (v(i+1,j)-v(i,j))/dx)
        dvw = mw * ((u(i-1,j+1)-u(i-1,j))/dy + (v(i,j)-v(i-1,j))/dx)
        
        mn = visc(i,j+1)
        ms = visc(i,j)
        dvn = 2.0d0 * mn * (v(i,j+1)-v(i,j)) / dy
        dvs = 2.0d0 * ms * (v(i,j)-v(i,j-1)) / dy
        
        rhop = 0.5d0 * (dens(i,j+1)+dens(i,j))
        
        re = rh(i)
        rw = rh(i-1)
        rp = r(i)
        
        vstar(i,j) = vstar(i,j) + dt/rhop * &
          (1.0d0/rp*(re*dve-rw*dvw)/dx + (dvn-dvs)/dy)
      enddo
      enddo
      
      return
      end subroutine vel_diff
      
      
      subroutine vel_jump( &
        ustar, vstar, &
        fx, fy, &
        dens, visc, &
        nx, ny, dx, dy, dt)
      use geom_module, only: iumax, jvmax
      use ns_module, only: gx, gy, rro
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: dx, dy, dt
      double precision, intent(out) :: ustar(nx+2,ny+2), vstar(nx+2,ny+2)
      double precision, intent(in) :: fx(1:nx+2,1:ny+2), fy(1:nx+2,1:ny+2)
      double precision, intent(in) :: dens(1:nx+2,1:ny+2), visc(1:nx+2,1:ny+2)
      !
      integer :: i, j
      double precision :: rp
      
      !do i = 2, nx
      do i = 2, iumax
      do j = 2, ny+1
        rp = 0.5d0 * (dens(i+1,j)+dens(i,j))
        ustar(i,j) = ustar(i,j) + dt * ( &
          fx(i,j) / rp &
          - (1.0d0-rro/rp) * gx)
      enddo
      enddo
      
      do i = 2, nx+1
      !do j = 2, ny
      do j = 2, jvmax
        rp = 0.5d0 * (dens(i,j+1)+dens(i,j))
        vstar(i,j) = vstar(i,j) + dt * ( &
          fy(i,j) / rp &
          - (1.0d0-rro/rp) * gy)
      enddo
      enddo
      
      return
      end subroutine vel_jump

