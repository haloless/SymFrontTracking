


      subroutine surf_tension( &
        Nf, xf, yf, tx, ty, &
        fx, fy, &
        nx, ny, dx, dy, coordsys)
      use ns_module, only: sigma
      implicit none
      integer, intent(in) :: Nf
      double precision, intent(in) :: xf(Nf+2), yf(Nf+2)
      double precision, intent(out) :: tx(Nf+2), ty(Nf+2)
      integer, intent(in) :: coordsys
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: dx, dy
      double precision, intent(out) :: fx(nx+2,ny+2), fy(nx+2,ny+2)
      !
      integer :: l, i, j
      double precision :: dxf, dyf, ds
      double precision :: xx, yy
      double precision :: nfx, nfy
      double precision :: nvecx, nvecy, invr
      double precision, parameter :: small = 1.0d-6
      
      ! tangent vectors
      do l = 2, Nf+1
        !dxf = xf(l+1) - xf(l)
        !dyf = yf(l+1) - yf(l)
        dxf = 0.5d0 * (xf(l+1) - xf(l-1))
        dyf = 0.5d0 * (yf(l+1) - yf(l-1))
        ds = sqrt(dxf**2 + dyf**2)
        tx(l) = dxf / ds
        ty(l) = dyf / ds
      enddo
      !tx(Nf+2) = tx(2)
      !ty(Nf+2) = ty(2)
      
      ! distribute to the grid
      fx(:,:) = 0.0d0
      fy(:,:) = 0.0d0
      
      do l = 2, Nf
        xx = 0.5d0 * (xf(l)+xf(l+1))
        yy = 0.5d0 * (yf(l)+yf(l+1))
        
        !nfx = sigma * (tx(l)-tx(l-1))
		!nfy = sigma * (ty(l)-ty(l-1))
        nfx = sigma * (tx(l+1)-tx(l))
        nfy = sigma * (ty(l+1)-ty(l))
        
        if (coordsys == 1) then
          dxf = xf(l+1) - xf(l)
          dyf = yf(l+1) - yf(l)
          ds = sqrt(dxf**2 + dyf**2)
          nvecx = dyf / ds
          nvecy = -dxf / ds
          
          invr = abs(nvecx) / (xx+1.0e-8)
          
          nfx = nfx - sigma*nvecx*invr*ds
          nfy = nfy - sigma*nvecy*invr*ds
        endif
        
        call distrib_face(xx,yy, nfx,nfy, fx,fy, nx,ny, dx,dy)
      enddo
      
      if (.true.) then ! enforce symmetry BC
        do j = 2, ny+1
          fx(1,j) = 0.0d0
          fx(nx+1,j) = 0.0d0
        enddo
        do i = 2, nx
          fx(i,2) = fx(i,2) + fx(i,1)
          fx(i,ny+1) = fx(i,ny+1) + fx(i,ny+2)
        enddo
        
        do i = 2, nx+1
          fy(i,1) = 0.0d0
          fy(i,ny+1) = 0.0d0
        enddo
        do j = 2, ny
          fy(2,j) = fy(2,j) + fy(1,j)
          fy(nx+1,j) = fy(nx+1,j) + fy(nx+2,j)
        enddo
      endif
      
      
      return
      end subroutine surf_tension
      





