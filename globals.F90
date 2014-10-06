

      module const_module
      implicit none
      !
      double precision, parameter :: PI = 4.0d0 * atan(1.0d0)
      double precision, parameter :: zero = 0.0d0
      double precision, parameter :: one = 1.0d0
      
      !
      integer, parameter :: COORDSYS_CART = 0
      integer, parameter :: COORDSYS_RZ = 1
      !
      integer, parameter :: BCTYPE_PER = 0
      integer, parameter :: BCTYPE_SYM = 1
      
      end module const_module
      
      module geom_module
      implicit none
      !
      !integer, parameter :: SDIM = 2
      !
      integer, save :: coordsys
      !
      double precision, save :: Lx, Ly
      integer, save :: nx, ny
      double precision, save :: dx, dy
      !
      integer, save :: bc_xlo, bc_xhi, bc_ylo, bc_yhi
      !
      integer, save :: iumax, jumax
      integer, save :: ivmax, jvmax
      
      
      end module geom_module
      
      
      module ns_module
      implicit none
      !
      !integer, save :: coordsys
      !integer, parameter :: COORDSYS_CART = 0
      !integer, parameter :: COORDSYS_RZ = 1
      !
      
      !double precision, save :: Lx, Ly
      double precision, save :: gx, gy
      !
      double precision, save :: rho1, rho2
      double precision, save :: mu1, mu2
      double precision, save :: sigma
      double precision, save :: rro
      !
      !integer, save :: nx, ny
      !double precision, save :: dx, dy
      double precision, save :: dt
      !
      integer, save :: nstep
      integer, save :: maxiter
      double precision, save :: maxerr
      double precision, save :: beta
      
      
      
      
      integer, save :: plot_int
      
      !
      !double precision, parameter :: PI = acos(-1.0d0)
      end module ns_module
      
      module front_module
      implicit none
      !
      integer, parameter :: maxnf = 4096
      end module front_module
      
      
      module prob_module
      implicit none
      !
      integer, save :: probtype
      
      integer, parameter :: PROBTYPE_DROPLET = 0
      double precision, save :: xc, yc, rad
      
      integer, parameter :: PROBTYPE_RT_INSTAB = 1
      double precision, save :: eta, omega, theta
      
      
      end module prob_module
      







