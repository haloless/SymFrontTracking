c-----------------------------------------------------------------------c
c ffffffffff  tttttttttttttt     cccccccccc    2222222     dddddddd     c
c ffffffffff  tttttttttttttt   cccccccccccc  c2222222222   dddddddddd   c
c ffff             tttt       cccccc         2222    2222  dddd   dddd  c
c ffff             tttt       cccc                   2222  dddd    dddd c
c ffffffff         tttt       cccc                 22222   dddd    dddd c
c ffffffff         tttt       cccc               22222     dddd    dddd c
c ffff             tttt       cccc             22222       dddd    dddd c
c ffff             tttt       cccccc         22222         dddd   dddd  c
c ffff             tttt        cccccccccccc  222222222222  dddddddddd   c
c ffff             tttt          cccccccccc  222222222222  dddddddd     c
c-----------------------------------------------------------------------c
c                                                                       c
c  A two-dimensional, Front Tracking Code based on the Immerged         c
c  Boundary Technique The code is based on the method developed by      c
c  S.O. Unverdi in his Thesis and is discussed in Unverdi and           c
c  Tryggvason, JCP, 100 (1992), 25                                      c
c                                                                       c
c  The code can handle both rigid walls as well as periodic boundaries  c
c-----------------------------------------------------------------------c
c  A new code for the method developed in Ozen's thesis. Written in May c
c  and June 1991 at NASA-Lewis by gt---                                 c
c  Version 1.0 a two-dimensional version of the above....oct 1991       c
c  Version 1.1 modified in May 1992                                     c
c  Version 2.0 modified in June 1992                                    c
c  Version 2.1 modified in Sept 1992                                    c 
c  Version 2.2 modified in Nov. 1992  <=== Current copy                 c
c-----------------------------------------------------------------------c
c----------                   NOTICE                  ------------
c Although the code is fully functional for several problems, some of 
c its general features are only ``half-way there''. Therefore, it is 
c possible to set up a problem that may look right but where some 
c aspects are treated incorrectly.
c-----------------------------------------------------------------------c
      program main
cc      parameter( nxp2=8, nyp2=8) 
cc      parameter( nxp2=18, nyp2=18) 
cc      parameter( nxp2=34, nyp2=18)
cc      parameter( nxp2=18, nyp2=34)
cc      parameter( nxp2=66, nyp2=34) 
       parameter( nxp2=34, nyp2=34) 
c      parameter( nxp2=66, nyp2=66)
cc        parameter( nxp2=130, nyp2=130)
      parameter( maxel=15000, maxpt=15000)
c---------------grid code-----------------------------------------------c
      dimension u(nxp2,nyp2), v(nxp2,nyp2), wallv(4)
      dimension ut(nxp2,nyp2),vt(nxp2,nyp2)
      dimension uo(nxp2,nyp2),vo(nxp2,nyp2),roo(nxp2,nyp2)
      dimension r(nxp2,nyp2), ro(nxp2,nyp2),p(nxp2,nyp2)
      dimension tmp1(nxp2,nyp2), tmp2(nxp2,nyp2)
      dimension tmp3(nxp2,nyp2), tmp4(nxp2,nyp2)
      real m(nxp2,nyp2),m1,m2
c---------------grid/front code-----------------------------------------c
      dimension temp1(nxp2+2,nyp2+2), temp2(nxp2+2,nyp2+2)
      dimension rtmp(nxp2-1,nyp2-1)
c---------------front tracking code-------------------------------------c
      integer ptcon,elcon,bptcon,belcon,np,ffp,lfp,fep,ne,ffe,lfe,fee
      dimension icp(maxel,2), ine(maxel,2), pt(maxpt,2), t1(maxel,2)
      dimension ptcon(maxpt),elcon(maxel),iptmp(maxpt),  t2(maxel,2)
      dimension bptcon(maxpt),belcon(maxel), pto(maxpt,2)
      dimension cv(maxel,2),cnt(maxel,2), elprop(maxel,3)
      dimension fp(20,6),ip(20)
c-----------------------------------------------------------------------c 
c234567890123456789012345678901234567890123456789012345678901234567890123
 
      character*20 fname, xname, statfile, datafile, filmfile,command
      logical restart, movie

      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl
      common/bounds/ibdry,jbdry,nxlast,nylast
      common/front/np,ffp,lfp,fep, ne,ffe,lfe,fee,amax,amin
      common/press/itmax,xerr,beta
      common/flprop/r1,r2,m1,m2,gx,gy,rro 
      common/fparam/xmv,mv,fxl,fyl
      common/source/num_srce,fsource(20),isource_pt(20)
 
      data noutf,noutv,noutp,nbkg,nbkf,noutr,nouto,nouts,noutu
     1    /    1,    1,    1,   1,   1,    1,    1,     1,   1/
      nprh=1
cc  These parameters should be inputed 
 
      call greadin(maxstep,endtime,r0,nxp2,nyp2)
      write(*,*)' type in bdry and bdry velocities'
      read(*,*)ibdry,jbdry,wallv(1),wallv(2),wallv(3),wallv(4)
      read(*,*)npvel,nppr,nback

cc      read(*,*)nprh,nppr,npom,npsf,npuel,npvel
cc      { Determine what to print out}
         write(*,*)' TYPE name of output files'
         read (5,'(a10)')fname
         lfnm=ltr(fname,20)+1 
         
       command='mkdir '
       command(7:7+lfnm)=fname
       call system(command)
       xname=fname
       fname(lfnm:lfnm+1)='/'
ccc       fname(lfnm+1:lfnm+lfnm)=xname
       lfnm=ltr(fname,20)+1
                
      datafile=fname
      datafile(lfnm:lfnm+4)='data'
      open(9,file=datafile,status='new')
        write(9,*)nxp2,nyp2,xl,yl,ibdry,jbdry
        write(9,*)wallv(1),wallv(2),wallv(3),wallv(4)
      close(9,status='keep') 
                                
      statfile=fname
      statfile(lfnm:lfnm+4)='stat'
      open(9,file=statfile,status='new')  
      write(9,*)' time, area, xc, yc, div, ekin, xmom, ymom'
           fxl=xl
           fyl=yl
           xmv=30.0
           mv=30

      call gsetup(ut,vt,ro,p,tmp1,tmp2,tmp3,nxp2,nyp2)

       iord=2
cc                  { Order of time integration, 1 or 2 ****************}

c set parameters for front restructuring
         amax=0.8/hxi
         amin=0.2/hxi
         am=0.5*(amax+amin)
         xclose=am
cc      { Initial spacing of points}
         write(*,400)amax,amin
400      format(' front restr., amax, amin:',2e13.5)

         write(*,*)'IS this a restart, true or false?'
         read(*,*)restart
      if(restart)then
         call gbackin(time,u,v,nxp2,nyp2)
         call fbackin(pt,icp,ine,ptcon,elcon,bptcon,belcon,
     1                                      elprop,maxpt,maxel)
      else
          time=0.0
          call ginit(u,v,tmp1,tmp2,rtmp,r,m,nxp2,nyp2) 
          call finit(nfronts,fp,ip,pt,icp,ine,ptcon,
     1         elcon,bptcon,belcon,elprop,maxpt,maxel,am)
      end if

        call setdens(r0,nfronts,fp,ip,r,nxp2,nyp2)
c        call gprint(r,1,time,nxp2,nyp2,fname,noutr,lfnm,nxp2,nyp2) 
        call dens2(r,rtmp,temp1,temp2,nxp2,nyp2, 
     1     pt,ptcon,icp,elcon,elprop,maxpt,maxel)
c        call gprint(r,1,time,nxp2,nyp2,fname,noutr,lfnm,nxp2,nyp2) 

      call grro(rro,r,nxp2,nyp2)
      call gvisc(m,r,nxp2,nyp2)

       write(*,*)' TYPE in first printout time'
       read(*,1000)ptime
1000   format(e12.5)
                
c find smallest and largest viscosity and set viscous time step
           bm=m1/r1
           if((m2/r2) .gt. bm)bm=m2/r2
           sm=m2/r2
           if(sm .gt. (m1/r1))sm=m1/r1
           dtvisc=1.0/(4.0*bm*hxi*hxi)
 
      do 100 nstep=1,maxstep
cc                      { START Time Integration} 
c set the size of next time step:
           call gvelmax(u,v,vmax,nxp2,nyp2)
           dt=0.5*dtvisc
           if(vmax .gt. 1.0e-04)then
             dtadv=0.5*2.0*sm/vmax
             if(dtadv .lt. dt)dt=dtadv
           end if
           write(*,2000)nstep,time,vmax,dt
2000       format(' Timestep: ',i5,' Time: ',e12.5,
     1                   ' vmax: ',e12.5,' dt: ',e12.5)

c--------------printout if wanted---------------                        
c calculate front statistics if wanted (not necessary):                 
          call fstatis(pt,icp,ine,ptcon,elcon,bptcon,belcon,            
     1                       cv,elprop,maxpt,maxel,a1,a2,a3)            
          call gstat(time,r,u,v,nxp2,nyp2,a4,a5,a6,a7)                  
          write(9,*)time,a1,a2,a3,a4,a5,a6,a7

c printout if wanted: 
      if( (time+0.5*dt) .gt. ptime)then
          write(*,3000)time
3000      format(' Printout at Time: ',e12.5)
          call fprint(pt,ptcon,iptmp,icp,elcon,maxpt,maxel,
     1                                     fname,noutf,lfnm)
cc  { The front}

          if(nprh .eq. 1)
     1     call gprint(r,1,time,nxp2,nyp2,fname,noutr,lfnm,nxp2,nyp2)
cc                                { Density}
          if(nppr .eq. 1)
     1     call gprint(p,2,time,nxp2,nyp2,fname,noutp,lfnm,nxp2,nyp2)
cc                                { Pressure}

        if(npvel .eq. 1)then
         call gprint(u,5,time,nxp2,nyp2,fname,noutu,lfnm,nxlast+1,nyp2)
         call gprint(v,6,time,nxp2,nyp2,fname,noutv,lfnm,nxp2,nylast+1)
        end if

          write(*,*)' TYPE in next printout time'
          read(*,*)ptime
      end if
c--------------printout done--------------------

       if(mod(nstep,nback).eq.0)then
        call gbackout(time,u,v,nxp2,nyp2,fname,nbkg,lfnm)
        call fbackout(pt,icp,ine,ptcon,elcon,bptcon,belcon,
     1                       elprop,maxpt,maxel,fname,nbkf,lfnm)
       end if
      If(time .gt. endtime)go to 101

c restructure the front:
      call fregrid(pt,icp,ine,ptcon,elcon,bptcon,belcon,
     1                                elprop,maxpt,maxel)

       if(iord.eq.2)then
         do 120 i=1,nxp2
         do 120 j=1,nyp2
         uo(i,j)=u(i,j)
         vo(i,j)=v(i,j) 
         roo(i,j)=r(i,j)
120      continue
         lp=ffp
          do 130 kk=1,np
          pto(lp,1)=pt(lp,1)
          pto(lp,2)=pt(lp,2)
          lp=ptcon(lp)
130       continue
       end if

      do 110 iter=1,iord

c save old density 
      do 5 i=1,nxp2
      do 5 j=1,nyp2
       ro(i,j)=r(i,j)
5     continue
      call grro(rro,r,nxp2,nyp2) 

c advect the density by front tracking: 
      call xmovefr(dt,u,v,temp1,temp2,nxp2,nyp2,
     1     pt,ptcon,maxpt) 

c generating the density field:
      call dens2(r,rtmp,temp1,temp2,nxp2,nyp2, 
     1     pt,ptcon,icp,elcon,elprop,maxpt,maxel) 


c find unprojected velocities (incl. advection, diffusion and bdy forces):   
      call gadvect(dt,u,v,ut,vt,m,r,ro,nxp2,nyp2,wallv) 

c add surface tension: 
       call fcurv(cv,t1,t2,cnt,pt,icp,ine,ptcon,elcon,maxel,maxpt)

       call xsurft(tmp1,tmp2,temp1,temp2,nxp2,nyp2,cv,
     1     cnt,pt,icp,elcon,elprop,maxpt,maxel)
     
c add surface tension to unprojected velocities: 
       do 80 i=2,nxlast
       do 80 j=2,nyp1
80       ut(i,j)=ut(i,j)+dt*tmp1(i,j)/(0.5*(r(i+1,j)+r(i,j))) 

       do 81 i=2,nxp1
       do 81 j=2,nylast
81       vt(i,j)=vt(i,j)+dt*tmp2(i,j)/(0.5*(r(i,j+1)+r(i,j)))

          call gbdry(ut,nxp2,nyp2)
          call gbdry(vt,nxp2,nyp2) 

c find the pressure. First calculate the source term:

      do 90 i=2,nxp1
      do 90 j=2,nyp1
      tmp1(i,j)= 0.5*(hxi*(ut(i,j)-ut(i-1,j))+
     1                  hyi*(vt(i,j)-vt(i,j-1)))/dt
90    continue

      call gpressur(r,p,tmp1,tmp2,tmp3,tmp4,nxp2,nyp2) 
                   
c update the velocities:
      do 91 i=2,nxlast 
      do 91 j=2,nyp1
91    u(i,j)=ut(i,j)-2.0*hxi*dt*(p(i+1,j)-p(i,j))/(r(i+1,j)+r(i,j))

      do 92 i=2,nxp1
      do 92 j=2,nylast
92    v(i,j)=vt(i,j)-2.0*hyi*dt*(p(i,j+1)-p(i,j))/(r(i,j+1)+r(i,j))

c correct boundaries:
          call gbdry(u,nxp2,nyp2)
          call gbdry(v,nxp2,nyp2) 

c update the viscosity, using the density field
      call gvisc(m,r,nxp2,nyp2)

110   continue
       if(iord.eq.2)then
         do 140 i=1,nxp2
         do 140 j=1,nyp2
         u(i,j)=0.5*(u(i,j)+uo(i,j))
         v(i,j)=0.5*(v(i,j)+vo(i,j)) 
         r(i,j)=0.5*(r(i,j)+roo(i,j))
140      continue
         lp=ffp
          do 150 kk=1,np
          pt(lp,1)=0.5*(pt(lp,1)+pto(lp,1))
          pt(lp,2)=0.5*(pt(lp,2)+pto(lp,2))
          lp=ptcon(lp)
150       continue
                
          call gvisc(m,r,nxp2,nyp2)
       end if  

      time=time+dt

100   continue 
cc                             { END Time Integration}
101   continue

      close(9,status='keep')

        call gbackout(time,u,v,nxp2,nyp2,fname,nbkg,lfnm)
        call fbackout(pt,icp,ine,ptcon,elcon,bptcon,belcon,
     1                       elprop,maxpt,maxel,fname,nbkf,lfnm)
      stop
      end   
c-----------------------------------------------------------------------c
c-----------------------------------------------------------------------c 
c subroutine to set grid density from initial front.
c ip and fp         ip   fp: 1      2      3                  4
c cylinder:           type=1,radius,xc,     yc,                density jump
c wave left to right: type=3,amplit,total length in x,mean y, density jump
c wave right to left: type=2,amplit,total length in x,mean y, density jump
c 2 density is set below surface; 3 density is set above surface

      subroutine setdens(r0,nfronts,fp,ip,r,nxp2,nyp2) 
      dimension r(nxp2,nyp2), fp(20,6), ip(20)
      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl 
      common/bounds/ibdry,jbdry,nxlast,nylast
      common/fparam/xmv,mv,fxl,fyl

      If(nfronts.gt.20)then
        write(*,*)'ERROR: Increase dimension of frprop'
        stop
      end if  
      pi=4.*ATAN(1.)

      do 10 j=1,nyp2
      do 10 i=1,nxp2
10    r(i,j)=r0

      do 20 l=1,nfronts

      do 30 j=1,nyp2
      do 30 i=1,nxp2 
       x=float(i-1)/hxi
       y=float(j-1)/hyi
         if(ip(l).eq.1)then
             xx=abs(x-fp(l,2))
             if(ibdry.eq.0)then
               if(xx.gt.(0.5*xl))xx=xl-xx
             end if
             yy=abs(y-fp(l,3))
             if(jbdry.eq.0)then
               if(yy.gt.(0.5*yl))yy=yl-yy
             end if

            xfl=sqrt(xx**2+yy**2)-fp(l,1)
          end if

         if(ip(l).eq.2)
     1    xfl=y-(fp(l,3)+fp(l,1)*sin(2.0*pi*x/fp(l,2)))
         if(ip(l).eq.3)
     1      xfl=(fp(l,3)+fp(l,1)*sin(2.0*pi*x/fp(l,2)))-y

      if(xfl.le.0.0)r(i,j)=r(i,j)+fp(l,4)
30    continue

20    continue

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccc subroutines for the grid code cccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c-----------------------------------------------------------------------c 
      subroutine gstat(time,r,u,v,nxp2,nyp2,div,ekin,xmom,ymom)
      dimension u(nxp2,nyp2), v(nxp2,nyp2) 
      dimension r(nxp2,nyp2)
      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl 
      common/bounds/ibdry,jbdry,nxlast,nylast
      dimension work(10000),bd(1)
      
      div=0.0
      xmom=0.0
      ymom=0.0 
      ekin=0.0 
      adivx=0.0
      adivy=0.0 
      dmax=0.0 
      iim=0
      jjm=0
      do 20 i=2,nxp1
      do 20 j=2,nyp1
       dd=abs(hxi*(u(i,j)-u(i-1,j))+hyi*(v(i,j)-v(i,j-1)) )
       if(dd.gt.dmax)then
        dmax=dd
        iim=i
        jjm=j
       end if
       div=div+dd
       xmom=xmom+u(i,j)*0.5*(r(i,j)+r(i+1,j))     
       ymom=ymom+v(i,j)*0.5*(r(i,j)+r(i,j+1))
       ekin=ekin+0.5*r(i+1,j+1)*( (0.5*(u(i+1,j)+u(i,j)))**2+ 
     1                     (0.5*(v(i,j+1)+v(i,j)))**2) 
       adivx=adivx+abs( hxi*(u(i,j)-u(i-1,j)) )
       adivy=adivy+abs( hyi*(v(i,j)-v(i,j-1)) )
20    continue 
        write(*,*)' DIV2: ',div,xmom,ymom,ekin

      return
      end
c-----------------------------------------------------------------------c 
      subroutine gadvect(dt,u,v,ut,vt,m,r,ro,nxp2,nyp2,wallv)
c subroutine to calculate the advection and diffusion terms on a
c staggered grid. The routine returns ut,vt and wt for points 2 to nxp1,
c and so on.

      dimension u(nxp2,nyp2), v(nxp2,nyp2), wallv(4)
      dimension ut(nxp2,nyp2),vt(nxp2,nyp2)
      dimension r(nxp2,nyp2), ro(nxp2,nyp2)
      real m(nxp2,nyp2),m1,m2 
      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl
      common/bounds/ibdry,jbdry,nxlast,nylast
      common/flprop/r1,r2,m1,m2,gx,gy,rro 

      do 910 j=1,nyp2
      do 910 i=1,nxp2
       ut(i,j)=0.0
       vt(i,j)=0.0
910   continue

      if(ibdry.eq.1)then
cc       { set boundary velocities for solid wall}
        do 10 j=1,nyp1
cc         { first for the first derivative}
          v(1,j)=2.0*wallv(1)-v(2,j)
10        v(nxp2,j)=2.0*wallv(3)-v(nxp1,j)
      end if
      if(jbdry.eq.1)then
        do 20 i=1,nxp1
          u(i,1)=2.0*wallv(2)-u(i,2)
20        u(i,nyp2)=2.0*wallv(4)-u(i,nyp1)
      end if
             
      do 100 i=2,nxlast
cc     { Advection terms}
      do 100 j=2,nyp1
      ut(i,j) =
     1  -hxi*(ro(i+1,j)*( (0.5*(u(i+1,j)+u(i,j)))**2)-
     1                ro(i,j)*( (0.5*(u(i,j)+u(i-1,j)))**2) ) 
     2  -hyi*( 0.25*(ro(i,j)+ro(i+1,j)+ro(i+1,j+1)+ro(i,j+1))
     2     *0.5*(u(i,j)+u(i,j+1))*0.5*(v(i,j)+v(i+1,j))-
     2         0.25*(ro(i,j)+ro(i+1,j)+ro(i+1,j-1)+ro(i,j-1))
     2     *0.5*(u(i,j)+u(i,j-1))*0.5*(v(i,j-1)+v(i+1,j-1)))
100   continue

      do 101 i=2,nxp1
      do 101 j=2,nylast
      vt(i,j) =
     2  -hxi*( 0.25*(ro(i,j)+ro(i+1,j)+ro(i+1,j+1)+ro(i,j+1))
     2     *0.5*(u(i,j)+u(i,j+1))*0.5*(v(i,j)+v(i+1,j))-
     2         0.25*(ro(i,j)+ro(i,j+1)+ro(i-1,j+1)+ro(i-1,j))
     2     *0.5*(u(i-1,j)+u(i-1,j+1))*0.5*(v(i,j)+v(i-1,j)))
     1  -hyi*(ro(i,j+1)*( (0.5*(v(i,j+1)+v(i,j)))**2)-
     1                 ro(i,j)*( (0.5*(v(i,j)+v(i,j-1)))**2) ) 
101   continue


      if(ibdry.eq.1)then
cc      { set boundary velocities for solid wall}
        do 30 j=1,nyp1
cc         { then for the second derivative}
          v(1,j)=(v(3,j)-6.0*v(2,j)+8.0*wallv(1))/3.0
30        v(nxp2,j)=(v(nx,j)-6.0*v(nxp1,j)+8.0*wallv(3))/3.0
      end if
      if(jbdry.eq.1)then
        do 40 i=1,nxp1
          u(i,1)=(u(i,3)-6.0*u(i,2)+8.0*wallv(2))/3.0
40        u(i,nyp2)=(u(i,ny)-6.0*u(i,nyp1)+8.0*wallv(4))/3.0
      end if

      do 200 i=2,nxlast
cc       { viscous terms}
      do 200 j=2,nyp1
      ut(i,j) = ut(i,j)
     4  +hxi*2.*(m(i+1,j)*hxi*(u(i+1,j)-u(i,j)) -
     4           m(i,j)  *hxi*(u(i,j)-u(i-1,j)) )
     5  +hyi*( 0.25*(m(i,j)+m(i+1,j)+m(i+1,j+1)+m(i,j+1))*
     5    (hyi*(u(i,j+1)-u(i,j)) + hxi*(v(i+1,j)-v(i,j)) ) -
     5         0.25*(m(i,j)+m(i+1,j)+m(i+1,j-1)+m(i,j-1))*
     5    (hyi*(u(i,j)-u(i,j-1))+ hxi*(v(i+1,j-1)- v(i,j-1))) )
200   continue
      do 201 i=2,nxp1
      do 201 j=2,nylast
      vt(i,j) = vt(i,j)
     5  +hxi*( 0.25*(m(i,j)+m(i+1,j)+m(i+1,j+1)+m(i,j+1))*
     5    (hyi*(u(i,j+1)-u(i,j)) + hxi*(v(i+1,j)-v(i,j)) ) -
     5         0.25*(m(i,j)+m(i,j+1)+m(i-1,j+1)+m(i-1,j))*
     5    (hyi*(u(i-1,j+1)-u(i-1,j))+ hxi*(v(i,j)- v(i-1,j))) )
     4  +hyi*2.*(m(i,j+1)*hyi*(v(i,j+1)-v(i,j)) -
     4           m(i,j)  *hyi*(v(i,j)-v(i,j-1)) )
201   continue      

c find unprojected velocities
      do 300 i=2,nxlast
      do 300 j=2,nyp1
       ut(i,j)=( 0.5*(ro(i+1,j)+ro(i,j))*u(i,j)+dt*(ut(i,j) + 
     1      (rro - 0.5*(ro(i+1,j)+ro(i,j)) )*gx) )
     2           /(0.5*(r(i+1,j)+r(i,j) ))
300   continue
      do 301 i=2,nxp1
      do 301 j=2,nylast
       vt(i,j)=( 0.5*(ro(i,j+1)+ro(i,j))*v(i,j)+dt*(vt(i,j) +
     1      (rro - 0.5*(ro(i,j+1)+ro(i,j)) )*gy) )
     2           /(0.5*(r(i,j+1)+r(i,j) )) 
301   continue    

       call gbdry(ut,nxp2,nyp2)
cc   { Set boundaries for periodic case}
       call gbdry(vt,nxp2,nyp2)
cc   { Needed in pressure calculations}
cc   THESE ARE ALSO SET IN MAIN PROGRAM, CHECK IF}
cc                                  { THESE CAN BE DELETED} 
      return
      end
c-----------------------------------------------------------------------c
      subroutine gbdry(a,nxp2,nyp2)
c subroutine to set the 1 and the p2 elements of a
c periodic array
      dimension a(nxp2,nyp2) 
      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl
      common/bounds/ibdry,jbdry,nxlast,nylast

      if(ibdry.eq.0)then
        do 10 j=1,nyp2
          a(1,j)=a(nxp1,j)
          a(nxp2,j)=a(2,j)
10      continue
      end if

      if(jbdry.eq.0)then
        do 20 i=1,nxp2
          a(i,1)=a(i,nyp1)
          a(i,nyp2)=a(i,2)
20      continue
      end if
      return
      end    
c-----------------------------------------------------------------------c
      subroutine gbdry2(a,nxp2,nyp2)
c subroutine to set the p1 and the p2 elements of a
c periodic array (mostly used in relation to the front
c tracking

      dimension a(nxp2,nyp2)
      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl
      common/bounds/ibdry,jbdry,nxlast,nylast

        do 10 j=1,nyp1-1
          a(nxp1,j)=a(1,j)
          a(nxp2,j)=a(2,j)
10      continue

        do 20 i=1,nxp2
          a(i,nyp1)=a(i,1)
          a(i,nyp2)=a(i,2)
20      continue

      return
      end
c-----------------------------------------------------------------------c
      subroutine gpressur(r,p,tmp1,tmp2,tmp3,tmp4,nxp2,nyp2)
c subroutine to setup, and solve the pressure equation.

      dimension r(nxp2,nyp2), p(nxp2,nyp2)
      dimension tmp1(nxp2,nyp2), tmp2(nxp2,nyp2)
      dimension tmp3(nxp2,nyp2), tmp4(nxp2,nyp2)
      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl
      common/bounds/ibdry,jbdry,nxlast,nylast
      common/press/itmax,xerr,beta

cc       { set boundary velocities for solid wall}
      if(ibdry.eq.1)then
        do 5 j=1,nyp2
          p(1,j)=0.0
          p(nxp2,j)=0.0
          r(1,j)=1.0e+09
          r(nxp2,j)=1.0e+09
5       continue
      end if
      if(jbdry.eq.1)then 
        do 6 i=1,nxp1
          p(i,1)=0.0
          p(i,nyp2)=0.0
          r(i,1)=1.0e+09
          r(i,nyp2)=1.0e+09
6       continue
      end if
 
cc     { calculate the coefficient for p(i,j)}
      do 20 i=2,nxp1
      do 20 j=2,nyp1
      tmp2(i,j)=1.0/( hxi*hxi/(r(i+1,j)+r(i,j))+
     1                  hxi*hxi/(r(i-1,j)+r(i,j))+
     2 hyi*hyi/(r(i,j+1)+r(i,j))+hyi*hyi/(r(i,j-1)+r(i,j)) )
20    continue

c iterate for the pressure, Gauss-Seidler Iteration with SOR

      do 100 iter=1,itmax 

      do 30 i=2,nxp1
      do 30 j=2,nyp1
       tmp3(i,j)=p(i,j)
30    continue   

c-----Red and Black SOR--------------------------
c set bdry values for tmp4
      do 39 i=1,nxp2
      do 39 j=1,nyp2
       tmp4(i,j)=p(i,j)
39    continue

      do 40 i=2,nx,2
      do 40 j=2,ny,2
      tmp4(i,j)=(1.0-beta)*p(i,j)+beta* tmp2(i,j)*(
     1  hxi*hxi*(p(i+1,j)/(r(i+1,j)+r(i,j))+
     2           p(i-1,j)/(r(i-1,j)+r(i,j)))+
     3  hyi*hyi*(p(i,j+1)/(r(i,j+1)+r(i,j))+
     4           p(i,j-1)/(r(i,j-1)+r(i,j))) - tmp1(i,j))
40    continue
      do 42 i=3,nxp1,2
      do 42 j=3,nyp1,2
      tmp4(i,j)=(1.0-beta)*p(i,j)+beta* tmp2(i,j)*(
     1  hxi*hxi*(p(i+1,j)/(r(i+1,j)+r(i,j))+
     2           p(i-1,j)/(r(i-1,j)+r(i,j)))+
     3  hyi*hyi*(p(i,j+1)/(r(i,j+1)+r(i,j))+
     4           p(i,j-1)/(r(i,j-1)+r(i,j))) - tmp1(i,j))
42    continue
 
      do 45 i=2,nx,2
      do 45 j=3,nyp1,2
      p(i,j)=(1.0-beta)*p(i,j)+beta* tmp2(i,j)*(
     1  hxi*hxi*(tmp4(i+1,j)/(r(i+1,j)+r(i,j))+
     2           tmp4(i-1,j)/(r(i-1,j)+r(i,j)))+
     3  hyi*hyi*(tmp4(i,j+1)/(r(i,j+1)+r(i,j))+
     4           tmp4(i,j-1)/(r(i,j-1)+r(i,j))) - tmp1(i,j))
45    continue
      do 47 i=3,nxp1,2
      do 47 j=2,ny,2
      p(i,j)=(1.0-beta)*p(i,j)+beta* tmp2(i,j)*(
     1  hxi*hxi*(tmp4(i+1,j)/(r(i+1,j)+r(i,j))+
     2           tmp4(i-1,j)/(r(i-1,j)+r(i,j)))+
     3  hyi*hyi*(tmp4(i,j+1)/(r(i,j+1)+r(i,j))+
     4           tmp4(i,j-1)/(r(i,j-1)+r(i,j))) - tmp1(i,j))
47    continue

      do 48 i=2,nx,2
      do 48 j=2,ny,2
      p(i,j)=tmp4(i,j)
48    continue
      do 49 i=3,nxp1,2
      do 49 j=3,nyp1,2
      p(i,j)=tmp4(i,j)
49    continue
c--------------------------------------------------------------------
      call gbdry(p,nxp2,nyp2)

c check for convergence
      err=0.0
      do 50 i=2,nxp1
      do 50 j=2,nyp1
        err=err+abs(p(i,j)-tmp3(i,j))
50    continue 
        err=err/float(nx*ny)
      if(err .lt. xerr)go to 101

100   continue
101   continue
        write(*,1000)err,iter
1000    format(' Pressure Solution. Error = ',e12.5,
     1                    ' after ',i5,' iterations')

cc       { restore boundary density} 
cc        { This is a somewhat clumsy way}
cc        { and should be reexamined}
      if(ibdry.eq.1)then
        do 7 j=1,nyp2
          r(1,j)=r(2,j)
7         r(nxp2,j)=r(nxp1,j)
      end if
      if(jbdry.eq.1)then 
        do 8 i=1,nxp1
          r(i,1)=r(i,2)
8         r(i,nyp2)=r(i,nyp1)
      end if 

      return
      end 
c-----------------------------------------------------------------------c         
      subroutine grro(rro,r,nxp2,nyp2)
      dimension r(nxp2,nyp2) 
      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl
      rro=0.0
      do 10 i=2,nxp1
      do 10 j=2,nyp1
       rro=rro+r(i,j)
10    continue
       rro=rro/float((nxp1-1)*(nyp1-1))
       write(*,100)rro
100    format(' average density: ',e13.5)
      return
      end
c-----------------------------------------------------------------------c
      subroutine greadin(maxstep,endtime,r0,nxp2,nyp2)
      real m1,m2 
      character*100 junk 
      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl
      common/press/itmax,xerr,beta  
      common/flprop/r1,r2,m1,m2,gx,gy,rro

      write(*,*)' First read in five lines of text'
        do 20 i=1,5
        read(*,1)junk
20      write(*,1)junk
1       format(a100)
      write(*,*)' type in the initial data'
      write(*,*)' type in: maxstep,endtime,xl,yl'
      read(*,*)maxstep,endtime,xl,yl
      write(*,*)' type in: r0,r1,r2,m1,m2,gx,gy'
      read(*,*)r0,r1,r2,m1,m2,gx,gy
      write(*,*)' type in parameters for pressure sol.:',
     1                 ' itmax, xerr, beta'
      read(*,*)itmax,xerr,beta

      write(*,300)maxstep,xl,yl,nxp2,nyp2,endtime
      write(*,301)r1,r2,m1,m2,gx,gy
      write(*,302)itmax,xerr,beta

300   format(' Initial Data: ',/,
     1  ' Maximum time steps:......... ',i5,/,
     2  ' Size of Box, xl, yl:......... ',2e12.5,/,
     3  ' Size of grid, nxp2, nyp2:.... ',2i5,/, 
     4  ' Final Time:.................. ',e12.5)
301   format(' Material Properties: ',/,
     1  ' Densities, r1, r2:........... ',2e12.5,/,
     2  ' Viscosities, m1, m2:......... ',2e12.5,/,
     4  ' Gravity, gx,gy:.............. ',2e13.5)
302   format(' Pressure Solution: ',/,
     1  ' Max. Iterations:............. ',i5,/,
     2  ' Max. Error:.................. ',e12.5,/, 
     3  ' SOR Parameters, beta:........ ',e12.5)

      return
      end
c-----------------------------------------------------------------------c
      subroutine gsetup(ut,vt,ro,p,tmp1,tmp2,tmp3,nxp2,nyp2)
      dimension ut(nxp2,nyp2),vt(nxp2,nyp2)
      dimension tmp1(nxp2,nyp2), tmp2(nxp2,nyp2)
      dimension tmp3(nxp2,nyp2),ro(nxp2,nyp2)
      dimension p(nxp2,nyp2)
      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl
      common/bounds/ibdry,jbdry,nxlast,nylast
   
      nxp1=nxp2-1
      nyp1=nyp2-1
      nx=nxp2-2
      ny=nyp2-2
      nxlast=nxp1   
      nylast=nyp1
      if(ibdry.eq.1)nxlast=nx
      if(jbdry.eq.1)nylast=ny

      hxi=float(nxp1-1)/xl
      hyi=float(nyp1-1)/yl
                    
      do 10 i=1,nxp2
      do 10 j=1,nyp2
       ut(i,j)=0.0
       vt(i,j)=0.0
       ro(i,j)=0.0 
       p(i,j)=0.0
       tmp1(i,j)=0.0
       tmp2(i,j)=0.0
       tmp3(i,j)=0.0
10    continue  

      return
      end
c-----------------------------------------------------------------------c        
      subroutine ginit(u,v,tmp1,tmp2,rtmp,r,m,nxp2,nyp2)
      dimension u(nxp2,nyp2), v(nxp2,nyp2),tmp1(nxp2,nyp2)
      dimension tmp2(nxp2,nyp2), rtmp(nxp2-1,nyp2-1)
      dimension r(nxp2,nyp2)
      dimension work(10000),bd(1)
      real m(nxp2,nyp2), m1,m2
      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl 
      common/bounds/ibdry,jbdry,nxlast,nylast
      common/flprop/r1,r2,m1,m2,gx,gy,rro
 
c           Initial velocity field should be set here
      do 10 i=1,nxp2
      do 10 j=1,nyp2 
       y=float(j-1)/hyi
       u(i,j)=0.0
       v(i,j)=0.0
c        u(i,j)=2.0*y-1.0
10    continue   

100    continue
       call gbdry(u,nxp2,nyp2)
       call gbdry(v,nxp2,nyp2)

      return
      end
c-----------------------------------------------------------------------c
      subroutine gvisc(m,r,nxp2,nyp2)
      dimension r(nxp2,nyp2)
      real m(nxp2,nyp2),m1,m2 
      common/flprop/r1,r2,m1,m2,gx,gy,rro
 
      if(m1.ne.m2)then
       do 10 i=1,nxp2
       do 10 j=1,nyp2
10      m(i,j)= (m1*(r(i,j)-r2)-m2*(r(i,j)-r1) )/(r1-r2)
      else
       do 20 i=1,nxp2
       do 20 j=1,nyp2
20      m(i,j)=m1
      end if

      return
      end
c-----------------------------------------------------------------------c 
      subroutine gprint(a,itype,time,nxp2,nyp2,fname,n,lm,nxp,nyp)
      dimension a(nxp2,nyp2)
      character*20 outfile, fname, varname(6)
      character*1 num(10)  
      character*3 des(6)
      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl 
      data (num(i),i=1,10)/'0','1','2','3','4','5','6','7','8','9'/ 
      data (des(i),i=1,6)/'rho','prs','vrt','srf','uvl','vvl'/   
      data (varname(i),i=1,6)/'density','pressure','vorticity',
     1     'streamfunct','u-velocity','v-velocity'/ 

      outfile=fname 
      m1=mod(n,10)+1
      m0=mod((n/10),10)+1
      outfile(lm:lm+3)=des(itype)
      outfile(lm+3:lm+4)='0'      
      outfile(lm+4:lm+5)=num(m0) 
      outfile(lm+5:lm+6)=num(m1)
      open(8,file=outfile,status='new') 
      nyh=nyp2/2
      write(8,300)varname(itype),time
      write(8,310)nxp,nyp
      write(8,320)((float(i-1)/hxi),i=1,nxp)
      write(8,320)((float(j-1)/hyi),j=1,nyp)
      do 10 j=1,nyp
      write(8,320)(a(i,j),i=1,nxp)
10    continue
300   format(' ',a20, e12.5)
310   format(2i5)
320   format(8e13.5)
      close(8,status='keep') 
      n=n+1
      return
      end          

c-----------------------------------------------------------------------c 
 
      function ltr(string,l)
      character*1 string(l)
      character*1 blnk
      data blnk/' '/
      do 100 i=l,1,-1
       l1=i
       if(string(i).ne.blnk)go to 200
100   continue
      l1=0
200   continue
      ltr=l1
      return
      end

c-----------------------------------------------------------------------c 
      subroutine gvelmax(u,v,vmax,nxp2,nyp2)
c subroutine to find maximum velocity.
      dimension u(nxp2,nyp2), v(nxp2,nyp2)
      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl

      vmax=0.0
      do 10 i=2,nxp1
      do 10 j=2,nyp1
       vv=(0.5*(u(i,j)+u(i-1,j)))**2+(0.5*(v(i,j)+v(i,j-1)))**2
ccc u(i,j)**2+v(i,j)**2
       if(vv .gt. vmax)then
        vmax=vv  
        im=i
        jm=j
       end if
10    continue 
       velm=sqrt(vmax)
cc      write(*,*)' maxvel,i,j',velm,im,jm

      return
      end
c-----------------------------------------------------------------------c      
c below are routines to back up the results, and start from a backup file.
      subroutine gbackin(time,u,v,nxp2,nyp2) 
      dimension u(nxp2,nyp2), v(nxp2,nyp2)
      character*10 infile
      write(*,*)' type in name of the Grid Backup file'
      read(*,100)infile
100   format(1a10)
      open(8,file=infile,status='old')
      read(8,*)time
      do 10 i=1,nxp2
      do 10 j=1,nyp2
       read(8,*)u(i,j),v(i,j)
10    continue 
      close(8,status='keep') 
      write(*,110)time
110   format(' Read from backup file. Time: ',e13.5)
      return
      end
c-----------------------------------------------------------------------c
      subroutine gbackout(time,u,v,nxp2,nyp2,fname,n,lm)
      dimension u(nxp2,nyp2), v(nxp2,nyp2)
      character*20 outfile, fname
      character*1 num(10), nsuf(10) 
      data (num(i),i=1,10)/'0','1','2','3','4','5','6','7','8','9'/ 
      outfile=fname 
      m1=mod(n,10)+1
      m0=mod((n/10),10)+1
      outfile(lm:lm+3)='bkg'
      outfile(lm+3:lm+4)=num(m0) 
      outfile(lm+4:lm+5)=num(m1)
      open(8,file=outfile,status='new') 
      write(8,*)time 
      do 10 i=1,nxp2
      do 10 j=1,nyp2
       write(8,*)u(i,j),v(i,j)
10    continue 
      close(8,status='keep') 
      n=n+1
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccc subroutines that connect the grid and the front code cccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----------------------------------------------------------------------c

      subroutine xmovefr(dt,u,v,tmp1,tmp2,nxp2,nyp2,
     1     pt,ptcon,maxpt)

      dimension u(nxp2,nyp2), v(nxp2,nyp2)
      dimension tmp1(nxp2+2,nyp2+2), tmp2(nxp2+2,nyp2+2)
      dimension pt(maxpt,2), ptcon(maxpt)
      integer ptcon,np,ffp,lfp,fep, ne,ffe,lfe,fee

      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl
      common/bounds/ibdry,jbdry,nxlast,nylast
      common/front/np,ffp,lfp,fep, ne,ffe,lfe,fee,amax,amin 
      common/fparam/xmv,mv,fxl,fyl
      
      pi=4.*ATAN(1.)
      pd2=0.5*pi
      cnst=0.25**2 

c generate a velocity field at the pressure notes
      do 5 i=2,nxp1
      do 5 j=2,nyp1
       tmp1(i+1,j+1)=0.5*(u(i,j)+u(i-1,j))
       tmp2(i+1,j+1)=0.5*(v(i,j)+v(i,j-1))
5     continue 

      call gbd3(tmp1,tmp2,ibdry,jbdry,nxp2,nyp2)

c find new interface position
      k=ffp
      do 20 ll=1,np

c find grid box
       xt=amod((pt(k,1)+0.5/hxi)+xmv*xl,xl)
       yt=amod((pt(k,2)+0.5/hyi)+xmv*yl,yl)

       ir=1+int(nx*xt/xl)
       jr=1+int(ny*yt/yl)  

c interpolate the velocities 
       up=0.0
       vp=0.0

       do 30 i1=1,4
       do 30 j1=1,4
        ii=ir-2+i1   
        jj=jr-2+j1

         drx = 1. + cos((xt*hxi - float(ii-1))*pd2)
         dry = 1. + cos((yt*hyi - float(jj-1))*pd2)
   
        up = up+tmp1(ii+1,jj+1)*drx*dry*cnst
        vp = vp+tmp2(ii+1,jj+1)*drx*dry*cnst

30      continue  
                  
          pt(k,1)=pt(k,1)+dt*up   
          pt(k,2)=pt(k,2)+dt*vp  

        k=ptcon(k)
20      continue 

      return
      end
c-----------------------------------------------------------------------c
      subroutine gbd3(tmp1,tmp2,ibd,jbd,nxp2,nyp2)

      dimension tmp1(nxp2+2,nyp2+2),tmp2(nxp2+2,nyp2+2)
      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl
c for xmovefr
c 0--periodic; 1--solid wall

      if(jbd .eq. 0)then 
       do 60 i=1,nxp2+2
cc              { Periodic in y}
         tmp1(i,1)=tmp1(i,nyp1)
         tmp1(i,2)=tmp1(i,nyp2)
         tmp1(i,nyp2+1)=tmp1(i,3)
         tmp1(i,nyp2+2)=tmp1(i,4)

         tmp2(i,1)=tmp2(i,nyp1)
         tmp2(i,2)=tmp2(i,nyp2)
         tmp2(i,nyp2+1)=tmp2(i,3)
         tmp2(i,nyp2+2)=tmp2(i,4)
60     continue

      else 

       do 62 i=1,nxp2+2
cc              { solid wall in y}
         tmp1(i,2)=-tmp1(i,3)
         tmp1(i,1)=-tmp1(i,4)
         tmp1(i,nyp2+2)=-tmp1(i,nyp1)
         tmp1(i,nyp2+1)=-tmp1(i,nyp2)
         tmp2(i,2)=-tmp2(i,3)
         tmp2(i,1)=-tmp2(i,4)
         tmp2(i,nyp2+2)=-tmp2(i,nyp1)
         tmp2(i,nyp2+1)=-tmp2(i,nyp2)
62     continue
      end if  

      if(ibd .eq. 0)then 

       do 70 j=1,nyp2+2
cc              { Periodic in x}
         tmp1(1,j)=tmp1(nxp1,j)
         tmp1(2,j)=tmp1(nxp2,j)
         tmp1(nxp2+1,j)=tmp1(3,j)
         tmp1(nxp2+2,j)=tmp1(4,j)

         tmp2(1,j)=tmp2(nxp1,j)
         tmp2(2,j)=tmp2(nxp2,j)
         tmp2(nxp2+1,j)=tmp2(3,j)
         tmp2(nxp2+2,j)=tmp2(4,j)
70     continue
      else 

       do 72 j=1,nyp2+2
cc              { solid wall in x}
         tmp1(2,j)=-tmp1(3,j)
         tmp1(1,j)=-tmp1(4,j)
         tmp1(nxp2+2,j)=-tmp1(nxp1,j)
         tmp1(nxp2+1,j)=-tmp1(nxp2,j)

         tmp2(2,j)=-tmp2(3,j)
         tmp2(1,j)=-tmp2(4,j)
         tmp2(nxp2+2,j)=-tmp2(nxp1,j)
         tmp2(nxp2+1,j)=-tmp2(nxp2,j)
 72     continue
      end if        

      return
      end
c-----------------------------------------------------------------------c

      subroutine dens2(r,rtmp,tmp1,tmp2,nxp2,nyp2, 
     1     pt,ptcon,icp,elcon,elprop,maxpt,maxel)  

      dimension tmp1(nxp2+2,nyp2+2), tmp2(nxp2+2,nyp2+2)
      dimension rtmp(nxp2-1,nyp2-1)
      dimension r(nxp2,nyp2)
      dimension iii(2000),jjj(2000)

      dimension pt(maxpt,2), ptcon(maxpt)
      dimension icp(maxel,2),elcon(maxel),elprop(maxel,3)

      integer ptcon,elcon,np,ffp,lfp,fep, ne,ffe,lfe,fee

      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl
      common/bounds/ibdry,jbdry,nxlast,nylast
      common/front/np,ffp,lfp,fep, ne,ffe,lfe,fee,amax,amin 
      common/fparam/xmv,mv,fxl,fyl
      common/flprop/r1,r2,m1,m2,gx,gy,rro  
       xx=(xmv+0.5)*xl    
       yy=(xmv+0.5)*yl

c generate gradient vector
 
      sc=1./2.
      pi=4.*ATAN(1.)
      pd2=0.5*pi
      cnstx=(.25*hxi)**2 
      cnsty=(.25*hyi)**2

      do 10 i=1,nxp2+2
      do 10 j=1,nyp2+2  
       tmp1(i,j)=0.0
       tmp2(i,j)=0.0
10    continue 

      k=ffe
      do 30 ll=1,ne 
       p0x=pt(icp(k,1),1)
       p0y=pt(icp(k,1),2)
       p1x= p0x+amod((pt(icp(k,2),1)-p0x+xx),xl)-0.5*xl
       p1y= p0y+amod((pt(icp(k,2),2)-p0y+yy),yl)-0.5*yl

       xc=0.5*(p1x+p0x)
       yc=0.5*(p1y+p0y)
       x1=p1x-p0x
       y1=p1y-p0y

c calculate components of n.dA
      rx=-y1 *elprop(k,1)
      ry= x1 *elprop(k,1)
                      
c generate a grid value   
       xt=amod((xc+0.5/hxi)+xmv*xl,xl)
       yt=amod((yc+0.5/hyi)+xmv*yl,yl)

       ir=1+int(nx*xt/xl)
       jr=1+int(ny*yt/yl)
                    
       do 40 i1=1,4
       do 40 j1=1,4
        ii=ir-2+i1   
        jj=jr-2+j1

         drx = 1. + cos((xt*hxi - float(ii-1))*pd2)
         dry = 1. + cos((yt*hyi - float(jj-1))*pd2)

        tmp1(ii+1,jj+1)=tmp1(ii+1,jj+1)+rx*drx*dry*cnstx
        tmp2(ii+1,jj+1)=tmp2(ii+1,jj+1)+ry*drx*dry*cnsty
40      continue

      k=elcon(k)
30    continue   

c correct boundaries:
      call gbd1(tmp1,tmp2,ibdry,jbdry,nxp2,nyp2)

       do 900 i=1,nxp2
       do 900 j=1,nyp2
900        rtmp(i,j)=0.0

       do 80 i=2,nxp1
       do 80 j=2,nyp1
        rtmp(i,j)= 0.5*(hxi*(tmp1(i+1+1,j+1)-tmp1(i-1+1,j+1))+
     1        hyi*(tmp2(i+1,j+1+1)-tmp2(i+1,j-1+1)))
80     continue  

c set connection
      node=0
       do 100 i=2,nxp1
       do 100 j=2,nyp1
         cc=0.0
         do 105 io=1,3
         do 105 jo=1,3
105      cc=cc+rtmp(i+io-2,j+jo-2)
           if(abs(cc) .gt. 1.0e-03)then
            node=node+1
            iii(node)=i
            jjj(node)=j
           end if
100      continue 
       numnd=node

c then iterate

c-----------the following should be read in
      itmax=1000
      bet=1.2
      errmax=1.0e-05 
c------------------------------------------
      cf=0.5/(hxi*hxi+hyi*hyi)
      do 200 iter=1,itmax

       err=0.0 
         do 210 node=1,numnd
           i=iii(node)
           j=jjj(node)
           rold=r(i,j)
           r(i,j)=(1.0-bet)*r(i,j)+bet*cf*
     1        (hxi*hxi*(r(i+1,j)+r(i-1,j))+
     2       hyi*hyi*(r(i,j+1)+r(i,j-1))-rtmp(i,j))
           err=err+abs(rold-r(i,j))
210       continue

c set boundary terms...
      if(ibdry .eq. 0)then
        do 110 j=1,nyp2
          r(nxp2,j)=r(2,j)
          r(1,j)=r(nxp1,j)
110      continue
      else
        do 115 j=1,nyp2
          r(nxp2,j)=r(nxp1,j)
          r(1,j)=r(2,j)
115      continue
      end if

      if(jbdry .eq. 0)then
        do 120 i=1,nxp2
          r(i,nyp2)=r(i,2)
          r(i,1)=r(i,nyp1)
120      continue 
      else
        do 125 i=1,nxp2
          r(i,nyp2)=r(i,nyp1)
          r(i,1)=r(i,2)
125      continue 
      end if

c check convergence......
      if(err.lt.(errmax*float(numnd)))go to 300
200   continue
300   continue
         write(*,*)' DENSITY: ',err,' after ',iter,' iterations',
     1           ' using ',numnd,' nodes'

      return
      end
c-----------------------------------------------------------------------c
      subroutine gbd1(tmp1,tmp2,ibd,jbd,nxp2,nyp2)
c   for dens2. 0---periodic; 1---solid wall

      dimension tmp1(nxp2+2,nyp2+2),tmp2(nxp2+2,nyp2+2)
      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl

      if(jbd .eq. 0)then 
       do 60 i=1,nxp2+2
cc              { Periodic in y}
         tmp1(i,3)=tmp1(i,3)+tmp1(i,nyp2+1)
         tmp1(i,4)=tmp1(i,4)+tmp1(i,nyp2+2)
         tmp1(i,nyp1)=tmp1(i,nyp1)+tmp1(i,1)
         tmp1(i,nyp2)=tmp1(i,nyp2)+tmp1(i,2)
         tmp2(i,3)=tmp2(i,3)+tmp2(i,nyp2+1)
         tmp2(i,4)=tmp2(i,4)+tmp2(i,nyp2+2)
         tmp2(i,nyp1)=tmp2(i,nyp1)+tmp2(i,1)
         tmp2(i,nyp2)=tmp2(i,nyp2)+tmp2(i,2)
60     continue

       do 61 i=1,nxp2+2
         tmp1(i,1)=tmp1(i,nyp1)
         tmp1(i,2)=tmp1(i,nyp2)
         tmp1(i,nyp2+1)=tmp1(i,3)
         tmp1(i,nyp2+2)=tmp1(i,4)

         tmp2(i,1)=tmp2(i,nyp1)
         tmp2(i,2)=tmp2(i,nyp2)
         tmp2(i,nyp2+1)=tmp2(i,3)
         tmp2(i,nyp2+2)=tmp2(i,4)
61     continue
      else 

       do 62 i=1,nxp2+2
cc              { solid wall in y}
         tmp1(i,3)=tmp1(i,3)+tmp1(i,2)
         tmp1(i,4)=tmp1(i,4)+tmp1(i,1)
         tmp1(i,nyp1)=tmp1(i,nyp1)+tmp1(i,nyp2+2)
         tmp1(i,nyp2)=tmp1(i,nyp2)+tmp1(i,nyp2+1)

         tmp2(i,3)=tmp2(i,3)-tmp2(i,2)
         tmp2(i,4)=tmp2(i,4)-tmp2(i,1)
         tmp2(i,nyp1)=tmp2(i,nyp1)-tmp2(i,nyp2+2)
         tmp2(i,nyp2)=tmp2(i,nyp2)-tmp2(i,nyp2+1)
62     continue

       do 63 i=1,nxp2+2
         tmp1(i,2)=tmp1(i,3)
         tmp1(i,1)=tmp1(i,4)
         tmp1(i,nyp2+2)=tmp1(i,nyp1)
         tmp1(i,nyp2+1)=tmp1(i,nyp2)
         tmp2(i,2)=-tmp2(i,3)
         tmp2(i,1)=-tmp2(i,4)
         tmp2(i,nyp2+2)=-tmp2(i,nyp1)
         tmp2(i,nyp2+1)=-tmp2(i,nyp2)
63     continue
      end if  

      if(ibd .eq. 0)then 
       do 70 j=1,nyp2+2
cc              { Periodic in x}
         tmp1(3,j)=tmp1(3,j)+tmp1(nxp2+1,j)
         tmp1(4,j)=tmp1(4,j)+tmp1(nxp2+2,j)
         tmp1(nxp1,j)=tmp1(nxp1,j)+tmp1(1,j)
         tmp1(nxp2,j)=tmp1(nxp2,j)+tmp1(2,j)

         tmp2(3,j)=tmp2(3,j)+tmp2(nxp2+1,j)
         tmp2(4,j)=tmp2(4,j)+tmp2(nxp2+2,j)
         tmp2(nxp1,j)=tmp2(nxp1,j)+tmp2(1,j)
         tmp2(nxp2,j)=tmp2(nxp2,j)+tmp2(2,j)
70     continue

       do 71 j=1,nyp2+2
         tmp1(1,j)=tmp1(nxp1,j)
         tmp1(2,j)=tmp1(nxp2,j)
         tmp1(nxp2+1,j)=tmp1(3,j)
         tmp1(nxp2+2,j)=tmp1(4,j)

         tmp2(1,j)=tmp2(nxp1,j)
         tmp2(2,j)=tmp2(nxp2,j)
         tmp2(nxp2+1,j)=tmp2(3,j)
         tmp2(nxp2+2,j)=tmp2(4,j)
71     continue
      else 

       do 72 j=1,nyp2+2
cc              { solid wall in x}
         tmp1(3,j)=tmp1(3,j)-tmp1(2,j)
         tmp1(4,j)=tmp1(4,j)-tmp1(1,j)
         tmp1(nxp1,j)=tmp1(nxp1,j)-tmp1(nxp2+2,j)
         tmp1(nxp2,j)=tmp1(nxp2,j)-tmp1(nxp2+1,j)

         tmp2(3,j)=tmp2(3,j)+tmp2(2,j)
         tmp2(4,j)=tmp2(4,j)+tmp2(1,j)
         tmp2(nxp1,j)=tmp2(nxp1,j)+tmp2(nxp2+2,j)
         tmp2(nxp2,j)=tmp2(nxp2,j)+tmp2(nxp2+1,j)
72     continue

       do 73 j=1,nyp2+2
         tmp1(2,j)=-tmp1(3,j)
         tmp1(1,j)=-tmp1(4,j)
         tmp1(nxp2+2,j)=-tmp1(nxp1,j)
         tmp1(nxp2+1,j)=-tmp1(nxp2,j)

         tmp2(2,j)=tmp2(3,j)
         tmp2(1,j)=tmp2(4,j)
         tmp2(nxp2+2,j)=tmp2(nxp1,j)
         tmp2(nxp2+1,j)=tmp2(nxp2,j)
 73     continue
      end if        

      return
      end
c-----------------------------------------------------------------------c

      subroutine xsurft(fx,fy,tmp1,tmp2,nxp2,nyp2,cv,
     1     cnt,pt,icp,elcon,elprop,maxpt,maxel) 
      dimension fx(nxp2,nyp2),fy(nxp2,nyp2)
      dimension tmp1(nxp2+2,nyp2+2), tmp2(nxp2+2,nyp2+2)
      integer elcon,np,ffp,lfp,fep, ne,ffe,lfe,fee
      dimension icp(maxel,2), pt(maxpt,2), elcon(maxel)
      dimension cv(maxel,2), cnt(maxel,2), elprop(maxel,3)

      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl
      common/bounds/ibdry,jbdry,nxlast,nylast
      common/front/np,ffp,lfp,fep, ne,ffe,lfe,fee,amax,amin 
      common/fparam/xmv,mv,fxl,fyl
      common/flprop/r1,r2,m1,m2,gx,gy,rro  

      do 10 j=1,nyp2+2
      do 10 i=1,nxp2+2
      tmp1(i,j)=0.0
      tmp2(i,j)=0.0
10    continue

      pi=4.*ATAN(1.)
      pd2=0.5*pi
      cnstx=(.25*hxi)**2
      cnsty=(.25*hyi)**2
      dxh=0.5/hxi
      dyh=0.5/hyi   

      k=ffe
      do 30 ll=1,ne

      xc=cnt(k,1)
      yc=cnt(k,2)
       
c generate a grid value:
       xt=amod((xc+dxh)+xmv*xl,xl)
       yt=amod((yc+dyh)+xmv*yl,yl)
         xth=amod((xc)+xmv*xl,xl)
         yth=amod((yc)+xmv*yl,yl)

       ir=1+int(nx*xt/xl)
       jr=1+int(ny*yt/yl)
         irh=1+int(nx*xth/xl)
         jrh=1+int(ny*yth/yl)
                    
       do 40 i1=1,4
       do 40 j1=1,4
        ii=ir-2+i1   
        jj=jr-2+j1
          iih=irh-2+i1   
          jjh=jrh-2+j1

         drx = 1. + cos((xt*hxi - float(ii-1))*pd2)
         dry = 1. + cos((yt*hyi - float(jj-1))*pd2)
           drxh = 1. + cos((xth*hxi - float(iih-1))*pd2)
           dryh = 1. + cos((yth*hyi - float(jjh-1))*pd2)

        tmp1(iih+1,jj+1)=tmp1(iih+1,jj+1)+
     1                 drxh*dry*cnstx*cv(k,1)*elprop(k,2)
        tmp2(ii+1,jjh+1)=tmp2(ii+1,jjh+1)+
     1                 drx*dryh*cnsty*cv(k,2)*elprop(k,2)

40      continue 

      k=elcon(k)
30    continue   
               
c correct the ends of the box

      call gbd2(tmp1,tmp2,ibdry,jbdry,nxp2,nyp2)

c copy into fx,fy
      do 100 i=1,nxp2
      do 100 j=1,nyp2
      fx(i,j)=tmp1(i+1,j+1)
      fy(i,j)=tmp2(i+1,j+1)
100   continue
      return
      end 
c-----------------------------------------------------------------------c

      subroutine gbd2(tmp1,tmp2,ibd,jbd,nxp2,nyp2)
c   for xsurf. 0---periodic; 1---solid wall

      dimension tmp1(nxp2+2,nyp2+2),tmp2(nxp2+2,nyp2+2)
      common/grid/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl

      if(jbd .eq. 0)then 
       do 60 i=1,nxp2+2
cc              { Periodic in y}
         tmp1(i,3)=tmp1(i,3)+tmp1(i,nyp2+1)
         tmp1(i,4)=tmp1(i,4)+tmp1(i,nyp2+2)
         tmp1(i,nyp1)=tmp1(i,nyp1)+tmp1(i,1)
         tmp1(i,nyp2)=tmp1(i,nyp2)+tmp1(i,2)

         tmp2(i,1)=tmp2(i,1)+tmp2(i,nyp1)
         tmp2(i,2)=tmp2(i,2)+tmp2(i,nyp2)
         tmp2(i,3)=tmp2(i,3)+tmp2(i,nyp2+1)

60     continue

       do 61 i=1,nxp2+2
         tmp1(i,1)=tmp1(i,nyp1)
         tmp1(i,2)=tmp1(i,nyp2)
         tmp1(i,nyp2+1)=tmp1(i,3)
         tmp1(i,nyp2+2)=tmp1(i,4)

         tmp2(i,nyp1)=tmp2(i,1)
         tmp2(i,nyp2)=tmp2(i,2)
         tmp2(i,nyp2+1)=tmp2(i,3)
61     continue
      else 

       do 62 i=1,nxp2+2
cc              { solid wall in y}
         tmp1(i,3)=tmp1(i,3)+tmp1(i,2)
         tmp1(i,4)=tmp1(i,4)+tmp1(i,1)
         tmp1(i,nyp1)=tmp1(i,nyp1)+tmp1(i,nyp2+2)
         tmp1(i,nyp2)=tmp1(i,nyp2)+tmp1(i,nyp2+1)

         tmp2(i,2)=0.0
         tmp2(i,3)=tmp2(i,3)-tmp2(i,1)
         tmp2(i,nyp2)=0.0
         tmp2(i,nyp1)=tmp2(i,nyp1)-tmp2(i,nyp2+1)

62     continue

       do 63 i=1,nxp2+2
         tmp1(i,2)=tmp1(i,3)
         tmp1(i,1)=tmp1(i,4)
         tmp1(i,nyp2+2)=tmp1(i,nyp1)
         tmp1(i,nyp2+1)=tmp1(i,nyp2)
         tmp2(i,1)=-tmp2(i,3)
         tmp2(i,nyp2+1)=-tmp2(i,nyp1)

63     continue
      end if  

      if(ibd .eq. 0)then 
       do 70 j=1,nyp2+2
cc              { Periodic in x}
         tmp1(1,j)=tmp1(1,j)+tmp1(nxp1,j)
         tmp1(2,j)=tmp1(2,j)+tmp1(nxp2,j)
         tmp1(3,j)=tmp1(3,j)+tmp1(nxp2+1,j)

         tmp2(3,j)=tmp2(3,j)+tmp2(nxp2+1,j)
         tmp2(4,j)=tmp2(4,j)+tmp2(nxp2+2,j)
         tmp2(nxp1,j)=tmp2(nxp1,j)+tmp2(1,j)
         tmp2(nxp2,j)=tmp2(nxp2,j)+tmp2(2,j)
70     continue

       do 71 j=1,nyp2+2
         tmp1(nxp1,j)=tmp1(1,j)
         tmp1(nxp2,j)=tmp1(2,j)
         tmp1(nxp2+1,j)=tmp1(3,j)

         tmp2(1,j)=tmp2(nxp1,j)
         tmp2(2,j)=tmp2(nxp2,j)
         tmp2(nxp2+1,j)=tmp2(3,j)
         tmp2(nxp2+2,j)=tmp2(4,j)
71     continue
      else 

       do 72 j=1,nyp2+2
cc              { solid wall in x}
         tmp1(2,j)=0.0
         tmp1(3,j)=tmp1(3,j)-tmp1(1,j)
         tmp1(nxp2,j)=0.0
         tmp1(nxp1,j)=tmp1(nxp1,j)-tmp1(nxp2+1,j)

         tmp2(3,j)=tmp2(3,j)+tmp2(2,j)
         tmp2(4,j)=tmp2(4,j)+tmp2(1,j)
         tmp2(nxp1,j)=tmp2(nxp1,j)+tmp2(nxp2+2,j)
         tmp2(nxp2,j)=tmp2(nxp2,j)+tmp2(nxp2+1,j)
72     continue

       do 73 j=1,nyp2+2
         tmp1(1,j)=-tmp1(3,j)
         tmp1(nxp2+1,j)=-tmp1(nxp1,j)

         tmp2(2,j)=tmp2(3,j)
         tmp2(1,j)=tmp2(4,j)
         tmp2(nxp2+2,j)=tmp2(nxp1,j)
         tmp2(nxp2+1,j)=tmp2(nxp2,j)
 73     continue
      end if        

      return
      end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccc subroutines for the front code ccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine finit(nfronts,fp,ip,pt,icp,ine,ptcon,
     1      elcon,bptcon,belcon,elprop,maxpt,maxel,am)
c a subroutine to initialize nsphr spheres by dividing the
c sphere into four ``strips'' and setting each strip separately. 

      dimension icp(maxel,2), ine(maxel,2), pt(maxpt,2)
      integer ptcon,elcon,bptcon,belcon,np,ffp,lfp,fep,ne,ffe,lfe,fee 
      dimension ptcon(maxpt),elcon(maxel),elprop(maxel,3)
      dimension bptcon(maxpt),belcon(maxel)
        dimension fp(20,6),ip(20)
      common/front/np,ffp,lfp,fep, ne,ffe,lfe,fee,amax,amin 
      common/fparam/xmv,mv,fxl,fyl 
      common/source/num_srce,fsource(20),isource_pt(20)


      pi=4.*ATAN(1.)
      write(*,*)' Type in number of fronts'
      read(*,*)nfronts
      
      num_srce=0
      nptot=0
      netot=0
      do 150 is=1,nfronts
      write(*,*)' type in type,rad,xc,yc and prop1, prop2 for a sphere'
      read(*,*)itp,radin,xc,yc,prop1,prop2
          ip(is)=itp
          fp(is,1)=radin
          fp(is,2)=xc
          fp(is,3)=yc
          fp(is,4)=prop1
          fp(is,5)=0.0
          fp(is,6)=0.0
  
      if(itp.eq.1)then
cc                         { itp=1 for a cylinder}
        nps=int(2.0*pi*radin/am) 
        write(*,*)'sphere: ',is,' number of pts: ',nps
        dth=2.0*pi/float(nps) 
        rad=radin 
        do 100 i=1,nps 
         th=dth*(float(i)-0.5) 
         pt(i+nptot,1)=xc + rad*cos(th)
         pt(i+nptot,2)=yc + rad*sin(th) 
100     continue 
      end if
c-------------------------------------------- 
      if(itp.eq.2)then
        nps=int(xc/am)
        write(*,*)'sphere: ',is,' number of pts: ',nps
        do 105 i=1,nps 
cc                           { a horizontal surface } 
         pt(i+nptot,1)=xc-am*(float(i)-0.5)
cc                          { yc is the mean elivation and rad is }
         pt(i+nptot,2)=yc+radin*cos(2.0*pi*pt(i+nptot,1)/xc)
cc                                        { the perturb. ampl.}
105     continue  
      end if
c--------------------------------------------------- 
      if(itp.eq.3)then
        nps=int(xc/am)
        write(*,*)'sphere: ',is,' number of pts: ',nps
        do 110 i=1,nps 
cc                           { a horizontal surface } 
         pt(i+nptot,1)=am*(float(i)-0.5)
cc                      { yc is the mean elivation and rad is }
         pt(i+nptot,2)=yc+radin*cos(2.0*pi*pt(i+nptot,1)/xc)
cc                                        { the perturb. ampl.}
110     continue  
      end if  

c---------------------------------------------------
      if(itp.eq.5)then
        nps=1
        write(*,*)'source: ',is,' number of pts: ',nps
c xc,yx: coordinates, r: strength
         pt(1+nptot,1)=xc
         pt(1+nptot,2)=yc 
         num_srce=num_srce+1
         fsource(num_srce)=radin
         isource_pt(num_srce)=1+nptot
        np=nptot+nps
        nptot=np
        go to 150
      end if
c---------------------------------------------------
 
       do 120 i=1,nps 
           icp(i+netot,1)=i    +nptot
           icp(i+netot,2)=i+1  +nptot
           ine(i+netot,1)=i-1  +netot
           ine(i+netot,2)=i+1  +netot 
             elprop(i+netot,1)=prop1 
             elprop(i+netot,2)=prop2
120    continue  
           icp(nps+netot,2)=1  +nptot
           ine(1+netot,1)=nps  +netot
           ine(nps+netot,2)=1  +netot

      ne=nps                   +netot 
      netot=ne
      np=nps                   +nptot 
      nptot=np
     
150   continue

      ffp=1
      lfp=np
      fep=np+1

      ffe=1
      lfe=ne
      fee=ne+1

c set connectivity

      do 200 i = 1,maxpt-1
       ptcon(i)=i+1
       bptcon(i+1)=i
200   continue
       ptcon(lfp)=1
       bptcon(1)=lfp

      do 220 i = 1,maxel-1
       elcon(i)=i+1
       belcon(i+1)=i
220   continue
       elcon(lfe)=1
       belcon(1)=lfe

      return
      end
c-----------------------------------------------------------------------c

      subroutine fprint(pt,ptcon,iptmp,icp,elcon,maxpt,maxel,
     1                                         fname,noutf,lm)
c subroutine to print frontfile suitable for movie.byu 
      integer ptcon,elcon,np,ffp,lfp,fep, ne,ffe,lfe,fee 
      dimension icp(maxel,2), pt(maxpt,2) 
      dimension ptcon(maxpt),elcon(maxel),iptmp(maxpt)  
      character*20 outfile, fname
      character*1 num(10)
      common/front/np,ffp,lfp,fep, ne,ffe,lfe,fee,amax,amin 
      common/fparam/xmv,mv,fxl,fyl
      data (num(i),i=1,10)/'0','1','2','3','4','5','6','7','8','9'/ 
      outfile=fname 
      m1=mod(noutf,10)+1
      m0=mod((noutf/10),10)+1
      outfile(lm:lm+3)='frd'
      outfile(lm+3:lm+4)='0'
      outfile(lm+4:lm+5)=num(m0)  
      outfile(lm+5:lm+6)=num(m1)
      open(8,file=outfile,status='new') 

        iii=1
        write(8,254) iii,np,ne,iii
        write(8,253) iii,ne
254     format(4I6)
253     format(2I6)

        k=ffp
        do 277 kk=1,np
          write(8,255) pt(k,1),pt(k,2)
255       format(' ',2e13.5)
          k=ptcon(k)
277     continue

        k=ffp
        do 290 kk=1,np
          iptmp(k)=kk
          k=ptcon(k)
290     continue

        k=ffe
        do 288 kk=1,ne
         i=icp(k,1)
         ii=icp(k,2)
         write(8,266) iptmp(i),-iptmp(ii)
266      format(' ',2I6)
         k=elcon(k)
288    continue
         close(8,status = 'keep')
         noutf=noutf+1
      return
      end 
c-----------------------------------------------------------------------c
      subroutine fcurv(cv,t1,t2,cnt,pt,icp,ine,ptcon,elcon,maxel,maxpt)

      dimension icp(maxel,2), ine(maxel,2), pt(maxpt,2) 
      dimension ptcon(maxpt),elcon(maxel),t1(maxel,2)
      dimension cv(maxel,2),cnt(maxel,2), t2(maxel,2)
      integer ptcon,elcon,np,ffp,lfp,fep, ne,ffe,lfe,fee
      common/front/np,ffp,lfp,fep, ne,ffe,lfe,fee,amax,amin
      common/fparam/xmv,mv,fxl,fyl

      k=ffe
      do 10 kk=1,ne 
       xx=(xmv+0.5)*fxl    
       yy=(xmv+0.5)*fyl
       p0x=pt(icp(k,1),1)
       p0y=pt(icp(k,1),2)
       p1x= p0x+amod((pt(icp(k       ,2),1)-p0x+xx),fxl)-0.5*fxl
       p1y= p0y+amod((pt(icp(k       ,2),2)-p0y+yy),fyl)-0.5*fyl
       pm1x=p0x+amod((pt(icp(ine(k,1),1),1)-p0x+xx),fxl)-0.5*fxl
       pm1y=p0y+amod((pt(icp(ine(k,1),1),2)-p0y+yy),fyl)-0.5*fyl
       pp1x=p0x+amod((pt(icp(ine(k,2),2),1)-p0x+xx),fxl)-0.5*fxl
       pp1y=p0y+amod((pt(icp(ine(k,2),2),2)-p0y+yy),fyl)-0.5*fyl
       cnt(k,1)=(-pm1x+9.0*p0x+9.0*p1x-pp1x)*0.0625
       cnt(k,2)=(-pm1y+9.0*p0y+9.0*p1y-pp1y)*0.0625

       s1x=(-2.0*pm1x-3.0*p0x+6.0*p1x-pp1x)/6.0
       s1y=(-2.0*pm1y-3.0*p0y+6.0*p1y-pp1y)/6.0
       s2x=(pm1x-6.0*p0x+3.0*p1x+2.0*pp1x)/6.0
       s2y=(pm1y-6.0*p0y+3.0*p1y+2.0*pp1y)/6.0
       s1=sqrt(s1x**2+s1y**2)
       s2=sqrt(s2x**2+s2y**2)
       t1(k,1)=s1x/s1
       t1(k,2)=s1y/s1  
       t2(k,1)=s2x/s2
       t2(k,2)=s2y/s2 
      k=elcon(k)
10    continue
            
      k=ffe
      do 20 kk=1,ne
       cv(k,1)=0.5*(t2(k,1)+t1(ine(k,2),1)-t1(k,1)-t2(ine(k,1),1) )
       cv(k,2)=0.5*(t2(k,2)+t1(ine(k,2),2)-t1(k,2)-t2(ine(k,1),2) ) 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c       ds=sqrt((pt(icp(k,2),1)-pt(icp(k,1),1))**2 +
c     1                (pt(icp(k,2),2)-pt(icp(k,1),2))**2) 
c       cu=sqrt( cv(k,1)**2 +cv(k,2)**2)
c       curv=cu/ds
c       xxh=0.5*(pt(icp(k,1),1)+pt(icp(k,2),1)) 
c       yyh=0.5*(pt(icp(k,1),2)+pt(icp(k,2),2))
c       write(*,100)k,cu,ds,curv,cnt(k,1),xxh,cnt(k,2),yyh 
c100    format(' CC: ',i5,7e13.5)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      k=elcon(k)
20    continue 

       return
       end                                 
c-----------------------------------------------------------------------c
      subroutine fregrid(pt,icp,ine,ptcon,elcon,bptcon,belcon,
     1                                      elprop,maxpt,maxel)
      integer ptcon,elcon,bptcon,belcon,np,ffp,lfp,fep,ne,ffe,lfe,fee 
      dimension icp(maxel,2), ine(maxel,2), pt(maxpt,2) 
      dimension ptcon(maxpt),elcon(maxel),elprop(maxel,3)
      dimension bptcon(maxpt),belcon(maxel)
      common/front/np,ffp,lfp,fep, ne,ffe,lfe,fee,amax,amin 
      common/fparam/xmv,mv,fxl,fyl

c first delete small elements:
11    m=ffe
      do 10 kk=1,ne
       p0x=pt(icp(m,1),1)
       p0y=pt(icp(m,1),2)
       dx=amod((pt(icp(m,2),1)-p0x+(xmv+0.5)*fxl),fxl)-0.5*fxl
       dy=amod((pt(icp(m,2),2)-p0y+(xmv+0.5)*fyl),fyl)-0.5*fyl
       ds=sqrt(dx**2+dy**2)

         if(ds .lt. amin)then
         call fdelet(m,pt,icp,ine,ptcon,elcon,bptcon,belcon,
     2                                   elprop,maxpt,maxel)
         go to 11
         end if 
        if(ds .gt. amax)then 
        call faddel(m,pt,icp,ine,ptcon,elcon,bptcon,belcon,
     2                                    elprop,maxpt,maxel)
         go to 11
         end if
      m=elcon(m)
10    continue

c check for roughness
c     call fsmooth(pt,icp,elcon,maxpt,maxel,ine)        { Not done here}

      return
      end
c-----------------------------------------------------------------------c
      subroutine fdelet(m,pt,icp,ine,ptcon,elcon,bptcon,belcon,
     1                                        elprop,maxpt,maxel)
      integer ptcon,elcon,bptcon,belcon,np,ffp,lfp,fep,ne,ffe,lfe,fee 
      dimension icp(maxel,2), ine(maxel,2), pt(maxpt,2) 
      dimension ptcon(maxpt),elcon(maxel),elprop(maxel,3)
      dimension bptcon(maxpt),belcon(maxel)
      common/front/np,ffp,lfp,fep, ne,ffe,lfe,fee,amax,amin 
      common/fparam/xmv,mv,fxl,fyl
      integer p2

              write(*,1001)m
1001          format(' ****** delete an element, m: ',i5)
                                                           
c set the coordinates of the new point, incl. curvature effects.
       xx=(xmv+0.5)*fxl    
       yy=(xmv+0.5)*fyl
       p0x=pt(icp(m,1),1)
       p0y=pt(icp(m,1),2)
       p1x= p0x+amod((pt(icp(m       ,2),1)-p0x+xx),fxl)-0.5*fxl
       p1y= p0y+amod((pt(icp(m       ,2),2)-p0y+yy),fyl)-0.5*fyl
       pm1x=p0x+amod((pt(icp(ine(m,1),1),1)-p0x+xx),fxl)-0.5*fxl
       pm1y=p0y+amod((pt(icp(ine(m,1),1),2)-p0y+yy),fyl)-0.5*fyl
       pp1x=p0x+amod((pt(icp(ine(m,2),2),1)-p0x+xx),fxl)-0.5*fxl
       pp1y=p0y+amod((pt(icp(ine(m,2),2),2)-p0y+yy),fyl)-0.5*fyl

       pt(icp(m,1),1)=(-pm1x+9.0*p0x+9.0*p1x-pp1x)*0.0625
       pt(icp(m,1),2)=(-pm1y+9.0*p0y+9.0*p1y-pp1y)*0.0625
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c set element properties,
c       elprop(ine(m,1),1)=elprop(ine(m,1),1)+0.5*elprop(m,1) 
c       elprop(ine(m,2),1)=elprop(ine(m,2),1)+0.5*elprop(m,1) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c eliminate p2: 
      p2=icp(m,2)
      if(p2.eq.ffp)ffp=ptcon(p2)
      if(p2.eq.lfp)lfp=bptcon(p2)
      np=np-1
      ptcon(bptcon(p2))=ptcon(p2)
      bptcon(ptcon(p2))=bptcon(p2)
      bptcon(fep)=p2
      ptcon(p2)=fep
      fep=p2

c eliminate m 
      if(m.eq.ffe)ffe=elcon(m)
      if(m.eq.lfe)lfe=belcon(m)
      ne=ne-1
      elcon(belcon(m))=elcon(m)
      belcon(elcon(m))=belcon(m)
      belcon(fee)=m
      elcon(m)=fee
      fee=m
    
c set connections
      ine(ine(m,2),1)=ine(m,1)
      ine(ine(m,1),2)=ine(m,2)
      icp(ine(m,2),1)=icp(m,1)

      return
      end
c-----------------------------------------------------------------------c
      subroutine faddel(m,pt,icp,ine,ptcon,elcon,bptcon,belcon,
     1                                       elprop,maxpt,maxel)
      integer ptcon,elcon,bptcon,belcon,np,ffp,lfp,fep,ne,ffe,lfe,fee 
      dimension icp(maxel,2), ine(maxel,2), pt(maxpt,2) 
      dimension ptcon(maxpt),elcon(maxel),elprop(maxel,3)
      dimension bptcon(maxpt),belcon(maxel)
      common/front/np,ffp,lfp,fep, ne,ffe,lfe,fee,amax,amin 
      common/fparam/xmv,mv,fxl,fyl

              write(*,1001)m
1001          format(' ****** add an element, m: ',i5)

c set new point
      newpt=fep
      fep=ptcon(fep)
      ptcon(newpt)=ptcon(lfp)
      ptcon(lfp)=newpt
      bptcon(newpt)=lfp
      bptcon(ffp)=newpt
      lfp=newpt
      np=np+1

c set the coordinates of the new point, incl. curvature effects.
       xx=(xmv+0.5)*fxl    
       yy=(xmv+0.5)*fyl
       p0x=pt(icp(m,1),1)
       p0y=pt(icp(m,1),2)
       p1x= p0x+amod((pt(icp(m       ,2),1)-p0x+xx),fxl)-0.5*fxl
       p1y= p0y+amod((pt(icp(m       ,2),2)-p0y+yy),fyl)-0.5*fyl
       pm1x=p0x+amod((pt(icp(ine(m,1),1),1)-p0x+xx),fxl)-0.5*fxl
       pm1y=p0y+amod((pt(icp(ine(m,1),1),2)-p0y+yy),fyl)-0.5*fyl
       pp1x=p0x+amod((pt(icp(ine(m,2),2),1)-p0x+xx),fxl)-0.5*fxl
       pp1y=p0y+amod((pt(icp(ine(m,2),2),2)-p0y+yy),fyl)-0.5*fyl
       pt(newpt,1)=(-pm1x+9.0*p0x+9.0*p1x-pp1x)*0.0625
       pt(newpt,2)=(-pm1y+9.0*p0y+9.0*p1y-pp1y)*0.0625

c set new element: 
      newel=fee
      fee=elcon(fee)
      elcon(newel)=elcon(lfe)
      elcon(lfe)=newel
      belcon(newel)=lfe
      belcon(ffe)=newel
      lfe=newel
      ne=ne+1

      icp(newel,2)=icp(m,2)
      icp(newel,1)=newpt
      icp(m,2)=newpt

      ine(ine(m,2),1)=newel
      ine(newel,2)=ine(m,2)
      ine(newel,1)=m
      ine(m,2)=newel

c set element properties
c       elprop(newel,1)=0.5*elprop(m,1) 
c       elprop(m,1)=0.5*elprop(m,1)   
       elprop(newel,1)=elprop(m,1)
cc    { jump in density}
       elprop(newel,2)=elprop(m,2)
cc    { surface tension}
      return
      end
c-----------------------------------------------------------------------c
      subroutine fbackin(pt,icp,ine,ptcon,elcon,bptcon,belcon,
     1                                      elprop,maxpt,maxel)
      integer ptcon,elcon,bptcon,belcon,np,ffp,lfp,fep,ne,ffe,lfe,fee 
      dimension icp(maxel,2), ine(maxel,2), pt(maxpt,2) 
      dimension ptcon(maxpt),elcon(maxel),elprop(maxel,3)
      dimension bptcon(maxpt),belcon(maxel) 
      dimension itmppt(15000), itmpel(15000)
      character*10 infile
      common/front/np,ffp,lfp,fep, ne,ffe,lfe,fee,amax,amin 
      common/fparam/xmv,mv,fxl,fyl
         do 2 i=1,maxpt
2        itmppt(i)=0
         do 4 i=1,maxel
4        itmpel(i)=0

      write(*,100)
100   format(' type in name of the Front Backup file')
      read(*,200)infile
200   format(1a10)
      open(8,file=infile,status='old')
      read(8,*)np,ffp,lfp,fep, ne,ffe,lfe,fee
      lp=ffp
      do 10 kk=1,np
      read(8,*)lp,pt(lp,1),pt(lp,2),ptcon(lp),bptcon(lp) 
      itmppt(lp)=1
      lp=ptcon(lp)
10    continue
      le=ffe
      do 20 kk=1,ne
      read(8,*)le,icp(le,1),icp(le,2),
     1   ine(le,1),ine(le,2),elcon(le),belcon(le), 
     2   elprop(le,1),elprop(le,2)
      itmpel(le)=1
      le=elcon(le)
20    continue  
c establish a continous connectivity among the empty spaces for points and elements:
      iem=fep
      do 30 i=1,maxpt
      if(itmppt(i).eq.1)go to 30
      if(iem.eq.i)go to 30
            ptcon(iem)=i
            iem=i 
30    continue 
      iem=fee
      do 40 i=1,maxel 
      if(itmpel(i).eq.1)go to 40
      if(iem.eq.i)go to 40
            elcon(iem)=i
            iem=i 
40    continue

      close(8,status='keep')
      return
      end
c-----------------------------------------------------------------------c
      subroutine fbackout(pt,icp,ine,ptcon,elcon,bptcon,belcon,
     1                             elprop,maxpt,maxel,fname,n,lm)
      integer ptcon,elcon,bptcon,belcon,np,ffp,lfp,fep,ne,ffe,lfe,fee 
      dimension icp(maxel,2), ine(maxel,2), pt(maxpt,2) 
      dimension ptcon(maxpt),elcon(maxel),elprop(maxel,3)
      dimension bptcon(maxpt),belcon(maxel)
      character*20 outfile, fname
      character*1 num(10)
      common/front/np,ffp,lfp,fep, ne,ffe,lfe,fee,amax,amin
      common/fparam/xmv,mv,fxl,fyl
      data (num(i),i=1,10)/'0','1','2','3','4','5','6','7','8','9'/ 
      outfile=fname 
      m1=mod(n,10)+1
      m0=mod((n/10),10)+1
      outfile(lm:lm+3)='bkf'
      outfile(lm+3:lm+4)=num(m0) 
      outfile(lm+4:lm+5)=num(m1)
      open(8,file=outfile,status='new') 
      write(8,*)np,ffp,lfp,fep, ne,ffe,lfe,fee
      lp=ffp
      do 10 kk=1,np
      write(8,*)lp,pt(lp,1),pt(lp,2),ptcon(lp),bptcon(lp)
      lp=ptcon(lp)
10    continue
      le=ffe
      do 20 kk=1,ne
      write(8,*)le,icp(le,1),icp(le,2),
     1   ine(le,1),ine(le,2),elcon(le),belcon(le),
     2   elprop(le,1),elprop(le,2)
      le=elcon(le)
20    continue
      close(8,status='keep') 
      n=n+1
      return
      end
c-----------------------------------------------------------------------c
      subroutine fstatis(pt,icp,ine,ptcon,elcon,bptcon,belcon,
     1                        cv,elprop,maxpt,maxel,area,xc,yc)
      integer ptcon,elcon,bptcon,belcon,np,ffp,lfp,fep,ne,ffe,lfe,fee 
      dimension icp(maxel,2), ine(maxel,2), pt(maxpt,2) 
      dimension ptcon(maxpt),elcon(maxel),cv(maxel,2)
      dimension bptcon(maxpt),belcon(maxel), elprop(maxel,3)
      integer p1,p2
      common/front/np,ffp,lfp,fep, ne,ffe,lfe,fee,amax,amin
      common/fparam/xmv,mv,fxl,fyl  
      area=0.0 
c      cx=0.0
c      cy=0.0
      xm=0.0
      ym=0.0
      k=ffe 
      do 10 kk=1,ne
      k=elcon(k) 
      p1=icp(k,1)
      p2=icp(k,2)   
      area=area+0.25*((pt(p2,2)-pt(p1,2))*(pt(p2,1)+pt(p1,1))-
     1                (pt(p2,1)-pt(p1,1))*(pt(p2,2)+pt(p1,2)) )    

       dy=pt(p2,2)-pt(p1,2)
       dx=pt(p2,1)-pt(p1,1)
       xa=0.5*(pt(p2,1)+pt(p1,1))
       ya=0.5*(pt(p2,2)+pt(p1,2))
      xm=xm+0.5*dy*xa**2
      ym=ym-0.5*dx*ya**2

c      cx=cx+cv(k,1)
c      cy=cy+cv(k,2)
10    continue 
      xc=xm/area
      yc=ym/area
      write(*,*)'Area: ',area,' xc,yc: ',xc,yc
      return
      end
c-----------------------------------------------------------------------c