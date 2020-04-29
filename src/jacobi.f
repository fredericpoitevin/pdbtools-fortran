        subroutine jacobi(a,n,np,d,v,nrot)
c
c...    computes all eigenvalues and eigenvectors of a real symmetric matrix A,
c...    of size N and physical size Np;
c...    on output, elements of A above the diagonal are destroyed
c...    D returns the eigenvalues of A in its first N elements;
c...    V is a matrix whose columns contain the eigenvectors of A.
c...    Nrot returns the number of Jacobi rotations required.
c
c...    Taken directly from Numerical Recipes, p 347.
c
	integer nmax
        parameter(nmax=100)
c
	integer n,nrot,np,i,ip,iq,j
	real*8 sm,thres,g,h,tau,theta,s,c,t
        real*8  a(np,np),d(np),v(np,np),b(nmax),z(nmax)
c
        do 12 ip=1,n
           do 11 iq=1,n
              v(ip,iq)=0.
11         continue
           v(ip,ip)=1.
12      continue
c
        do 13 ip=1,n
           b(ip)=a(ip,ip)
           d(ip)=b(ip)
           z(ip)=0.
13      continue
c
        nrot=0
        do 24 i=1,50
           sm=0.
           do 15 ip=1,n-1
                do 14 iq=ip+1,n
                   sm=sm+abs(a(ip,iq))
14              continue
15         continue
           if(sm.eq.0) goto 100
           if(i.lt.4) then
                thres=0.2*sm/n**2
           else
                thres=0.
           endif
c
c
           do 22 ip=1,n-1
                do 21 iq=ip+1,n
                   g=100.*abs(a(ip,iq))
c
c...    after 4 sweeps, skip the rotation if the off diagonal rotation is small
c
                   if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.
     1                  (abs(d(iq))+g.eq.abs(d(iq)))) then
                        a(ip,iq)=0.
                   elseif(abs(a(ip,iq)).gt.thres) then
                        h=d(iq)-d(ip)
                        if(abs(h)+g.eq.abs(h)) then
                           t=a(ip,iq)/h
                        else
                           theta=0.5*h/a(ip,iq)
                           t=1./(abs(theta)+sqrt(1.+theta**2))
                           if(theta.lt.0) t=-t
                        endif
c
c
                   c=1./sqrt(1+t**2)
                   s=t*c
                   tau=s/(1.+c)
                   h=t*a(ip,iq)
                   z(ip)=z(ip)-h
                   z(iq)=z(iq)+h
                   d(ip)=d(ip)-h
                   d(iq)=d(iq)+h
                   a(ip,iq)=0.
c
                   do 16 j=1,ip-1
                        g=a(j,ip)
                        h=a(j,iq)
                        a(j,ip)=g-s*(h+g*tau)
                        a(j,iq)=h+s*(g-h*tau)
16                 continue
c
                   do 17 j=ip+1,iq-1
                        g=a(ip,j)
                        h=a(j,iq)
                        a(ip,j)=g-s*(h+g*tau)
                        a(j,iq)=h+s*(g-h*tau)
17                 continue
c
c
                   do 18 j=iq+1,n
                        g=a(ip,j)
                        h=a(iq,j)
                        a(ip,j)=g-s*(h+g*tau)
                        a(iq,j)=h+s*(g-h*tau)
18                 continue
c
                   do 19 j=1,n
                        g=v(j,ip)
                        h=v(j,iq)
                        v(j,ip)=g-s*(h+g*tau)
                        v(j,iq)=h+s*(g-h*tau)
19                 continue
c
                   nrot=nrot+1
c
                   endif
21              continue
22         continue
c
           do 23 ip=1,n
                b(ip)=b(ip)+z(ip)
                d(ip)=b(ip)
                z(ip)=0.
23         continue
c
24      continue
c
100     continue
c
        return
        end
c       
        subroutine showrotaxis(u,comref,commov)
c       
        include 'param.h'
        integer ip,np,k
        real*8 u(3,3),comref(3),crd(3),crdmin(3),commov(3)
c       
1       format("ATOM",4x,i3,"  CA  VAL A ",i3,4x,3f8.3)
c       
        np=100
        do 99 k = 1,3
          crdmin(k) = comref(k)-(real(np)/2.d0)*u(k,3)
99      continue
c       
        open(unit=1,file='comaxis.pdb',status='unknown')
        write(1,1) 1,1,commov(1),commov(2),commov(3)
        write(1,1) 2,2,comref(1),comref(2),comref(3)
        do 100 ip = 1,np
          do 101 k = 1,3
            crd(k) = crdmin(k)
     1             + (real(ip)-1.d0)*u(k,3)
101       continue
          write(1,1) ip+2,ip+2,crd(1),crd(2),crd(3)
100     continue
        close(unit=1)
c       
        return
        end
c       
        subroutine diagopdb(n,np,d,u)
c       
        integer n,np,k
        real*8 d(np),u(np,np),d_old(np)
c       
1       format("ATOM",6x,i1,"  CA  VAL A   ",i1,4x,3f8.3)
c       
        do 100 k = 1,np
          d_old(k) = d(k)
          d(k) = 1.d0/d(k)
100     continue
c       
        open(unit=1,file='ellips.pdb',status='unknown')
          write(1,1) 1,1,0.d0,0.d0,0.d0
          write(1,1) 2,2,d(1)*u(1,1),d(1)*u(2,1),d(1)*u(3,1)
          write(1,1) 3,3,d(2)*u(1,2),d(2)*u(2,2),d(2)*u(3,2)
          write(1,1) 4,4,d(3)*u(1,3),d(3)*u(2,3),d(3)*u(3,3)
        close(unit=1)
c       
        do 200 k = 1,np
          d(k) = d_old(k)
200     continue
c       
        return
        end
c      
        subroutine diagorder(np,d,u)
c       
        integer idir,np,kmax,kmin,k,j,kmid
        real*8 d(3),u(3,3),dkmax,dkmin
        real*8 c(3),v(3,3)
c       
        dkmax = -1.d3
        dkmin = -dkmax
        do 100 k = 1,3
          d(k) = 1.d0/sqrt(d(k))
          if(d(k).gt.dkmax) then
            dkmax = d(k)
            kmax = k
          endif
          if(d(k).lt.dkmin) then
            dkmin = d(k)
            kmin = k
          endif
100     continue
        do 200 k = 1,3
          if(k.eq.kmin.or.k.eq.kmax) goto 200
          kmid = k
200     continue
c       
        idir = 0
        if(kmax.eq.1) then
          if(kmin.eq.2) idir = 1
        elseif(kmax.eq.2) then
          if(kmin.eq.3) idir = 1
        elseif(kmax.eq.3) then
          if(kmin.eq.1) idir = 1
        endif
c       
        if(idir.eq.0) then
          c(1) = dkmax
          c(2) = d(kmid)
          c(3) = dkmin
        else
          c(1) = d(kmid)
          c(2) = dkmax
          c(3) = dkmin
        endif
        do 300 k = 1,3
          if(idir.eq.0) then
            v(k,1) = u(k,kmax)
            v(k,2) = u(k,kmid)
            v(k,3) = u(k,kmin)
          else
            v(k,1) = u(k,kmid)
            v(k,2) = u(k,kmax)
            v(k,3) = u(k,kmin)
          endif
300     continue
c       
        do 400 k = 1,3
          d(k) = c(k)
          do 401 j = 1,3
            u(j,k) = v(j,k)
401       continue
400     continue
c       
        return
        end
c       
        subroutine rotate(crd,nat,u)
c       
        include 'param.h'
c       
        integer k,iat,nat,j
        real*8 crd(3*natmax),u(3,3),rot(3*natmax)
c       
        do 100 iat = 1,nat
          do 101 k = 1,3
            rot(3*(iat-1)+k) = 0.d0
            do 102 j = 1,3
              rot(3*(iat-1)+k) = rot(3*(iat-1)+k)
     1          + crd(3*(iat-1)+j)*u(j,k)
102         continue
101       continue
100     continue
c       
        do 200 iat = 1,nat
          do 201 k = 1,3
            crd(3*(iat-1)+k) = rot(3*(iat-1)+k)
201       continue
200     continue
c       
        return
        end
