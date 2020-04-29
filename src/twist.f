c       
c       
        include 'param.h'
        integer k
        integer ierr,narg,nat,nch,ires(natmax),ch(natmax),info
        integer natmp
        integer iat,idx,ih,nh,iok
        real*8 crdmov(3*natmax),crdref(3*natmax),a(3,3),d(3),u(3,3)
        real*8 crdnew(3*natmax),commov(3),comref(3),rmsd,rotmat(3,3)
        real*8 tr,pi,dr
        real*8 zmin,zmax,zh(100),cc,nrm1,nrm2,x1(3),x2(3),un
        character*64 fmov,fref,foutput
        character aat(natmax)*3,rat(natmax)*3
        character lleft(natmax)*13,lright(natmax)*24
        character lleft_s(natmax)*22
c       
1       format(" ",/,"ERROR while executing the program...")
2       format(" ",/,"... Usage : ./twist arg1 arg2 arg3",/,
     1         "          arg1: to be moved PDB file name",/,
     2         "          arg2: reference PDB file name",/,
     3         "          arg3: output file name",/,
     4         " ")
3       format(" ",/,"Happy ending !")
4       format(" ",/,"... Not the same number of atoms in two PDBs")
c       
        write(6,'(" ",/,"> Initializing...")')
        ierr = 0
        pi=acos(-1.d0)
c       
c       # Retrieve arguments #
c       #
        write(6,'(" ",/,"> Retrieving arguments...")')
        narg = iargc()
        if(narg.lt.3) then
          ierr = 1
          goto 999
        endif
        call getarg(1,fmov)
        call getarg(2,fref)
        call getarg(3,foutput)
c       
c       # Store coordinates #
c       #
        write(6,'(" ",/,"> Reading input PDB files")')
        call readpdb_simple(fmov,nat,ires,crdmov,lleft_s,lright)
        natmp=nat
        call readpdb_simple(fref,nat,ires,crdref,lleft_s,lright)
        if(nat.ne.natmp) then
          ierr = 2
          goto 999
        endif
c       
        zmin=1.d6
        zmax=-zmin
        do 100 iat = 1,nat
          idx=3*(iat-1)+3
          if(crdref(idx).lt.zmin) zmin = crdref(idx)
          if(crdref(idx).gt.zmax) zmax = crdref(idx)
100     continue
c       
        nh=20
        do 200 ih = 1,nh
          zh(ih) = zmin + (zmax-zmin)*(real(ih)/real(nh))
          cc = 0.d0
          nrm1 = 0.d0
          nrm2 = 0.d0
          un = 0.d0
          do 201 iat = 1,nat
            do 202 k = 1,3
              idx=3*(iat-1)+k
              x1(k) = crdref(idx) 
              x2(k) = crdmov(idx)  
202         continue
            iok=0
            if(ih.eq.1) then
              if(x1(3).le.zh(ih)) iok=1
            else
              if(x1(3).gt.zh(ih-1).and.x1(3).le.zh(ih)) iok=1
            endif
            if(iok.eq.1) then
              do 203 k = 1,2
                cc    = cc + x1(k)*x2(k)
                nrm1 = nrm1 + x1(k)**2
                nrm2 = nrm2 + x2(k)**2
203           continue
              un = un + x1(1)*x2(2) - x1(2)*x2(1)
            endif
201       continue
          cc = cc/sqrt(nrm1*nrm2)
          un = un/sqrt(nrm1*nrm2)
          cc = 1.8d2*acos(cc)/pi
          if(un.lt.0.d0) cc = -cc
          write(6,'(2(e12.5,1x))') zh(ih),cc
200     continue
c       
c       # End of the program #
c       #
c       
999     continue
        if(ierr.ne.0) then
          write(6,1)
          if(ierr.eq.1) write(6,2)
          if(ierr.eq.2) write(6,4)
        else
          write(6,3)
        endif
        stop
        end
