c       
c       
        include 'param.h'
        integer iat,istep,nstep,idx
        integer k
        integer ierr,narg,nat,nch,ires(natmax),ch(natmax),info
        integer natmp
        real*8 theta
        real*8 x(3),crd0(3*natmax),dcrd(3*natmax)
        real*8 crdmov(3*natmax),crdref(3*natmax),a(3,3),d(3),u(3,3)
        real*8 crdnew(3*natmax),commov(3),comref(3),rmsd,rotmat(3,3)
        real*8 tr,pi,dr
        character*64 fmov,fref,foutput
        character aat(natmax)*3,rat(natmax)*3
        character lleft(natmax)*13,lright(natmax)*24
        character lleft_s(natmax)*22
c       
1       format(" ",/,"ERROR while executing the program...")
2       format(" ",/,"... Usage : ./orient arg1 arg2",/,
     1         "          arg1: to be moved PDB file name",/,
     2         "          arg2: reference PDB file name",/,
     3         "          arg3: output PDB file name",/,
     4         " ")
3       format(" ",/,"Happy ending !")
4       format(" ",/,"... Not the same number of atoms in two PDBs")
5       format(a22,i4,4x,3f8.3,a26)
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
c        call readpdb(fmov,nat,nch,aat,rat,ires,ch,
c     1               crdmov,lleft,lright)
        natmp=nat
        call readpdb_simple(fref,nat,ires,crdref,lleft_s,lright)
c        call readpdb(fref,nat,nch,aat,rat,ires,ch,
c     1               crdref,lleft,lright)
        if(nat.ne.natmp) then
          ierr = 2
          goto 999
        endif
c
c        do 100 iat = 1,nat
c          do 101 k = 1,3
c            idx=3*(iat-1)+k
c            crd0(idx) = (crdmov(idx)+crdref(idx))/2.d0
c            dcrd(idx) = (crdmov(idx)-crdref(idx))/2.d0
c101       continue
c100     continue
c       
        nstep=20
        open(unit=1,file=foutput,status='unknown')
        do 200 istep = 1,nstep/2
          theta = 2.d0*pi*real(istep)/real(nstep)
          do 201 iat = 1,nat
            do 202 k = 1,3
              idx=3*(iat-1)+k
              x(k) = crdmov(idx)*(1+cos(theta))/2.d0
     1             + crdref(idx)*(1-cos(theta))/2.d0
c              x(k)= dcrd(idx)*cos(theta) + crd0(idx)
202         continue
            write(1,5) lleft_s(iat),ires(iat),
     1                   x(1),x(2),x(3),lright(iat)
201       continue
          write(1,'("ENDMDL")')
200     continue
        close(unit=1)
c              
c       # End of the program #
c       # Give some insights if necessary #
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
