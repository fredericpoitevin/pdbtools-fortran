c       
c       F. Poitevin     November 2013
c       This tool superimpose a structure onto a target
c       and outputs the corresponding trans-rot elements
c       
        include 'param.h'
        integer k
        integer ierr,narg,nat,nch,ires(natmax),ch(natmax),info
        integer natmp
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
c       # Center on center of mass #
c       #
        write(6,'(" ",/,"> Computing center of mass")')
        call com(crdmov,commov,nat)
        call com(crdref,comref,nat)
c       
c       # Now computing the rotation matrix
c       #
        write(6,'(" ",/,"> Superimposing to retrieve rotation matrix")')
        call compute(crdref,crdmov,crdnew,nat,rmsd,rotmat,ierr)
        write(6,'("ROTMAT> RMSD = ",e12.5," Angstroem")') rmsd
c       
c       # Now retrieve rotation axis and angle
c       #
        call jacobi(rotmat,3,3,d,u,info)
c       
        tr = d(1)+d(2)+d(3)
        tr = (tr-1.d0)/2.d0
        tr = acos(tr)
        tr = tr*180.d0/pi
        write(6,*) d(1),d(2),d(3)
        write(6,'("ROTMAT> ROTATION ANGLE = ",e12.5," degrees")') tr
c       
        call showrotaxis(u,comref,commov)
        dr = 0.d0
        do 100 k = 1,3
          dr = dr + (commov(k)-comref(k))**2
100     continue
        dr = sqrt(dr)
        write(6,'("ROTMAT> TRANSLATION = ",e12.5," Angstroem")') dr
c       
        write(6,'(" ",/,"> Writing output PDB file")')
        call recenter(crdnew,nat,comref)
        call writepdb_simple(foutput,nat,ires,crdnew,lleft_s,lright)
c        call writepdb(foutput,nat,nch,aat,rat,ires,ch,
c     1                crdnew,lleft,lright)
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
