c       
c       F. Poitevin     October 2013
c       This tool orients the principal axis of a given structure
c       with the Z-axis, after centering its center of mass with
c       the origin.
c       
c       # Declare variables #
c       #
        include 'param.h'
        integer ierr,narg,nat,nch,ires(natmax),ch(natmax),info
        integer k
        real*8 crd(3*natmax),a(3,3),d(3),u(3,3),c(3)
        character*64 finput,foutput
        character aat(natmax)*3,rat(natmax)*3
        character lleft(natmax)*13,lright(natmax)*24
c       
c       # Define I/O formats #
c       #
1       format(" ",/,"ERROR while executing the program...")
2       format(" ",/,"... Usage : ./orient arg1 arg2",/,
     1         "          arg1: input PDB file name",/,
     2         "          arg2: output PDB file name",/,
     3         " ")
3       format(" ",/,"Happy ending !")
c       # Do some initializing #
c       #
        write(6,'(" ",/,"> Initializing...")')
        ierr = 0
c       
c       # Retrieve arguments #
c       #
        write(6,'(" ",/,"> Retrieving arguments...")')
        narg = iargc()
        if(narg.lt.2) then
          ierr = 1
          goto 999
        endif
        call getarg(1,finput)
        call getarg(2,foutput)
c       
c       # Store coordinates #
c       #
        write(6,'(" ",/,"> Reading input PDB file")')
        call readpdb(finput,nat,nch,aat,rat,ires,ch,crd,lleft,lright)
c       
c       # Center on center of mass #
c       #
        write(6,'(" ",/,"> Computing center of mass")')
        call com(crd,c,nat)
c       
c       # Compute matrix of inertia #
c       #
        write(6,'(" ",/,"> Computing matrix of inertia")')
        call maxinert(crd,crd,nat,a)
c       
c       # Diagonalize the matrix of inertia and retrieve the 
c       # rotation matrix (containing the eigenvectors) and 
c       # the principal moments of inertia (the eigenvalues)
        write(6,'(" ",/,"> Diagonalizing the matrix")')
        call jacobi(a,3,3,d,u,info)
        write(6,*) info
        write(6,'("DIAGO> The moments of inertia are:")')
        write(6,'(3(e12.5,1x))') d(1),d(2),d(3)
        write(6,'("DIAGO> i.e. the semi-principal diameters of the ",
     1  "ellipsoid are (after re-ordering) :")')
        call diagorder(3,d,u)
        write(6,'(3(e12.5,1x))') d(1),d(2),d(3)
        write(6,'("DIAGO> The corresponding inertial ellipsoid is ",
     1  "written as a PDB file in ellips.pdb (after inflation)")')
        call diagopdb(3,3,d,u)
c       
c       # Apply the rotation to align the structure with the
c       # reference frame
c       #
        write(6,'(" ",/,"> Rotating the centered structure:")')
        call rotate(crd,nat,u)
        write(6,'("ROTATE> writing output PDB file")')
        call writepdb(foutput,nat,nch,aat,rat,ires,ch,crd,lleft,lright)
c       
c       # End of the program #
c       # Give some insights if necessary #
c       #
999     continue
        if(ierr.ne.0) then
          write(6,1)
          if(ierr.eq.1) write(6,2)
        else
          write(6,3)
        endif
        stop
        end
