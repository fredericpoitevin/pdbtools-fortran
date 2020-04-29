c       
c       F. Poitevin     October 2013
c       This tool rotates a structure by a given angle around either
c       X, Y or Z axis
c       
c       # Declare variables #
c       #
        include 'param.h'
        integer ierr,narg,nat,nch,ires(natmax),ch(natmax),info
        integer k,iaxis
        real*8 crd(3*natmax),theta,c(3)
        real*8 u(3,3)
        character*64 finput,foutput,tmparg
        character aat(natmax)*3,rat(natmax)*3,axis*1
        character lleft(natmax)*13,lright(natmax)*24
c       
c       # Define I/O formats #
c       #
1       format(" ",/,"ERROR while executing the program...")
2       format(" ",/,"... Usage : ./orient arg1 arg2 arg3 arg4",/,
     1         "          arg1: input PDB file name",/,
     2         "          arg2: output PDB file name",/,
     3         "          arg3: X, Y or Z",/,
     4         "          arg4: angle value (in degrees)",/,
     5         " ")
3       format(" ",/,"... arg3 must be either X, Y or Z")
4       format(" ",/,"Happy ending !")
c       
c       # Do some initializing #
c       #
        write(6,'(" ",/,"> Initializing...")')
        ierr = 0
c       
c       # Retrieve arguments #
c       #
        write(6,'(" ",/,"> Retrieving arguments...")')
        narg = iargc()
        if(narg.lt.4) then
          ierr = 1
          goto 999
        endif
        call getarg(1,finput)
        call getarg(2,foutput)
        call getarg(3,axis)
          if(axis.eq."X") then
            iaxis = 1
          elseif(axis.eq."Y") then
            iaxis = 2
          elseif(axis.eq."Z") then
            iaxis = 3
          else
            ierr = 2
          endif
        call getarg(4,tmparg)
          read(tmparg,*) theta
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
c       # Define rotation matrix #
c       #
        write(6,'(" ",/,"> Defining rotation matrix")')
        call rotmat(iaxis,theta,u)
c       
c       # Apply it to the centered structure
c       #
        write(6,'(" ",/,"> Apply it to the centered structure")')
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
          if(ierr.eq.2) write(6,3)
        else
          write(6,4)
        endif
        stop
        end
