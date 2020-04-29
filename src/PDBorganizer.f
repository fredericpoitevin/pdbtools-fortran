c       
c       F.Poitevin      November 2013
c       This programs allows to reorder the chains and renumber the
c       residues according to the user's wish
c       
c       INPUTS:
c       - PDB file to be re-organized
c       - input file listing arguments needed for the run
c       
c       OUTPUTS:
c       - PDB file re-organized
c       
c       **************************
c       *BEGINNING OF THE PROGRAM*
c       **************************
c       
c       -> Define variables
        include 'param.h'
c       
        integer ierr,narg
        integer nat,nch,ires(natmax),ch(natmax)
        real*8 crd(3*natmax)
        character*64 finput,fpdbin,fpdbout
        character aat(natmax)*3,rat(natmax)*3
        character lleft(natmax)*13,lright(natmax)*24
c       
        common /files/ fpdbin,fpdbout
c       
c       -> Define formats
1       format(" ",/," <<< PDBorganizer 1.0 >>>",
     1             /,"  F.Poitevin - Nov. 2013 ",
     2             /," ")
2       format(" ",/,"An error occured...")
3       format("USAGE ERROR> ./PDBorganizer 'input filename'")
4       format("PROBLEM while reorganizing. Contact me...")
c       
c       -> Collect arguments for the run and initialize
        write(6,1)
        ierr=0
        narg=iargc()
        if(narg.lt.1) then
          ierr=1
          goto 999
        endif
        call getarg(1,finput)
        call pdborg_readinput(finput)
c       
c       -> Read input PDB file
        call readpdb(fpdbin,nat,nch,aat,rat,ires,ch,crd,lleft,lright)
c       
c       -> Re-organize chain and residue numbering
        call reorganize(nat,nch,ires,ch,ierr)
        if(ierr.ne.0) goto 999
c       
c       -> Write output PDB file
        call writepdb(fpdbout,nat,nch,aat,rat,ires,ch,crd,lleft,lright)
c       
999     continue
c       ********************
c       *END OF THE PROGRAM*
c       ********************
c       -> do some troubleshooting, if necessary
        if(ierr.ne.0) then
          write(6,2)
          if(ierr.eq.1) write(6,3)
          if(ierr.eq.2) write(6,4)
        endif
        write(6,1)
        stop
        end
