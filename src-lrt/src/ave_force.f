c       
        integer i,ierr,ires(5,3000)
        real*8 force(5,3000),aveforce(3000),stdforce(3000)
        character line*100,ch*1,temp(3000)*54
        character*64 finput,fpdbout,ftxtout
c       
1       format(a)
2       format(22x,i4,40x,e12.5)
3       format("chain ",a1," is ",i1," with ",i4," atoms")
4       format(a54,2f5.2)
5       format(i4,2(1x,e12.5))
c       
        ierr = 0
        narg = iargc()
        if(narg.lt.3) then
          ierr = 1
          write(6,*) "ERROR"
          goto 999
        endif
        call getarg(1,finput)
        call getarg(2,fpdbout)
        call getarg(3,ftxtout)
        i=0
        ich=0
        ch='X'
        open(unit=1,file=finput,status='unknown')
100       read(1,1,end=200) line
          if(line(1:4).ne."ATOM") goto 100
          if(line(22:22).ne.ch) then
            write(6,3) ch,ich,i
            i=0
            ich=ich+1
            ch=line(22:22)
          endif
          i=i+1
          read(line,2) ires(ich,i),force(ich,i)
          if(line(22:22).eq."A") then
            temp(i) = line(1:54)
          endif
          goto 100
200       continue
        close(unit=1)
        write(6,3) ch,ich,i
        nres = i
        nch = ich
c       
        do 300 i = 1,nres
          aveforce(i) = 0.d0
          do 301 ich = 1,nch
            aveforce(i) = aveforce(i) + force(ich,i)
301       continue
          aveforce(i) = aveforce(i)/real(nch)
300     continue
c       
        do 302 i = 1,nres
          stdforce(i) = 0.d0
          do 303 ich = 1,nch
            stdforce(i) = stdforce(i)
     1                  + (force(ich,i)-aveforce(i))**2
303       continue
          stdforce(i) = sqrt(stdforce(i)/real(nch))
302     continue
c       
        open(unit=1,file=fpdbout,status='unknown')
        open(unit=2,file=ftxtout,status='unknown')
          do 400 i = 1,nres
            write(1,4) temp(i),stdforce(i),aveforce(i)
            write(2,5) i+4,aveforce(i),stdforce(i)
400       continue
        close(unit=1)
        close(unit=1)
c       
999     continue
        stop
        end
