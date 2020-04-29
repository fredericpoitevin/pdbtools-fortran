c       
        integer i,iline,x(3000)
        real*8 y1(3000),z1(3000),y2(3000)
        character line*100
c       
1       format(a)
2       format(i5,1x,e12.5)
3       format(i5,2(1x,e12.5))
4       format(i5,3(1x,e12.5))
c       
        iline = 0
        open(unit=1,file='aveforce.txt',status='unknown')
100       read(1,1,end=200) line
          iline = iline + 1
          read(line,3) x(iline),y1(iline),z1(iline)
          goto 100
200       continue
        close(unit=1)
c       
        iline = 0
        open(unit=1,file='avevector.txt',status='unknown')
101       read(1,1,end=201) line
          iline = iline + 1
          read(line,2) x(iline),y2(iline)
          goto 101
201       continue
        close(unit=1)
c       
        open(unit=1,file='catave.txt',status='unknown')
        do 300 i = 1,iline
          write(1,4) x(i),y1(i),z1(i),y2(i)
300     continue
        close(unit=1)
c       
        stop
        end
