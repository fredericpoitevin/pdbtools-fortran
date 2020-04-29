c
        integer ires,iresold,iresnew,diff
        character*64 finput,foutput
        character line*80,ch*1,chold*1
        character*1 ires2*1,ires2old*1
        character lleft*21,lright*53
c       
1       format(a)
2       format(a21,a1,i4,a1,a53)
c
        call getarg(1,finput)
        call getarg(2,foutput)
c       
        chold='0'
        iresold=-1
        ires2old=' '
        open(unit=1,file=finput,status='unknown')
        open(unit=2,file=foutput,status='unknown')
100       read(1,1,end=200) line
           if(line(1:4).eq."ATOM") then
            read(line,2) lleft,ch,ires,ires2,lright
            if(ch.ne.chold) then
              chold=ch
              iresnew=ires
            else
              diff=ires-iresold
              if(diff.eq.0) then
                if(ires2.eq.' ') then
                  if(ires2old.ne.' ') then
                    iresnew=iresnew+1
                    ires2old=ires2
                  endif
                else
                  if(ires2.ne.ires2old) then
                    iresnew=iresnew+1
                    ires2old=ires2
                  endif
                endif
              else
                iresnew=iresnew+1
              endif
            endif
            iresold=ires
            write(6,*) ch,ires,ires2,iresnew
            write(2,2) lleft,ch,iresnew,' ',lright
          endif
          goto 100
200       continue
        close(unit=2)
        close(unit=1)
c
        stop
        end
