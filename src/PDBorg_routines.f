c       
        subroutine reorganize(nat,nch,ires,ch,ierr)
c       
        include 'param.h'
c       
        integer i,ierr,iat,nat,nch,ires(natmax),ch(natmax)
        integer nchain,nres,ich,jch,kch,start_res,end_res
        integer iatmp,natmp,newres(natmax),newch(natmax)
        character*1 chname(nchmax),chtest,chout
        character*45 seqchainin,seqchainout
c       
        common /chains/ nchain,seqchainin,seqchainout
        common /resids/ nres,start_res,end_res
c       
        data (chname(i),i=1,61)/
     1  'A','B','C','D','E','F','G','H','I','J','K','L','M','N',
     2  'O','P','Q','R','S','T','U','V','W','X','Z','0','1','2',
     3  '3','4','5','6','7','8','9','a','b','c','d','e','f','g',
     4  'h','i','j','k','l','m','n','o','p','q','r','s','t','u',
     5  'v','w','x','y','z'/
c       
        do 100 iat = 1,nat
c       
          if(iat.eq.1) then
            iatmp = 1
            newres(iatmp) = start_res
          else
            iatmp = iatmp + 1
            if(ch(iat).ne.ch(iat-1)) then
              newres(iatmp) = start_res
            else
              if(ires(iat).ne.ires(iat-1)) then
                newres(iatmp) = newres(iatmp-1) + 1
                if(newres(iatmp).gt.end_res) then
                  iatmp = iatmp-1
                  goto 100
                endif
              else
                newres(iatmp) = newres(iatmp-1)
              endif
            endif
          endif
c       
          do 101 ich = 1,nchain
            chtest = seqchainin(ich:ich)
            chout = seqchainout(ich:ich)
            do 102 jch = 1,nchmax
              if(chtest.eq.chname(jch)) then
                if(ch(iat).eq.jch) then
                  do 103 kch = 1,nchmax
                    if(chout.eq.chname(kch)) then
                      newch(iat) = kch
                      goto 100
                    endif
103               continue  
                endif
              endif
102         continue
101       continue
c       
100     continue
c       
        natmp = iatmp
        if(natmp.ne.nat) then 
          ierr = 2
          goto 999
        endif
c       
        do 200 iat = 1,nat
          ires(iat) = newres(iat)
          ch(iat) = newch(iat)
200     continue
c
999     continue
        return
        end
c____________________________________________________________________
c              
        subroutine pdborg_readinput(finput)
c       
        include 'param.h'
c       
        integer i,nkeys,imin,imax,start_res,end_res
        integer nres,ich,nchain
        character keys(13)*35,param*45
        character line*80,keyword*35,value*45
        character*64 fpdbin,fpdbout,finput
        character*45 seqchainin,seqchainout
c       
        common /files/ fpdbin,fpdbout
        common /chains/ nchain,seqchainin,seqchainout
        common /resids/ nres,start_res,end_res
c       
1       format(a)
2       format("> INITIALIZING:",/,
     1  "  - retrieving arguments and parameters")
c       
        data nkeys /6/
        data (keys(i),i=1,6)/
     &  'Initial PDB filename              :',
     &  'Final   PDB filename              :',
     &  'Chains considered in initial PDB  :',
     &  'New order of chains               :',
     &  'Sequence number of 1st residue    :',
     &  'Sequence number of last residue   :'
     &  /
c       
        write(6,2)
        open(unit=1,file=finput,status='old')
10       read(1,1,end=20) line
c       
          param = '                              '
          if(line(1:1).eq.'#') goto 10
          keyword = line(1:35)
          value   = line(36:80)
c       
          do 30 i = 1,45
            if(value(i:i).eq.'!') then
              imax = i-1
              goto 40
            endif
30        continue
          imax = 45
40        continue
          do 50 i = imax,1,-1
            if(value(i:i).ne.' ') goto 60
50        continue
60        continue
          imax = i
          do 70 i = 1,imax
            if(value(i:i).ne.' ') goto 80
70        continue
80        continue
          imin = i
          param(1:imax-imin+1) = value(imin:imax)
          if(param.eq.'') goto 10
c       
          do 100 i = 1,nkeys
            if(keyword.eq.keys(i)) goto 110
100       continue
110       continue
c       
          if(i.eq.1) then
            read(param,*) fpdbin
          elseif(i.eq.2) then
            read(param,*) fpdbout
          elseif(i.eq.3) then
            read(param,*) seqchainin
            nchain = imax-imin+1
          elseif(i.eq.4) then
            read(param,*) seqchainout
          elseif(i.eq.5) then
            read(param,*) start_res
          elseif(i.eq.6) then
            read(param,*) end_res
          endif
c       
          goto 10
20        continue
c       
        nres = end_res-start_res+1
c       
        write(6,'("INIT> Reorganizing ",a15,$)') fpdbin
        write(6,'(" into ",a15)') fpdbout
        write(6,'("INIT> Reordering chains from ",$)') 
          do 300 ich = 1,nchain
            write(6,'(a1,$)') seqchainin(ich:ich)
300       continue
        write(6,'(" to ",$)') 
          do 301 ich = 1,nchain
            write(6,'(a1,$)') seqchainout(ich:ich)
301       continue
          write(6,'(" ")')
        write(6,'("INIT> Number of residues per chain ",i5,$)') nres
        write(6,'(" starting at ",i5,$)') start_res
        write(6,'(" ending at ",i5)') end_res
        write(6,'(" ")')
c       
        close(unit=1)
c       
        return
        end
