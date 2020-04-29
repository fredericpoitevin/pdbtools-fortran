c
        subroutine writepdb(fpdb,nat,nch,aat,rat,ires,ch,crd,
     1  lleft,lright)
c       
        include 'param.h'
        integer i
        integer ierr,iat,nat,nch,ires(natmax),ch(natmax)
        real*8 crd(3*natmax)
        character fpdb*64,aat(natmax)*3,rat(natmax)*3
        character lleft(natmax)*13,lright(natmax)*24
        character chtest*1,chname(nchmax)*1
	data (chname(i),i=1,61)/
     1  'A','B','C','D','E','F','G','H','I','J','K','L','M','N',
     2  'O','P','Q','R','S','T','U','V','W','X','Z','0','1','2',
     3  '3','4','5','6','7','8','9','a','b','c','d','e','f','g',
     4  'h','i','j','k','l','m','n','o','p','q','r','s','t','u',
     5  'v','w','x','y','z'/
c       
1       format(a4,i7,2x,a3,1x,a3,1x,a1,i4,4x,3f8.3,a24)
2       format("WRITEPDB> writing to file ",a20)
c       
        write(6,2) fpdb
        open(unit=1,file=fpdb,status='unknown')
          do 100 iat = 1,nat
            chtest = chname(ch(iat))
            write(1,1) "ATOM",iat,aat(iat),rat(iat),chtest,ires(iat),
     1        crd(3*(iat-1)+1),crd(3*(iat-1)+2),crd(3*(iat-1)+3),
     2        lright(iat)
100       continue
        close(unit=1)
c       
        return
        end
c
c_________________________________________________________________
c       
        subroutine writepdb_simple(fpdb,nat,ires,crd,lleft_s,lright)
        include 'param.h'
c
        integer iat,nat,ires(natmax)
        real*8 crd(3*natmax)
        character fpdb*64,lleft_s(natmax)*22,lright(natmax)*24
c
1       format(a22,i4,4x,3f8.3,a24)
2       format("WRITEPDB> writing to file ",a20)
c
        write(6,2) fpdb
        open(unit=1,file=fpdb,status='unknown')
          do 100 iat = 1,nat
            write(1,1) lleft_s(iat),ires(iat),
     1        crd(3*(iat-1)+1),crd(3*(iat-1)+2),crd(3*(iat-1)+3),
     2        lright(iat)
100       continue
        close(unit=1)
c
        return
        end
c_________________________________________________________________
c       
        subroutine writepdb_com_anisou(fpdb,c,U)
        include 'param.h'
c
        integer i,j,V(3,3)
        real*8 c(3),U(3,3),Bfac
        character fpdb*64
c
1       format("WRITEPDB> writing to file ",a20)
c2       format("ATOM   ",i4,"  CA  GLY A",i4,4x,3f8.3,f6.3,f6.2)
2       format("ATOM   ",i4,"  CA  GLY A",i4,4x,3f8.3,f6.3,e12.5)
3       format("ANISOU ",i4,"  CA  GLY A",i4,2x,6(i7))
c
        Bfac = 0.d0
        do 100 i = 1,3
          Bfac = Bfac + U(i,i)
          do 200 j = 1,3
            V(i,j) = int(U(i,j))
200       continue
100     continue
        Bfac = Bfac
c       
        write(6,1) fpdb
        open(unit=1,file=fpdb,status='unknown')
          write(1,2) 1,1,c(1),c(2),c(3),1.d0,Bfac
          write(1,3) 1,1,V(1,1),V(2,2),V(3,3),V(1,2),V(1,3),V(2,3)
        close(unit=1)
c
        return
        end
c_________________________________________________________________
c       
        subroutine writepdbbfac(fpdb,nat,nch,aat,rat,ires,ch,crd,
     1  lleft,bfac)
c       
        include 'param.h'
        integer ierr,iat,nat,nch,ires(natmax),ch(natmax)
        integer iresiat,chiat
        integer i
        real*8 crd(3*natmax),bfac(nchmax,natmax)
        character fpdb*64,aat(natmax)*3,rat(natmax)*3
        character lleft(natmax)*13
        character chtest*1,chname(nchmax)*1
        data (chname(i),i=1,20)/
     1  'A','B','C','D','E','F','G','H','I','J','K','L','M','N',
     2  'O','P','Q','R','S','T'/
c       
c1       format(a13,a3,1x,a3,1x,a1,i4,4x,3f8.3,a24)
1       format(a13,a3,1x,a3,1x,a1,i4,4x,3f8.3,f6.2,f7.3)
2       format("WRITEPDB> writing to file ",a20)
c       
        write(6,2) fpdb
        open(unit=1,file=fpdb,status='unknown')
          do 100 iat = 1,nat
            chiat = ch(iat)
            iresiat = ires(iat)
            chtest = chname(ch(iat))
            write(1,1) lleft(iat),aat(iat),rat(iat),chtest,ires(iat),
     1        crd(3*(iat-1)+1),crd(3*(iat-1)+2),crd(3*(iat-1)+3),
     2        0.d0,bfac(chiat,iresiat)
100       continue
        close(unit=1)
c       
        return
        end
c_________________________________________________________________
c       
        subroutine writepdbbfac2(fpdb,nat,nch,aat,rat,ires,ch,crd,
     1  lleft,bfac)
c       
        include 'param.h'
        integer ierr,iat,nat,nch,ires(natmax),ch(natmax)
        integer iresiat,chiat,i
        real*8 crd(3*natmax),bfac(natmax)
        character fpdb*64,aat(natmax)*3,rat(natmax)*3
        character lleft(natmax)*13
        character chtest*1,chname(nchmax)*1
        data (chname(i),i=1,20)/
     1  'A','B','C','D','E','F','G','H','I','J','K','L','M','N',
     2  'O','P','Q','R','S','T'/
c       
c1       format(a13,a3,1x,a3,1x,a1,i4,4x,3f8.3,a24)
1       format(a13,a3,1x,a3,1x,a1,i4,4x,3f8.3,f6.2,f7.3)
2       format("WRITEPDB> writing to file ",a20)
c       
        write(6,2) fpdb
        open(unit=1,file=fpdb,status='unknown')
          do 100 iat = 1,nat
            chiat = ch(iat)
            iresiat = ires(iat)
            chtest = chname(ch(iat))
            write(1,1) lleft(iat),aat(iat),rat(iat),chtest,ires(iat),
     1        crd(3*(iat-1)+1),crd(3*(iat-1)+2),crd(3*(iat-1)+3),
     2        0.d0,bfac(iresiat)
100       continue
        close(unit=1)
c       
        return
        end
c_________________________________________________________________     
c             
        subroutine readpdb_simple(fpdb,nat,ires,crd,lleft_s,lright)
c       this subroutine is simpler than readpdb
        include 'param.h'
c
        integer nat,i_read_ok,ires(natmax)
        real*8 crd(3*natmax)
        character fpdb*64,line*80
        character lleft_s(natmax)*22,lright(natmax)*24
c
1       format(a)
2       format(a22,i4,4x,3f8.3,a26)
3       format("READPDB> number of atoms  ",i6)
        nat=0
        open(unit=1,file=fpdb,status='unknown')
100       read(1,1,end=200) line
          i_read_ok=0
          if(line(1:4).eq."ATOM") i_read_ok = 1
          if(i_read_ok.eq.0) goto 100
          if(line(17:17).eq." ".or.line(17:17).eq."A") i_read_ok = 1
          if(i_read_ok.eq.0) goto 100
          nat = nat + 1
          read(line,2) lleft_s(nat),ires(nat),
     1      crd(3*(nat-1)+1),crd(3*(nat-1)+2),crd(3*(nat-1)+3),
     2      lright(nat)
          goto 100
200       continue
        close(unit=1)
        write(6,3) nat
c
        return
        end 
c_________________________________________________________________
c 
        subroutine readpdb(fpdb,nat,nch,aat,rat,ires,ch,crd,
     1    lleft,lright)
c       
        include 'param.h'
c       
        integer ich
        integer i,nat,nch,ires(natmax),ch(natmax)
        real*8 crd(3*natmax)
        character chtest*1,chain*1,fpdb*64,aat(natmax)*3,line*80
        character rat(natmax)*3,lleft(natmax)*13,lright(natmax)*24
        character chname(nchmax)*1,alt*1
        data (chname(i),i=1,20)/
     1  'A','B','C','D','E','F','G','H','I','J','K','L','M','N',
     2  'O','P','Q','R','S','T'/
c       
1       format(a)
2       format(a13,a3,a1,a3,1x,a1,i4,4x,3f8.3,a24)
3       format("READPDB> number of atoms  ",i6)
4       format("READPDB> number of chains ",i6)
c
        nat = 0
        chain = 'X'
        nch = 0
        open(unit=1,file=fpdb,status='unknown')
100       read(1,1,end=200) line
          if(line(1:4).ne."ATOM") goto 100
          nat = nat + 1
          read(line,2) lleft(nat),aat(nat),alt,rat(nat),chtest,
     1      ires(nat),crd(3*(nat-1)+1),crd(3*(nat-1)+2),
     1      crd(3*(nat-1)+3),lright(nat)
          if(alt.eq.' '.or.alt.eq.'A') then
            goto 99
          else
            nat = nat - 1
            goto 100
          endif
99          continue
            if(chtest.ne.chain) then
              do 101 ich = 1,nchmax
                if(chtest.eq.chname(ich)) then
                  chain = chtest
                  goto 102
                endif
101           continue
102           continue
              nch = nch + 1
            endif
            ch(nat) = ich
          goto 100
200       continue
        close(unit=1)
c       
c       # Some output informations ...
        write(6,3) nat
        write(6,4) nch
c       
        return
        end
c       _________
        subroutine com(crd,c,nat)
c       
        include 'param.h'
c       
        integer k,iat,nat
        real*8 c(3),crd(3*natmax)
c       
1       format("COM> center of mass is:",/,3(e12.5,1x))
c       
        do 100 k = 1,3
          c(k) = 0.d0
100     continue
c       
        do 200 iat = 1,nat
          do 201 k = 1,3
            c(k) = c(k) + crd(3*(iat-1)+k)
201       continue
200     continue
c       
        do 300 k = 1,3
          c(k) = c(k)/real(nat)
300     continue
c       
        do 400 iat = 1,nat
          do 401 k = 1,3
            crd(3*(iat-1)+k) = crd(3*(iat-1)+k) - c(k)
401       continue
400     continue
c       
c       # Some output informations ...
        write(6,1) c(1),c(2),c(3)
c       
        return
        end
c       _________
        subroutine Dmax(crd,nat,D)
c
        include 'param.h'
c       
        integer i1,i2,k1,k2,nat
        real*8 crd(3*natmax),x1(3),x2(3)
        real*8 D,Dtmp
c
        D = 0.d0
        do 100 i1 = 1,nat-1
          do 101 k1 =  1,3
            x1(k1) = crd(3*(i1-1)+k1)
101       continue
          do 200 i2 = i1+1,nat
            Dtmp = 0.d0
            do 201 k2 = 1,3
              x2(k2) = crd(3*(i2-1)+k2)
              Dtmp = Dtmp
     1        + (x1(k1)-x2(k2))**2
201         continue
            if(Dtmp.gt.D) D = Dtmp
200       continue
100     continue
        D = sqrt(D)
c
        return
        end
c       _________
        subroutine recenter(crd,nat,com)
c       
        include 'param.h'
        integer k,iat,nat
        real*8 crd(3*natmax),com(3)
c       
        do 100 iat = 1,nat
          do 101 k = 1,3
            crd(3*(iat-1)+k) = crd(3*(iat-1)+k) + com(k)
101       continue
100     continue
c       
        return
        end
c       __________
        subroutine maxinert(crd1,crd2,nat,a)
c       
        include 'param.h'
c
        integer i,j,iat,nat
        real*8 crd1(3*natmax),crd2(3*natmax),x1,x2,a(3,3)
c       
1       format("INERTIA> The matrix of inertia is: ")
2       format(3(e12.5,1x))
c
        do 100 i = 1,3
          do 200 j = 1,3
            a(i,j) = 0.d0
            do 300 iat = 1,nat
              x1 = crd1(3*(iat-1)+i)
              x2 = crd2(3*(iat-1)+j)
              a(i,j) = a(i,j)
     1        + x1*x2
300         continue
200       continue
100     continue
c       
        write(6,1)
        do 400 i = 1,3
          do 500 j = 1,3
            a(i,j) = a(i,j)/real(nat)
500       continue
          write(6,2) a(i,1),a(i,2),a(i,3)
400     continue
c       
        return
        end
c       __________
        subroutine rotmat(iaxis,theta,u)
c       
        integer iaxis
        real*8 pi,theta,c,s,u(3,3)
c       
        pi = acos(-1.d0)
        theta = theta*pi/180.d0
        c = cos(theta)
        s = sin(theta)
c       
        if(iaxis.eq.1) then
          u(1,1) = 1.d0
          u(2,1) = 0.d0
          u(3,1) = 0.d0
          u(1,2) = 0.d0
          u(2,2) = c
          u(3,2) = s
          u(1,3) = 0.d0
          u(2,3) = -s
          u(3,3) = c
        elseif(iaxis.eq.2) then
          u(1,1) = c
          u(2,1) = 0.d0
          u(3,1) = -s
          u(1,2) = 0.d0
          u(2,2) = 1.d0
          u(3,2) = 0.d0
          u(1,3) = s
          u(2,3) = 0.d0
          u(3,3) = c
        else
          u(1,1) = c
          u(2,1) = s
          u(3,1) = 0.d0
          u(1,2) = -s
          u(2,2) = c
          u(3,2) = 0.d0
          u(1,3) = 0.d0
          u(2,3) = 0.d0
          u(3,3) = 1.d0
        endif
c       
        return
        end
