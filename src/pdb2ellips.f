c       
        include 'param.h'
c
        integer i,j
        integer nat,ires(natmax)
        real*8 crd(3*natmax),c(3),U(3,3),rg2,D
        character*64 finput,foutput
        character lleft_s(natmax)*22,lright(natmax)*24
c      
        call getarg(1,finput)
        call getarg(2,foutput)
c       
        call readpdb_simple(finput,nat,ires,crd,lleft_s,lright)
        call com(crd,c,nat)
c        call recenter(crd,nat,c)
        call maxinert(crd,crd,nat,U)
c
c       Update U from S matrix to M matrix, see:
c       An ellipsoidal approximation of protein shape
c       WR Taylor JM Thornton WG Turnell
c       https://doi.org/10.1016/0263-7855(83)80001-0       
c       First, rescale to actual diameter of the object
        call Dmax(crd,nat,D)
        call rescale_ellips(U,D)       
c        call StoM_ellips(U)
        call writepdb_com_anisou(foutput,c,U)
c        
        stop
        end
c       __________
        subroutine StoM_ellips(U)
c
        integer i,j
        real*8 U(3,3),rg2
c       
        rg2 = U(1,1) + U(2,2) + U(3,3)
        do 100 i = 1,3
          do 200 j = 1,3
            U(i,j) = -U(i,j)
200       continue
          U(i,i) = rg2 + U(i,i)
          write(6,*) U(i,1),U(i,2),U(i,3)
100     continue
        return
        end
c       __________
        subroutine rescale_ellips(U,D)
c       
        integer i,j
        real*8 U(3,3),rg2,D
c
        D = 1.d2*D**2
        rg2 = U(1,1) + U(2,2) + U(3,3)
        do 100 i = 1,3
          do 200 j = 1,3
            U(i,j) = D*U(i,j)/rg2
200       continue
100     continue
c       
        return
        end
