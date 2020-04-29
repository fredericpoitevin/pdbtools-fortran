c
	include 'param.h'
        integer iat,jat,k
	integer nat1,ires1(natmax)
        integer nat2,ires2(natmax)
	real*8 crd1(3*natmax),crd2(3*natmax)
        real*8 dist,dist_min
        character*64 fpdb1,fpdb2
	character lleft1(natmax)*22,lright1(natmax)*26
        character lleft2(natmax)*22,lright2(natmax)*26
c	
	call getarg(1,fpdb1)
	call getarg(2,fpdb2)
c	
	call readpdb_simple(fpdb1,nat1,ires1,crd1,lleft1,lright1)
	call readpdb_simple(fpdb2,nat2,ires2,crd2,lleft2,lright2)
c
	dist_min = 1.d6
	do 100 iat = 1,nat1
	  do 200 jat = 1,nat2
	    dist = 0.d0
	    do 300 k = 1,3
	      dist = dist
     1        + (crd1(3*(iat-1)+k) - crd2(3*(jat-1)+k))**2
300	    continue
	    dist = sqrt(dist)
	    if(dist.le.dist_min) then
	      dist_min = dist
	    endif
200	  continue
100	continue
c
        write(6,'(f6.2)') dist_min
c
999     continue
	stop
	end
c
