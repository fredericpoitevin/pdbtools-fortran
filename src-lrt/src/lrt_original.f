c									
c	LRT	Frederic Poitevin	July 2013			
c									
c	This program compute the principal components of a matrix Q	
c	containing one structure per column (3N*m).			
c	More precisely, each column contains the deviation of one 	
c	structure to the mean structure.				
c	Instead of inverting the covariance matrix C=QQt,		
c	a singular value decomposition of Q is performed:		
c		Q = UWVt						
c	where W (m*m) is a diagonal matrix with singular values		
c	      U (3N*m) contains the corresponding eigenvectors 		
c	      V (m*m) is an orthogonal matrix                  		
c	The principal components are stored in U.			
c	One can easily show that:					
c		C   = QQt = UW(UW)t = U(W**2)Ut				
c		C-1 = U(W**-2)Ut					
c									
c	Thus, if one wants to solve Cx=b, one just does:		
c		x = U(W**-2)Ut b					
c	If b is the difference vector between the average resting	
c	state and the activated one, and C is known from the 		
c	resting state, thus x is the force vector, normalized		
c	by kT								
c									
c	
	include 'lrt.prm'
c	
	integer i,j,nmod,nat,m,mmax,equiv(nmodmax)
	integer ierr,nat2
c	
	real*8 crd(3*natmax,nmodmax),crdave(3*natmax)
	real*8 crdsave(3*natmax,nmodmax),v(nmodmax,nmodmax)
        real*8 crdforce(3*natmax)
	real*8 w(nmodmax),crdtar(3*natmax)
	real*8 uw(3*natmax,nmodmax),force(3*natmax)
	real*8 norm(natmax)
c	
	character fmodel*64,lleft(natmax)*30,lright(natmax)*12
	character ftarget*64
c	
	ierr = 0
c	
	write(6,*) " "
	write(6,'("Linear Response Theory ...")')
c									
c		STORE ALL MODELS AND BUILD Q (as crd minus average)	
c									
	write(6,*) " "
	call getarg(1,fmodel)
	call getarg(2,ftarget)
	call readmodels(fmodel,nmod,nat,crd,lleft,lright)
	write(6,*) " "
	call buildcrd(nmod,nat,crd,crdsave,crdave)
	call readtarget(ftarget,nat2,crdtar)
	if(nat2.ne.nat) then
	  ierr = 1
	  goto 999
	endif
c									
c		PERFORM THE SVD OF Q=UWVt, where Q is replaced by U	
c		and saved in crdsave					
c									
	write(6,*) " "
	write(6,'(". Now, performing the SVD (reordering from here)")')
	m = 3*nat
	mmax = 3*natmax
	call svdcmp(crd,m,nmod,mmax,nmodmax,w,v)
c									
c		REORDER (by decreasing singular values)			
c									
c	call reorder(m,nmod,crd,w,v,equiv)
c									
c		CHECK THAT U AND V ARE ORTHOGONAL MATRICES		
c									
	call orthocheck(m,nmod,crd,v)
c									
c	OUTPUTS PCA related						
c									
	write(6,'("...Writing singular values in singval.log")')
	open(unit=1,file="singval.log",status='unknown')
	open(unit=2,file="proj.log",status='unknown')
	write(1,'(17x,"RANK i",6x,"Wi")')
	write(2,'("# Proj. along first four PC")')
	do 100 i = 1,nmod
	  write(1,'("SINGULAR VALUE : ",i6,1x,e12.5)') i,w(i)
100	continue
	do 101 i = 1,nmod
	  write(2,'(i3,4(1x,e12.5))') i,w(1)*v(i,1),w(2)*v(i,2),
     1	    w(3)*v(i,3),w(4)*v(i,4)
101	continue
	close(unit=2)
	close(unit=1)
c	
	write(6,*) " "
	call writepcs(nmod,nat,crd,crdave,w,lleft,lright)
	write(6,*) " "
c									
c	OUTPUTS LRT related						
c									
c		Build UW**-1 (C**-1 = (UW**-1)(UW**-1)t			
c		Careful: close to zero singular values will not be	
c		inverted, but set to zero				
c									
	call builduw(crd,m,nmod,w,uw)
	call computef(uw,m,nmod,crdave,crdtar,force)
	call fbyres(force,norm,nat)
	call writepdb(crdave,norm,lleft,nat,1)
          do 200 i = 1,nat
            do 201 k = 1,3
              crdforce(3*(i-1)+k) = 
     1        crdave(3*(i-1)+k)+force(3*(i-1)+k)
201         continue
200       continue
        call writepdb(crdforce,norm,lleft,nat,2)
c	
c									
999	continue
	if(ierr.eq.1) then
	  write(6,'("ERROR: target has not the same number of atoms")')
	endif
	write(6,*) " "
	write(6,'("END OF THE PROGRAM")')
	write(6,*) " "
	stop
	end
c									
c									
c	
	subroutine fbyres(force,norm,nat)
c	
	include 'lrt.prm'
c	
	integer nat,i,k
c	
	real*8 force(3*natmax),norm(natmax)
c	
1	format(i5,1x,e12.5)
c	
	do 100 i = 1,nat
	  norm(i) = 0.d0
	  do 101 k = 1,3
	    norm(i) = norm(i) + force(3*(i-1)+k)**2
101	  continue
	  norm(i) = sqrt(norm(i))
	  write(6,1) i,norm(i)
100	continue
c	
	return
	end
c									
	subroutine computef(uw,m,nmod,crdave,crdtar,force)
c	
	include 'lrt.prm'
c	
	integer m,nmod,i,j,k
c	
	real*8 uw(3*natmax,nmodmax),crdave(3*natmax)
	real*8 crdtar(3*natmax),force(3*natmax)
	real*8 temp(nmodmax)
c	
	do 100 i = 1,nmod
	    temp(i) = 0.d0
	    do 102 k = 1,m
	      temp(i) = temp(i)
     1	        + uw(k,i)*(crdtar(k)-crdave(k))
102	    continue
100	continue
c	
	do 200 i = 1,m
	  force(i) = 0.d0
	  do 201 j = 1,nmod
	    force(i) = force(i)
     1	      + uw(i,j)*temp(j)
201	  continue
200	continue
c	
	return
	end
c									
	subroutine builduw(crd,m,nmod,w,uw)
c	
	include 'lrt.prm'
c	
	integer m,nmod,i,j,k
c	
	real*8 crd(3*natmax,nmodmax),w(nmodmax)
	real*8 uw(3*natmax,nmodmax)
	real*8 fac,thresh
c	
	thresh = 1.d-12
c	
	do 100 i = 1,m
	  do 101 j = 1,nmod
	    uw(i,j) = 0.d0
	    do 102 k = 1,nmod
	      if(w(k).lt.thresh) then
	        fac = 0.d0
	      else
	        fac = 1.d0/w(k)
	      endif
	      uw(i,j) = uw(i,j)
     1	        + crd(i,k)*fac
102	    continue
101	  continue
100	continue
c	
	return
	end
c									
	subroutine readtarget(ftarget,nat,crd)
c	
	include 'lrt.prm'
c	
	integer nat
c	
	real*8 crd(3*natmax)
c	
	character ftarget*64,line*80
c	
1	format(a)
2	format(30x,3(f8.3))
	nat = 0
	write(6,*) " "
	write(6,'(". Reading : ",a64)') ftarget
	open(unit=1,file=ftarget,status='unknown')
100	  read(1,1,end=200) line
	  if(line(1:4).ne."ATOM") goto 100
	  nat = nat + 1
	  read(line,2) crd(3*(nat-1)+1),crd(3*(nat-1)+2),
     1	  crd(3*(nat-1)+3)
	  goto 100
200	  continue
	close(unit=1)
c	
	write(6,'("... Number of atoms : ",i5)') nat
	return
	end
c									
c	
	subroutine writepdb(crd,bfac,lleft,nat,iflag)
c	
	include 'lrt.prm'
c	
	integer i,nat
c	
	real*8 crd(3*natmax),bfac(natmax),occ
c	
	character lleft(natmax)*30
c	
1	format(a30,3(f8.3),2(f6.2))
c	
	occ = 1.d0
        if(iflag.eq.1) then
	open(unit=1,file='force.pdb',status='unknown')
        else
        open(unit=1,file='vector.pdb',status='unknown')
        endif
	do 100 i = 1,nat
	  write(1,1) lleft(i),crd(3*(i-1)+1),crd(3*(i-1)+2),
     1	             crd(3*(i-1)+3),occ,bfac(i)
100	continue
	close(unit=1)
c	
	return
	end
c									
	subroutine readmodels(fmodel,nmod,nat,crd,lleft,lright)
c	
	include 'lrt.prm'
c	
	integer nmod,nat,natkeep
c	
	real*8 crd(3*natmax,nmodmax)
c	
	character*64 fmodel
	character line*100,lleft(natmax)*30,lright(natmax)*12
c	
1	format(a)
2	format(30x,3(f8.3))
3	format(a30,3(f8.3),a12)
c	
	write(6,'(". Reading : ",a64)') fmodel
	nat = 0
	nmod = 1
	open(unit=1,file=fmodel,status='unknown')
100	  read(1,1,end=200) line
	  if(line(1:4).eq."ATOM".and.line(14:15).eq."CA") then
	    nat = nat + 1
	    if(nmod.eq.1) then
	    read(line,3) lleft(nat),crd(3*(nat-1)+1,nmod),
     1	                 crd(3*(nat-1)+2,nmod),crd(3*(nat-1)+3,nmod),
     2	                 lright(nat)
	    else
	    read(line,2) crd(3*(nat-1)+1,nmod),crd(3*(nat-1)+2,nmod),
     1	                 crd(3*(nat-1)+3,nmod)
	    endif
	  elseif(line(1:6).eq."ENDMDL") then
	    natkeep = nat
	    nat = 0
	    nmod = nmod + 1
	  endif
	  goto 100
200	  continue
	close(unit=1)
	nmod = nmod - 1
	nat = natkeep
	write(6,'("...Number of models          : ",i3)') nmod
	write(6,'("...Number of atoms per model : ",i5)') nat
c	
	return
	end
c									
c	
	subroutine buildcrd(nmod,nat,crd,crdsave,crdave)
c	
	include 'lrt.prm'
c	
	integer nmod,nat,i,k,imod,jmod
c	
	real*8 crd(3*natmax,nmodmax),crdave(3*natmax)
	real*8 crdsave(3*natmax,nmodmax),rmsd(nmodmax)
	real*8 rmsdpair(nmodmax,nmodmax)
c	
	write(6,'(". Computing average and subtracting...")')
c	
	do 100 i = 1,nat
	  do 101 k = 1,3
	    crdave(3*(i-1)+k) = 0.d0
	    do 102 imod = 1,nmod
	      crdave(3*(i-1)+k) = crdave(3*(i-1)+k)
     1	        + crd(3*(i-1)+k,imod)
102	    continue
	    crdave(3*(i-1)+k) = crdave(3*(i-1)+k)/real(nmod)
	    do 103 imod = 1,nmod
	      crd(3*(i-1)+k,imod) = crd(3*(i-1)+k,imod)
     1	                          - crdave(3*(i-1)+k)
	      crdsave(3*(i-1)+k,imod) = crd(3*(i-1)+k,imod)
103	    continue
101	  continue
100	continue
c	
	
	write(6,'("...Writing deviation to ave for each in rmsd.log")')
	write(6,'("...Writing rmsd between pairs in rmsdpair.log")')
c	
	open(unit=1,file="rmsdpair.log",status='unknown')
	open(unit=2,file="rmsd2ave.log",status='unknown')
	do 200 imod = 1,nmod
	  rmsd(imod) = 0.d0
	  do 201 i =1,nat
	    do 202 k = 1,3
	    rmsd(imod) = rmsd(imod) 
     1	      + (crdsave(3*(i-1)+k,imod))*(crdsave(3*(i-1)+k,imod))
202	    continue
201	  continue
	  rmsd(imod) = sqrt(rmsd(imod)/real(nat))
	  write(2,'(i3,1x,e12.5,1x," Angstroem")') imod,rmsd(imod)
c	
	  do 203 jmod = 1,nmod
	    if(jmod.le.imod) goto 203
	    rmsdpair(imod,jmod) = 0.d0
	    do 204 i = 1,nat
	      do 205 k = 1,3
	        rmsdpair(imod,jmod) = rmsdpair(imod,jmod)
     1	         + (crdsave(3*(i-1)+k,imod)-crdsave(3*(i-1)+k,jmod))**2
205	      continue
204	    continue
	    rmsdpair(imod,jmod) = sqrt(rmsdpair(imod,jmod)/real(nat))
	    write(1,'(i3,1x,i3,1x,e12.5)') imod,jmod,rmsdpair(imod,jmod)
203	  continue
c	
200	continue
	close(unit=2)
	close(unit=1)
c	
	return
	end
c									
c	
	subroutine reorder(m,nmod,crd,w,v,equiv)
c	
	include 'lrt.prm'
c	
	integer m,nmod,i,j,l,equiv(nmodmax),equivtmp
c	
	real*8 crd(3*natmax,nmodmax),w(nmodmax),v(nmodmax,nmodmax)
	real*8 wtmp,vtmp,xtmp
c	
	do 99 i = 1,nmod
	  equiv(i) = i
99	continue
	do 100 i = 1,nmod
	  do 101 j = 1,nmod
	    if(j.le.i) goto 101
	    if(w(j).gt.w(i)) then
	      wtmp = w(i)
	      w(i) = w(j)
	      w(j) = wtmp
	      equivtmp = equiv(i)
	      equiv(i) = equiv(j)
	      equiv(j) = equivtmp
	      do 102 l = 1,nmod
	        vtmp = v(l,i)
	        v(l,i) = v(l,j)
	        v(l,j) = vtmp
102	      continue
	      do 103 l = 1,m
	        xtmp = crd(l,i)
	        crd(l,i) = crd(l,j)
	        crd(l,j) = xtmp
103	      continue
	    endif
101	  continue
100	continue
c	
	return
	end
c									
c	
	subroutine orthocheck(m,nmod,crd,v)
c	
	include 'lrt.prm'
c	
	integer m,nmod,i,j,k,imat,ierr
c	
	real*8 crd(3*natmax,nmodmax),v(nmodmax,nmodmax)
	real*8 matU(nmodmax,nmodmax),test,thresh
c	
1	format("...CHECK SVD OK: Matrix U is indeed orthogonal")
2	format("...CHECK SVD WRONG! Matrix U is NOT orthogonal...")
3	format("...CHECK SVD OK: Matrix V is indeed orthogonal")
4	format("...CHECK SVD WRONG! Matrix V is NOT orthogonal...")
c	
	thresh = 1.d-8
	ierr = 0
c	
	do 99 imat = 1,2
	do 100 i = 1,nmod
	  do 101 j = 1,nmod
	    matU(i,j) = 0.d0
	    if(imat.eq.1) then
	    do 102 k = 1,m
	      matU(i,j) = matU(i,j) + crd(k,i)*crd(k,j)
102	    continue
	    else
	    do 103 k = 1,nmod
	      matU(i,j) = matU(i,j) + v(k,i)*v(k,j)
103	    continue
	    endif
	    if(i.eq.j) then
	      test = abs(matU(i,j)-1.d0)
	      if(test.gt.thresh) ierr = 1
	    else
	      test = abs(matU(i,j))
	      if(test.gt.thresh) ierr = 1
	    endif
101	  continue
100	continue
c	
	if(ierr.eq.0.and.imat.eq.1) write(6,1)
	if(ierr.eq.0.and.imat.eq.2) write(6,3)
	if(ierr.eq.1.and.imat.eq.1) write(6,2)
	if(ierr.eq.1.and.imat.eq.2) write(6,4)
99	continue
c	
	return
	end
c									
c	
	subroutine writepcs(nmod,nat,crd,crdave,w,lleft,lright)
c	
	include 'lrt.prm'
c	
	integer i,imod,nmod,nat
c	
	real*8 crd(3*natmax,nmodmax),crdave(3*natmax)
	real*8 x,y,z,w(nmodmax)
c	
	character lleft(natmax)*30,lright(natmax)*12,fpc*7
c	
1	format(a30,3(f8.3),a12)
2	format("REMARK PRINCIPAL COMPONENT ",i3)
c	
	do 99 imod = 1,5
	  if(imod.eq.1) fpc = "pc1.pdb"
          if(imod.eq.2) fpc = "pc2.pdb"
          if(imod.eq.3) fpc = "pc3.pdb"
          if(imod.eq.4) fpc = "pc4.pdb"
          if(imod.eq.5) fpc = "pc5.pdb"
          open(unit=2,file=fpc,status='unknown')
	  write(2,2) imod
	  close(unit=2)
99	continue
c	
	write(6,'(". Writing average and 5 1st components as PDB files")')
	open(unit=1,file="ave.pdb",status='unknown')
	do 100 i = 1,nat
	  write(1,1) lleft(i),crdave(3*(i-1)+1),crdave(3*(i-1)+2),
     &	             crdave(3*(i-1)+3),lright(i)
	  do 101 imod = 1,5
	  if(imod.eq.1) fpc = "pc1.pdb"
	  if(imod.eq.2) fpc = "pc2.pdb"
	  if(imod.eq.3) fpc = "pc3.pdb"
	  if(imod.eq.4) fpc = "pc4.pdb"
	  if(imod.eq.5) fpc = "pc5.pdb"
	  open(unit=2,file=fpc,status='unknown',access='append')
	  x = crdave(3*(i-1)+1) + w(imod)*crd(3*(i-1)+1,imod)
	  y = crdave(3*(i-1)+2) + w(imod)*crd(3*(i-1)+2,imod)
	  z = crdave(3*(i-1)+3) + w(imod)*crd(3*(i-1)+3,imod)
	  write(2,1) lleft(i),x,y,z,lright(i)
	  close(unit=2)
101	  continue
100	continue
	close(unit=1)
c	
	return
	end
