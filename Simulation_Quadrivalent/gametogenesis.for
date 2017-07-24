c	*************************************
	subroutine gametogensisa(a,p,ia,ib,g)
 
	implicit double precision (a-h,o-z)

	integer p(4),g(2)

	call drnun(1,x)
	if (x.le.a) then						!DR gamete
	  ix=itger(1,4)
	  iy=ix
	else									!Non-DR gamete
	  ix=itger(1,4)
	  iy=ix
	  do while (iy.eq.ix)
		iy=itger(1,4)
	  end do	!while
	end if

	g(1)=p(ix)
	g(2)=p(iy)

	ia=ix
	ib=iy

	return
	end
c	*************************************
	subroutine gametogensisr(r,p,ia,ib,g)

	implicit double precision (a-h,o-z)

	integer p(4),g(2)
	
	ja=ia
	ic1=ja
	do while (ic1.eq.ja)
	  ic1=itger(1,4)
	end do	!while

	jb=ib
	ic2=jb
	do while (ic2.eq.jb)
	  ic2=itger(1,4)
	end do	!while

	f1=(1.0d0-r)**2
	f2=1.0d0-r**2
	call drnun(1,x1)

	if(x1.le.f1) then						!neither recombinant				
	  ix=ja
	  iy=jb
	 
	else if (x1.le.f2) then					!one recombinant				
	  call drnun(1,x2)
	  if(x2.le.0.5d0) then			
		ix=ja
		iy=ic2
	  else
		ix=ic1
		iy=jb
	  end if

	else									!both recombinants
	  ix=ic1
	  iy=ic2
	  
	end if

	g(1)=p(ix)
	g(2)=p(iy)

	ia=ix
	ib=iy

	return
	end