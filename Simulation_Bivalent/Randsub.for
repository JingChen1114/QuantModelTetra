c       randsub.for --- Fortran-77 procedures to generate random variate(s)
c                      of a given distribution. These include
c
c       (0). iseed---- to generate an unrepeated random seed;
c       (1). unif ---- to generate uniformly distributed random variate(s)
c                      over (0,1);
c       (2). gauss---- to generate normally distributed random variate(s)
c                      for given mean u and standard deviation sigma;
c       (3). itger---- to generate a random integer n1<=n<=n2;
c       (4). sequ.---- to generate an unrepeated random sequence of integers
c                      between n1 and n2 inclusive;

c     ***************
	function unif()

	implicit double precision (a-h,o-z)

	parameter(im1=2147483563,im2=2147483399,am=1.0d0/im1,imm1=im1-1,
     &     ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,
     &     ntab=32,ndiv=1+imm1/ntab,esp=1.2e-7,rnmx=1.0d0-esp)
	     
	integer iv(ntab)

	common /seed/ idum
	save iv,iy,idum2
	data idum2/123456789/,iv/ntab*0/,iy/0/

	unif=0.0d0

	if (idum.le.0) then
	  idum=max(-idum,1)
	  idum2=idum
	  do j=ntab+8,1,-1
	    k=idum/iq1
	    idum=ia1*(idum-k*ir1)
	    if (idum.lt.0) idum=idum+im1
	    if (j.le.ntab) iv(j)=idum
	  end do
	  iy=iv(1)
	end if
	k=idum/iq1
	idum=ia1*(idum-k*iq1)-k*ir1
	if (idum.lt.0) idum=idum+im1
	k=idum2/iq2
	idum2=ia2*(idum2-k*iq2)-k*ir2
	if (idum2.lt.0) idum2=idum2+im2
	j=1+iy/ndiv
	iy=iv(j)-idum2
	if (iy.lt.1) iy=iy+imm1
	unif=min(am*iy,rnmx)

	return
	end
c     ***********************
	function gauss(u,sigma)
	
	implicit double precision (a-h,o-z)
	
	common /seed/ idum
	save iset,gset
	data iset/0/

	if (iset.eq.0) then
	  rsq=0.0d0
	  do while (rsq.ge.1.0d0.or.rsq.eq.0.0d0)
	    v1=2.0d0*unif()-1.0d0
	    v2=2.0d0*unif()-1.0d0
	    rsq=v1**2+v2**2
	  end do !while
	  fac=sqrt(-2.0d0*log(rsq)/rsq)
	  gset=v1*fac
	  gauss=v2*fac
	  iset=1
	else
	  gauss=gset
	  iset=0
	end if
	gauss=sigma*gauss+u

	return
	end
c     *********************
	function itger(n1,n2)

c     generate a random integer n1<=iteger<=n2.       

	implicit double precision (a-h,o-z)
	parameter (maxind=200)
	dimension w(2*maxind)
	common /seed/ idum

	n=n2-n1+1
	w(1)=1.0d0/n  !size of interval
	do i=2,n
	  w(i)=w(i-1)+w(1)
	end do
	x=unif()
	i=1
	do while (x.gt.w(i))
	  i=i+1
	end do !while
	itger=n1+i-1

	return
	end
c     *****************************
	subroutine intsqce(n1,n2,ind)

c     generate an unrepeated integer sequence between n1 and n2 inclusive.
	
	implicit double precision (a-h,o-z)
	parameter (maxind=200)
	dimension ind(maxind),iseq(maxind)
	common /seed/ idum

	n=n2-n1+1
	do i=1,n
	  iseq(i)=i
	end do

	k=n
	do i=1,n
	  id=itger(1,k)
	  ind(i)=iseq(id)+n1-1
	  if (id.lt.k) then
	    do j=id,k-1
		iseq(j)=iseq(j+1)
	    end do
	  end if
	  k=n-i
	end do

	return
	end
c     ******************
	subroutine iseed()
	
	common /seed/ idum
	
	open(1,status='old',file='seed.dat')

	read(1,*) idum

	return
	end
c     ******************
	subroutine iexit()

	common /seed/ idum

	close(1,status='delete')
	open(2,status='new',file='seed.dat')

	newseed=-idum
	write(2,*) newseed
	write(2,*)

	return
	end