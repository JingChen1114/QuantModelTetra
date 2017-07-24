c     *******************************************************************************
c     subroutine to calculate the orthogonal contrast scales for offspring population
      INCLUDE 'link_fnl_shared.h'
      
      subroutine orthogonalscales(fre,w)
      use lin_sol_gen_int
      implicit double precision(a-h,o-z)
      
      dimension fre(5),w(4,5),wi(5,1)
      dimension x(5,5),y(5,1)
      
      u=0.0d0
	u=-4.0d0*fre(5)-3.0d0*fre(4)-2.0d0*fre(3)-fre(2)
      do i=1,5
	  w(1,i)=u+i*1.0d0-1.0d0
      end do                                    !scales for monogenic
      x=0.0d0
      x(1,1)=1.0d0
	x(1,2)=-2.0d0
	x(1,3)=1.0d0
	x(2,2)=1.0d0
	x(2,3)=-2.0d0
	x(2,4)=1.0d0
	x(3,3)=1.0d0
	x(3,4)=-2.0d0
	x(3,5)=1.0d0
	do i=1,5
	  x(4,i)=fre(i)
	end do
	do i=1,5
	  x(5,i)=w(1,i)*fre(i)
	end do

      do i=1,3
	  y(i,1)=1.0d0
	end do
	do i=4,5
	  y(i,1)=0.0d0
	end do

	call lin_sol_gen(x,y,wi)
	do i=1,5
	  w(2,i)=wi(i,1)
	end do                                    !scales for digenic

	do i=1,5
	  do j=1,5
	    x(i,j)=0.0d0
	  end do
	end do
      x(1,1)=-1.0d0
	x(1,2)=3.0d0
	x(1,3)=-3.0d0
	x(1,4)=1.0d0
	x(2,2)=-1.0d0
	x(2,3)=3.0d0
	x(2,4)=-3.0d0
	x(2,5)=1.0d0
	do i=1,5
	  x(3,i)=fre(i)
	end do
	do i=1,5
	  do j=1,2
	    x(j+3,i)=w(j,i)*fre(i)
	  end do
	end do

      do i=1,2
	  y(i,1)=1.0d0
	end do
	do i=3,5
	  y(i,1)=0.0d0
	end do

	call lin_sol_gen(x,y,wi)
	do i=1,5
	  w(3,i)=wi(i,1)
	end do                                    !scales for trigenic

      do i=1,5
	  do j=1,5
	    x(i,j)=0.0d0
	  end do
	end do
	x(1,1)=1.0d0
	x(1,2)=-4.0d0
	x(1,3)=6.0d0
	x(1,4)=-4.0d0
	x(1,5)=1.0d0
	do i=1,5
	  x(2,i)=fre(i)
	end do
	do i=1,5
	  do j=1,3
	    x(j+2,i)=w(j,i)*fre(i)
	  end do
	end do
 
	y(1,1)=1.0d0
	do i=2,5
	  y(i,1)=0.0d0
	end do

	call lin_sol_gen(x,y,wi)
	do i=1,5
	  w(4,i)=wi(i,1)
	end do                                    !scales for quadrigenic

      return
	end