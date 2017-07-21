      subroutine gametogenesis(alpha,p,g,freg)
      
      implicit double precision(a-h,o-z)
      
      integer p(4),g(10,2)
      dimension freg(10)
      
      do i=1,4
        g(i,1)=p(i)
        g(i,2)=p(i)
        freg(i)=alpha/4.0d0
      end do
      
      g(5,1)=p(1)
      g(5,2)=p(2)
      freg(5)=(1.0d0-alpha)/6.0d0
      
      g(6,1)=p(1)
      g(6,2)=p(3)
      freg(6)=(1.0d0-alpha)/6.0d0
      
      g(7,1)=p(1)
      g(7,2)=p(4)
      freg(7)=(1.0d0-alpha)/6.0d0
      
      g(8,1)=p(2)
      g(8,2)=p(3)
      freg(8)=(1.0d0-alpha)/6.0d0
      
      g(9,1)=p(2)
      g(9,2)=p(4)
      freg(9)=(1.0d0-alpha)/6.0d0
      
      g(10,1)=p(3)
      g(10,2)=p(4)
      freg(10)=(1.0d0-alpha)/6.0d0
      
      return
      end
      