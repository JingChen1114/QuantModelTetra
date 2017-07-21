c     contains subroutine offspring and subroutine genotype
c     subroutine genotype assigns parental genotype configurations (1) - (12)
      subroutine offspring(ia,i1,ioff,fre_off)
      
      implicit double precision(a-h,o-z)
      
      integer p1(4),p2(4),g1(10,2),g2(10,2)
      integer zygote(100,4),ioff(5)
      
      dimension fre_off(5),freg1(10),freg2(10)
      dimension frezygote(100)
      
      alpha=(ia-1)*0.005d0
      
      call genotype(i1,p1,p2)       !parental genotype
      call gametogenesis(alpha,p1,g1,freg1) !meiosis parent 1
      call gametogenesis(alpha,p2,g2,freg2) !meiosis parent 2
      
      fre_off=0.0d0
      do i=1,10
        do j=1,10
          k=(i-1)*10+j
          zygote(k,1)=g1(i,1)
          zygote(k,2)=g1(i,2)
          zygote(k,3)=g2(j,1)
          zygote(k,4)=g2(j,2)
          frezygote(k)=freg1(i)*freg2(j)
        end do
      end do
      do i=1,100
        isum=0
        do k=1,4
          isum=isum+zygote(i,k)
        end do
        if(isum==4)then
          fre_off(1)=fre_off(1)+frezygote(i)
        else if(isum==5)then
          fre_off(2)=fre_off(2)+frezygote(i)
        else if(isum==6)then
          fre_off(3)=fre_off(3)+frezygote(i)
        else if(isum==7)then
          fre_off(4)=fre_off(4)+frezygote(i)
        else if(isum==8)then
          fre_off(5)=fre_off(5)+frezygote(i)
        end if
      end do
      
      ioff=0
      do i=1,5
        if(fre_off(i)>0.0d0)then
          ioff(i)=1
        end if
      end do
          
      return
      end
      
c     *******************************************************************
      subroutine genotype(i,p1,p2)
c     subroutine to assign parental genotype configuration (1) - (12)
c     where A allele is represented by code 2 and a allele by code 1      
      implicit double precision(a-h,o-z)
      
      integer p1(4),p2(4)
      
      if(i==1)then
        p1(1)=1
        p1(2)=1
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=1
        p2(3)=1
        p2(4)=1
      else if(i==2)then
        p1(1)=1
        p1(2)=1
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=1
        p2(4)=1
      else if(i==3)then
        p1(1)=1
        p1(2)=1
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=2
        p2(4)=1
      else if(i==4)then
        p1(1)=2
        p1(2)=1
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=1
        p2(4)=1
      else if(i==5)then
        p1(1)=2
        p1(2)=1
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=2
        p2(4)=1
      else if(i==6)then
        p1(1)=2
        p1(2)=1
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=2
        p2(4)=2
      else if(i==7)then
        p1(1)=2
        p1(2)=2
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=2
        p2(4)=1
      else if(i==8)then
        p1(1)=2
        p1(2)=2
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=2
        p2(4)=2
      else if(i==9)then
        p1(1)=2
        p1(2)=2
        p1(3)=2 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=2
        p2(4)=2
      else if(i==10)then
        p1(1)=2
        p1(2)=1
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=1
        p2(3)=1
        p2(4)=1  
      else if(i==11)then
        p1(1)=2
        p1(2)=2
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=1
        p2(4)=1
      else if(i==12)then
        p1(1)=2
        p1(2)=2
        p1(3)=2 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=2
        p2(4)=1
      
      end if
      
      return
      end
      