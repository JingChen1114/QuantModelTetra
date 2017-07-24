c     ***************************************************
c     subroutines used to generate gametes 
      subroutine gametogenesis1(p,ia,ib,g,ipair)
      
      implicit double precision(a-h,o-z)
      
      integer p(4),g(2),pair1(2),pair2(2)
      
      ipair=itger(1,3)
       
      if(ipair==1)then         !pairing pattern 1: 1,2&3,4
        pair1(1)=p(1)    
        pair1(2)=p(2)
        pair2(1)=p(3)
        pair2(2)=p(4)
      else if(ipair==2)then    !pairing pattern 2: 1,3&2,4
        pair1(1)=p(1)
        pair1(2)=p(3)
        pair2(1)=p(2)
        pair2(2)=p(4)
      else if(ipair==3)then    !pairing pattern 3: 1,4&2,3
        pair1(1)=p(1)
        pair1(2)=p(4)
        pair2(1)=p(2)
        pair2(2)=p(3)
      end if
      
      ix=itger(1,2)
      iy=itger(1,2)
      g(1)=pair1(ix)
      g(2)=pair2(iy)
      
      ia=ix
      ib=iy
      
      return
      end
c     ********************************************************
      subroutine gametogenesisr(r,p,ipair,ia,ib,g)
      
      implicit double precision(a-h,o-z)
      
      integer p(4),g(2),pair1(2),pair2(2)
      
      ja=ia
      ic1=ja
      do while(ic1.eq.ja)
        ic1=itger(1,2)
      end do    !while
      
      jb=ib
      ic2=jb
      do while(ic2.eq.jb)
        ic2=itger(1,2)
      end do    !while
      
      if(ipair==1)then         !pairing pattern 1: 1,2&3,4
        pair1(1)=p(1)    
        pair1(2)=p(2)
        pair2(1)=p(3)
        pair2(2)=p(4)
      else if(ipair==2)then    !pairing pattern 2: 1,3&2,4
        pair1(1)=p(1)
        pair1(2)=p(3)
        pair2(1)=p(2)
        pair2(2)=p(4)
      else if(ipair==3)then    !pairing pattern 3: 1,4&2,3
        pair1(1)=p(1)
        pair1(2)=p(4)
        pair2(1)=p(2)
        pair2(2)=p(3)
      end if
      
      f1=(1.0d0-r)**2
      f2=1.0d0-r**2
      call drnun(1,x1)
      if(x1.le.f1)then         !neither recombinant
        ix=ja
        iy=jb
      else if(x1.le.f2)then    !one recombinant
        call drnun(1,x2)
        if(x2.le.0.5d0)then
          ix=ja
          iy=ic2
        else
          ix=ic1
          iy=jb
        end if
      else
        ix=ic1
        iy=ic2
      end if
       
      g(1)=pair1(ix)
      g(2)=pair2(iy)
      
      ia=ix
      ib=iy
      
      return
      end
      
        