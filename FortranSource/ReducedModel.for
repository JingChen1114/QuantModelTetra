      INCLUDE 'link_fnl_shared.h'
      
      subroutine ReducedModel(g_esti,fre_off,igenotype,para,var)
      
      use lin_sol_gen_int
      
      implicit double precision(a-h,o-z)
      
      dimension g_esti(5),fre_off(5),para(5,1),var(4)
      dimension x(5,5),w(4,5),y(5,1),wi(5,1)
      dimension genotype2(2,1),genotype3(3,1),genotype4(4,1)
      dimension fre2(2),fre3(3),fre4(4)
      dimension x2(2,2),x3(3,3),x4(4,4),y3(3,1),y4(4,1)
      dimension wi3(3,1),wi4(4,1)
      dimension para2(2,1),para3(3,1),para4(4,1)
      dimension orthogonal2(2,2),orthogonal3(3,3),orthogonal4(4,4)
      
      if(igenotype==2)then
        orthogonal2=1.0d0
        i1=1
        do i=1,5
          if(fre_off(i)>0.0d0)then
            fre2(i1)=fre_off(i)
            i1=i1+1
          end if
        end do
        u=-1.0d0*fre2(2)
        do i=1,2
	    orthogonal2(2,i)=u+i*1.0d0-1.0d0
        end do
	  ik=1
	  do i=1,5
	    if(fre_off(i)>0.0d0)then
	      genotype2(ik,1)=g_esti(i)
	      ik=ik+1
	    end if
	  end do
	  
	  do i=1,2
          x2(i,1)=1.0d0
        end do
        do i=2,2
          do j=1,2
            x2(j,i)=orthogonal2(i,j)
          end do
        end do
        
        call lin_sol_gen(x2,genotype2,para2)
        
	  para=0.0d0
	  do i=1,2
	    para(i,1)=para2(i,1)
	  end do 
        var=0.0d0
        do i=2,2
          s1=0.0d0
          do j=1,2
            s1=s1+fre2(j)*orthogonal2(i,j)*orthogonal2(i,j)
          end do
          s2=0.0d0
          do j=1,2
            s2=s2+fre2(j)*genotype2(j,1)*orthogonal2(i,j)
          end do
          var(i)=s2*s2/s1
        end do
        
      else if(igenotype==3)then
        orthogonal3=1.0d0
        i1=1
        do i=1,5
          if(fre_off(i)>0.0d0)then
            fre3(i1)=fre_off(i)
            i1=i1+1
          end if
        end do
        u=-2.0d0*fre3(3)-fre3(2)
        do i=1,3
          orthogonal3(2,i)=u+i*1.0d0-1.0d0
        end do
        x3=0.0d0
        x3(1,1)=1.0d0
	  x3(1,2)=-2.0d0
	  x3(1,3)=1.0d0
        do i=1,3
	    x3(2,i)=fre3(i)
	  end do
	  do i=1,3
	    x3(3,i)=orthogonal3(2,i)*fre3(i)
	  end do 
	  y3(1,1)=1.0d0
	  do i=2,3
	    y3(i,1)=0.0d0
	  end do
	  call lin_sol_gen(x3,y3,wi3)
	  do i=1,3
	    orthogonal3(3,i)=wi3(i,1)
	  end do                                    !scales for digenic
	  ik=1
	  do i=1,5
	    if(fre_off(i)>0.0d0)then
	      genotype3(ik,1)=g_esti(i)
	      ik=ik+1
	    end if
	  end do
	  
	  do i=1,3
          x3(i,1)=1.0d0
        end do
        do i=2,3
          do j=1,3
            x3(j,i)=orthogonal3(i,j)
          end do
        end do
        
        call lin_sol_gen(x3,genotype3,para3)
        
	  para=0.0d0
	  do i=1,3
	    para(i,1)=para3(i,1)
	  end do 
        var=0.0d0
        do i=2,3
          s1=0.0d0
          do j=1,3
            s1=s1+fre3(j)*orthogonal3(i,j)*orthogonal3(i,j)
          end do
          s2=0.0d0
          do j=1,3
            s2=s2+fre3(j)*genotype3(j,1)*orthogonal3(i,j)
          end do
          i1=i-1
          var(i1)=s2*s2/s1
        end do
      
      else if(igenotype==4)then
        orthogonal4=1.0d0
        i1=1
        do i=1,5
          if(fre_off(i)>0.0d0)then
            fre4(i1)=fre_off(i)
            i1=i1+1
          end if
        end do
        u=-3.0d0*fre4(4)-2.0d0*fre4(3)-fre4(2)
        do i=1,4
          orthogonal4(2,i)=u+i*1.0d0-1.0d0
        end do
        x4=0.0d0
        x4(1,1)=1.0d0
	  x4(1,2)=-2.0d0
	  x4(1,3)=1.0d0
	  x4(2,2)=1.0d0
	  x4(2,3)=-2.0d0
	  x4(2,4)=1.0d0
        do i=1,4
	    x4(3,i)=fre4(i)
	  end do
	  do i=1,4
	    x4(4,i)=orthogonal4(2,i)*fre4(i)
	  end do 
	  y4(1,1)=1.0d0
	  y4(2,1)=1.0d0
	  do i=3,4
	    y4(i,1)=0.0d0
	  end do
	  call lin_sol_gen(x4,y4,wi4)
	  do i=1,4
	    orthogonal4(3,i)=wi4(i,1)
	  end do                                    !scales for digenic
	  x4=0.0d0
	  x4(1,1)=-1.0d0
	  x4(1,2)=3.0d0
	  x4(1,3)=-3.0d0
	  x4(1,4)=1.0d0
	  do i=1,4
	    x4(2,i)=fre4(i)
	  end do
	  do i=1,4
	    do j=2,3
	      x4(j+1,i)=orthogonal4(j,i)*fre4(i)
	    end do
	  end do
	  y4(1,1)=1.0d0
	  do i=2,4
	    y4(i,1)=0.0d0
	  end do
	  call lin_sol_gen(x4,y4,wi4)
	  do i=1,4
	    orthogonal4(4,i)=wi4(i,1)
	  end do                                    !scales for trigenic
	  ik=1
	  do i=1,5
	    if(fre_off(i)>0.0d0)then
	      genotype4(ik,1)=g_esti(i)
	      ik=ik+1
	    end if
	  end do
	  
	  do i=1,4
          x4(i,1)=1.0d0
        end do
        do i=2,4
          do j=1,4
            x4(j,i)=orthogonal4(i,j)
          end do
        end do
        
        call lin_sol_gen(x4,genotype4,para4)
        
	  para=0.0d0
	  do i=1,4
	    para(i,1)=para4(i,1)
	  end do 
        var=0.0d0
        do i=2,4
          s1=0.0d0
          do j=1,4
            s1=s1+fre4(j)*orthogonal4(i,j)*orthogonal4(i,j)
          end do
          s2=0.0d0
          do j=1,4
            s2=s2+fre4(j)*genotype4(j,1)*orthogonal4(i,j)
          end do
          var(i)=s2*s2/s1
        end do
     
      
      end if
      
      end
      