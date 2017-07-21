      subroutine chi_square(n,fre_off,trait,g_esti,var_esti,
     &                      ChiSquare)
      
      parameter(maxind=1000,pi=3.1415926d0)
      
      implicit double precision(a-h,o-z)
      
      dimension trait(maxind),fre_off(5),g_esti(5)
      dimension expected(5),observed(5)
      dimension density(5),probability(5)
      
      do i=1,5
        expected(i)=n*fre_off(i)
      end do
      
      observed=0.0d0
      do i=1,n
        density=0.0d0
        probablitiy=0.0d0
        do j=1,5
          if(g_esti(j).ne.0.0d0)then
            density(j)=((1.0d0/(sqrt(2.0d0*(var_esti**2)*pi))*(exp(
     &              -((trait(i)-g_esti(j))**2)/(2.0d0*(var_esti**2
     &              ))))))
          end if
        end do
        s=0.0d0
        do j=1,5
          s=s+density(j)
        end do
        do j=1,5
          probability(j)=density(j)/s
        end do
        do j=1,5
          observed(j)=observed(j)+probability(j)
        end do
      end do
      ChiSquare=0.0d0
      do j=1,5
        if(expected(j).ne.0.0d0)then
          ChiSquare=ChiSquare+((observed(j)-expected(j))**2)/expected(j)
        end if
      end do
      
      return
      end
          
        