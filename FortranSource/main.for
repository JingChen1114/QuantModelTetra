c     main.for a program written by Dr Jing Chen and Dr Lindsey Leach 
c     Orthogonal Contrast Based Models for Quantitative Genetic Analysis in Autotetraploid Species
c     Jing Chen, Fengjun Zhang, Lin Wang, Lindsey Leach and Zewei Luo 

      INCLUDE 'link_fnl_shared.h'
      
      use lin_sol_gen_int
      
      parameter(maxind=1000,pi=3.1415926d0)
      
      implicit double precision(a-h,o-z)
      
      integer p1(4),p2(4),ioff(5)
      
      dimension trait(maxind),fre_off(5)
      dimension g_esti1(5),g_esti2(5),posterior(maxind,5)
      dimension g(5),f(5),w(4,5),x(5,5),g1(5,1),para(5,1),var(4)
      
      open(3,file='data_in_FT.txt',status='old')
      open(4,file='data_out.txt',status='old')
      open(5,file='LOD_scores.txt',status='old')
      open(6,file='AIC.txt',status='old')
      open(7,file='Likelihood.txt',status='old')
      open(8,file='Chi_square.txt',status='old')
      open(9,file='Ioff.txt',status='old')
      
      read(3,*) n
      do i=1,n
        read(3,*) trait(i)
      end do
      
      do i1=1,12  !for each parental genotype configuration
        do ia=1,50  !for a grid of values for the coefficient of DR 
          if((i1==9).and.(ia==34))then
            write(*,*) 'stop'
          end if
          call offspring(ia,i1,ioff,fre_off)        !generate offpsring population
          igenotype=0
          do i=1,5
            igenotype=igenotype+ioff(i)
          end do   
          sum=0.0d0        !Initial value for genotypic value
          do i=1,n
            sum=sum+trait(i)
          end do
          g_esti1(3)=sum/(n*1.0d0)
          g_esti1(1)=g_esti1(3)*0.6d0
          g_esti1(2)=g_esti1(3)*0.8d0
          g_esti1(4)=g_esti1(3)*1.2d0
          g_esti1(5)=g_esti1(3)*1.4d0
            
          sum=0.0d0        !Initial value for environmental variance
          do i=1,n
            sum=sum+(trait(i)-g_esti1(3))**2
          end do
          var_esti1=sqrt(sum/(n-1))
              
          slikelihood1=0.0d0
          do i=1,n
            s=0.0d0
            do k=1,5
              s=s+fre_off(k)*exp(-((trait(i)-g_esti1(k))**2)/(2.0d0
     &          *(var_esti1**2)))/(var_esti1*sqrt(2.0d0*pi))
            end do
            slikelihood1=slikelihood1+log(s)
          end do
              
          s=1.0d0
          ic=1
          do while(s>0.000001d0)
            do i=1,5
              do j=1,n
                p=0.0d0
                do k=1,5
                  p=p+((1.0d0/(sqrt(2.0d0*(var_esti1**2)*pi))*(exp(
     &              -((trait(j)-g_esti1(k))**2)/(2.0d0*(var_esti1**2
     &              ))))))*fre_off(k)
                end do
            
                posterior(j,i)=(((1.0d0/(sqrt(2.0d0*pi*(var_esti1**2
     &                          )))*(exp(-((trait(j)-g_esti1(i))**2)/
     &                         (2.0d0 *(var_esti1**2))))))*
     &                         fre_off(i))/p
              end do
            end do
        
            do i=1,5
              ss=0.0d0
              do j=1,n
                ss=ss+posterior(j,i)*trait(j)
              end do
              tt=0.0d0
              do j=1,n
                tt=tt+posterior(j,i)
              end do
              if(tt==0.0d0)then
                g_esti2(i)=0.0d0
              else
                g_esti2(i)=ss/tt
              end if
            end do
        
            p=0.0d0
            do i=1,5
              do j=1,n
                p=p+posterior(j,i)*((trait(j)-g_esti2(i))**2)
              end do
            end do
            var_esti2=sqrt(p/(n*1.0d0))
        
            slikelihood2=0.0d0
            do i=1,n
              ss=0.0d0
              do k=1,5
                ss=ss+fre_off(k)*exp(-((trait(i)-g_esti2(k))**2)/
     &             (2.0d0*(var_esti2**2)))/(var_esti2*sqrt(2.0d0*pi))
              end do
              slikelihood2=slikelihood2+log(ss)
            end do
        
            s=slikelihood2-slikelihood1
        
            slikelihood1=slikelihood2
        
            do i=1,5
              g_esti1(i)=g_esti2(i)
            end do
            var_esti1=var_esti2   
        
            write(*,*) ic
            ic=ic+1
      
          end do
           
          call chi_square(n,fre_off,trait,g_esti1,var_esti1,
     &                    ChiSquare)
          
          sum=0.0d0        !mean
          do i=1,n
            sum=sum+trait(i)
          end do
          smean=sum/(n*1.0d0)
            
          sum=0.0d0        !variance
          do i=1,n
            sum=sum+(trait(i)-smean)**2
          end do
          svar=sqrt(sum/(n-1))
              
          slikelihood=0.0d0
          do i=1,n
            s=0.0d0
            do k=1,5
              s=s+fre_off(k)*exp(-((trait(i)-smean)**2)/(2.0d0
     &          *(svar**2)))/(svar*sqrt(2.0d0*pi))
            end do
            slikelihood=slikelihood+log(s)
          end do
              
          sLOD=slikelihood1-slikelihood
              
          AIC=2*igenotype-2*slikelihood1
          
          if(igenotype==5)then    
            do i=1,5
              g(i)=g_esti1(i)
              f(i)=fre_off(i)
            end do
      
            call orthogonalscales(f,w)
      
            var=0.0d0
            do i=1,4
              s1=0.0d0
              do j=1,5
                s1=s1+f(j)*w(i,j)*w(i,j)
              end do
              s2=0.0d0
              do j=1,5
                s2=s2+f(j)*g(j)*w(i,j)
              end do
              var(i)=s2*s2/s1
            end do
      
            do i=1,5
              x(i,1)=1.0d0
            end do
            do i=1,4
              do j=1,5
                k=i+1
                x(j,k)=w(i,j)
              end do
            end do
            do i=1,5
              g1(i,1)=g(i)
            end do
            call lin_sol_gen(x,g1,para)
              
            write(4,*) "LOD score:", sLOD
            write(4,*) "AIC",AIC
            write(4,*) "Parental genotype:", i1
            call genotype(i1,p1,p2)
            write(4,*) "P1:", p1
            write(4,*) "P2:", p2
            write(4,*) "Coefficient of DR:", (ia-1)*0.005d0
            write(4,*) "Offspring:"
            write(4,*) ioff
            write(4,*) "Offspring probability distribution:" 
            write(4,*) fre_off
            write(4,*) "Estimated genetic effects:"
            do i=1,5
              write(4,*) para(i,1)
            end do
            write(4,*) "Estimated genetic variance:"
            do i=1,4
              write(4,*) var(i)
            end do
            write(4,*) "Estimated environmental variance:"
            write(4,*) var_esti1**2
            write(4,*) "***********************************"
            write(4,*) "***********************************"
            write(5,*) sLOD
            write(6,*) AIC
            write(7,*) slikelihood1
            write(8,*) ChiSquare
            noff=0
            do i=1,5
              noff=noff+ioff(i)
            end do
            write(9,*) noff
            
          else if(igenotype.ne.5)then
            call ReducedModel(g_esti1,fre_off,igenotype,para,var)
            write(4,*) "LOD score:", sLOD
            write(4,*) "AIC",AIC
            write(4,*) "Parental genotype:", i1
            call genotype(i1,p1,p2)
            write(4,*) "P1:", p1
            write(4,*) "P2:", p2
            write(4,*) "Coefficient of DR:", (ia-1)*0.005d0
            write(4,*) "Offspring:"
            write(4,*) ioff
            write(4,*) "Offspring probability distribution:" 
            write(4,*) fre_off
            write(4,*) "Estimated genetic effects:"
            do i=1,igenotype
              write(4,*) para(i,1)
            end do
            write(4,*) "Estimated genetic variance:"
            do i=1,igenotype-1
              write(4,*) var(i)
            end do
            write(4,*) "Estimated environmental variance:"
            write(4,*) var_esti1**2
            write(4,*) "***********************************"
            write(4,*) "***********************************"
            write(5,*) sLOD
            write(6,*) AIC
            write(7,*) slikelihood1
            write(8,*) ChiSquare
            noff=0
            do i=1,5
              noff=noff+ioff(i)
            end do
            write(9,*) noff
             
          end if !igenotype
        end do
          
      end do
         
      end      