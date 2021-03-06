c     Purpose:
c             to simulate cosegregation between marker loci and QTL in a
c             tetrasomic outbred population in a simulation study
c 
c     Records of revision:
c     Date                Programmer
c     28/10/2014          Jing Chen
      INCLUDE 'link_fnl_shared.h'
      
      parameter(maxind=1000,mxloci=20)
      
      implicit double precision(a-h,o-z)
      
      integer ftype(mxloci,4),mtype(mxloci,4)
      integer gt1(2),gt2(2),pt1(4),pt2(4)
      integer g(maxind,mxloci,4),o(maxind,mxloci,8)
      integer icount(5)
      
      dimension alpha(mxloci),rec(mxloci),rec1(mxloci)
      dimension fre_qtl(5),w(4,5),gene_para(4)
      dimension g_theo(5),trait(maxind)
      
      common /seed/ idum
      call iseed()
      
      open(3,file='data_in.txt',status='old')
      open(4,file='data_out.txt',status='old')
      open(5,file='info_mapping.txt',status='old')
      
      read(3,*) n,nloci,iqtl         !n: # of offspring individual
                                     !nloci: # of loci (marker+QTL)
                                     !iqtl: the location of QTL
      do i=1,nloci
        read(3,*) (ftype(i,j),j=1,4) !paternal genotype on all loci
      end do
      do i=1,nloci
        read(3,*) (mtype(i,j),j=1,4) !maternal genotype on all loci
      end do
      read(3,*) alpha(1)             !coefficient of DR on the first locus
      read(3,*) (rec(i),i=1,nloci)   !recombination frequencies between loci
      read(3,*) h                    !heritability
      read(3,*) gmean                !mean
      read(3,*) (gene_para(i),i=1,4) !monogenic,digenic,trigenic,quadrigenic effects
      do i=2,nloci  !calculate coefficient of DR for the remaining loci
        i1=i-1
        alpha(i)=(alpha(i1)*((3.0d0-4.0d0*rec(i))**2)+2.0d0*rec(i)*	  
     &           (3.0d0-2.0d0*rec(i)))/9.0d0
      end do    

c     frequencies distribution of offspring genoypes given parental genotypes 1122*1122
      fre_qtl(5)=(1.0d0/6.0d0+alpha(iqtl)/3.0d0)**2
	fre_qtl(4)=(1.0d0/6.0d0+alpha(iqtl)/3.0d0)*(2.0d0/3.0d0-2.0d0*
     &           alpha(iqtl)/3.0d0)*2.0d0
	fre_qtl(3)=2.0d0*((1.0d0/6.0d0+alpha(iqtl)/3.0d0)**2)+
     &           (2.0d0/3.0d0-2.0d0*alpha(iqtl)/3.0d0)**2
	fre_qtl(2)=(1.0d0/6.0d0+alpha(iqtl)/3.0d0)*(2.0d0/3.0d0-2.0d0* 
     &           alpha(iqtl)/3.0d0)*2.0d0
	fre_qtl(1)=(1.0d0/6.0d0+alpha(iqtl)/3.0d0)**2
	
	call orthogonalscales(fre_qtl,w)
	
	do i=1,5
	  u=gmean
	  do j=1,4
          u=u+gene_para(j)*w(j,i)
	  end do
	  g_theo(i)=u   !genotypic values
      end do
      
      write(*,*) "theoretical genotypic value:"
	write(4,*) "theoretical genotypic value:"
	
	do i=1,5
	  write(*,*) i,g_theo(i)
	  write(4,*) i,g_theo(i)
	end do
	
	g_var=0.0d0
	do i=1,5
	  g_var=g_var+((g_theo(i)-gmean)**2)*fre_qtl(i) !genetic variance
	end do
	ev=sqrt(g_var/h-g_var) !residual variance
	
	write(*,*) "environmental variance:"
	write(*,*) ev
	write(4,*) "environmental variance:"
	write(4,*) ev
	
c     generate gametes and offspring individuals	
	icount=0
	o=0
	g=0
	sum=0.0d0
	do i=1,n
	  do j=1,nloci
	    r=rec(j)
	    do k=1,4
	      pt1(k)=ftype(j,k)
	      pt2(k)=mtype(j,k)
	    end do
	    if(j==1)then          !the first locus
	      a=alpha(j)
	      call gametogensisa(a,pt1,ia1,ib1,gt1)
	      call gametogensisa(a,pt2,ia2,ib2,gt2)
	    else                  !the following loci
	      r=rec(j)
	      call gametogensisr(r,pt1,ia1,ib1,gt1)
	      call gametogensisr(r,pt2,ia2,ib2,gt2)
	    end if
	    
	    do k=1,2
	      g(i,j,k)=gt1(k)
	      g(i,j,k+2)=gt2(k)
	    end do
	    
	    do k=1,4
	      if(g(i,j,k).gt.0) o(i,j,g(i,j,k))=1    !from gtype to ptype
	    end do
	    if(j==iqtl)then
	      m=0
	      do k=1,4
	        m=m+g(i,j,k)
	      end do
	      s=gauss(0.0d0,ev)
	      if(m==8)then
	        trait(i)=g_theo(5)+s            !2222
	        icount(5)=icount(5)+1
	        
	      else if(m==7)then
	        trait(i)=g_theo(4)+s            !2221
	        icount(4)=icount(4)+1
	        
	      else if(m==6)then
	        trait(i)=g_theo(3)+s            !2211
	        icount(3)=icount(3)+1
	        
	      else if(m==5)then
	        trait(i)=g_theo(2)+s            !2111
	        icount(2)=icount(2)+1
	        
	      else if(m==4)then
	        trait(i)=g_theo(1)+s            !1111
	        icount(1)=icount(1)+1
	        sum=sum+trait(i)
	      end if 
	       
	    end if
	    
	  end do
	end do

c     print out the simulated data	
	write(5,*) n,nloci-1
	do i=1,nloci
	  if(i.ne.iqtl)then
    	write(5,*) (ftype(i,j),j=1,4)
      end if
	end do
	do i=1,nloci
	  if(i.ne.iqtl)then
   	    write(5,*) (mtype(i,j),j=1,4)
   	  end if
	end do
	do i=1,nloci
	  if(i.ne.iqtl)then
   	    write(5,*) alpha(i)
   	  end if
   	end do
   	ir=1
	do i=1,nloci
	  if(i<iqtl)then
   	    rec1(ir)=rec(i)
   	    ir=ir+1
   	  else if(i.eq.(iqtl+1))then
   	    ii=iqtl+1
   	    rec1(ir)=rec(iqtl)+rec(ii)-4.0d0*rec(iqtl)*rec(ii)/3.0d0
   	    ir=ir+1
   	  else if(i>(iqtl+1))then
   	    rec1(ir)=rec(i)
   	    ir=ir+1
   	  end if
   	end do
   	write(5,*) (rec1(i),i=1,nloci-1)
   	
	do i=1,n
	  do j=1,nloci
	    if(j.ne.iqtl)then
    	  write(5,*) (o(i,j,k),k=1,8)
    	end if
	  end do
	end do
	write(5,*) (trait(i),i=1,n)

	
	call iexit()
	
	end
	      