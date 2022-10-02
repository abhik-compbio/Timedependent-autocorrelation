
	! code for calculation of Timedependent correlation function for protein dihedral angle
	! This Code can be used  only for auto correlation TDCF of dihedral angle.
	! Made by Abhik Ghosh Moulick

	implicit none
    real phip(150,80001),psip(150,80001),chi1p(150,80001),pi,c
	real corr1(0:80001),corr8(0:80001),corr15(0:80001),value,corr13(0:80001),corr14(0:80001)
    real corr2(0:80001),corr7(0:80001),corr3(0:80001),corr9(0:80001)
	real sigma_phi(150),sigma_psi(150),sigma_chi1(150)
	real av11,av12,av13
	real delt,norm1,norm8,norm15,norm2,norm3,norm7,norm9,norm13,norm14
	integer i,j,s,l,ll,M,p,t0,r1,nstep,steps,resnump,resnuml,chain,resnumb(20)
    integer kount(1000),kount1(1000),kount2(1000),jmax,del,pair
	character*20 resl(20),resp(1000)!,secstrp(150,80001),secstrl(150,80001)
    character*100 output1,output2
	data delt/10.0/,pi/3.14/,c/180.0/,value/0.000/     ! delt=8.0 ps 

    open(unit=7,file='dyncorr_respair.dat')
	read(7,*) r1
	do pair=1,1 !change
	read(7,*) output1,output2

	sigma_phi(r1)=0.00
	sigma_psi(r1)=0.00
	sigma_chi1(r1)=0.00
	
!	reading protein data
!   Input dataset contain dihedral angle value of all residue in each frame.

    open(unit=10,file='input-dataset.dat')!change
    resnump = 56    ! Total amino acid residue
    steps = 80000   ! Total frame number
	av11 = 0.0
	av12 = 0.0
	av13 = 0.0
	
	do i = 1,steps
	  do j = 1,resnump
	    read(10,1)nstep,resp(j),phip(j,i),psip(j,i),chi1p(j,i)
1	    format(i4,1x,a4,1x,3(f10.4,1x))
	    phip(j,i) = phip(j,i)*pi/c
	    psip(j,i) = psip(j,i)*pi/c
	    chi1p(j,i) = chi1p(j,i)*pi/c
	  enddo
	   av11 = av11 + sin(phip(r1,i))
	   av12 = av12 + sin(psip(r1,i))
	   av13 = av13 + sin(chi1p(r1,i))	   	    
    enddo
	av11 = av11/float(steps)
	av12 = av12/float(steps)
	av13 = av13/float(steps)
	
    close(10)

! protein-protein dihedral dynamic correlation
                                      
	open(unit=20,file=output1)
    open(unit=97, file=output2)


	write(*,*)pair
	do del=0,steps-1

	  M = steps-del
          corr1(del) = 0.0
		  corr2(del) = 0.0
          corr3(del) = 0.0
          corr7(del) = 0.0
          corr8(del) = 0.0
          corr9(del) = 0.0
          corr13(del) = 0.0
          corr14(del) = 0.0
          corr15(del)= 0.0
         
     do t0 = 1,M
		if (del.eq.0) then     !condition for normalization
		sigma_phi(r1)=sigma_phi(r1)+((sin(phip(r1,t0))-av11)**2)  
		sigma_psi(r1)=sigma_psi(r1)+((sin(psip(r1,t0))-av12)**2)
		sigma_chi1(r1)=sigma_chi1(r1)+((sin(chi1p(r1,t0))-av13)**2)
		end if

        corr1(del) = corr1(del) + (sin(phip(r1,t0))-av11)*(sin(phip(r1,t0+del))-av11)  !phi-phi autocorrelation
        corr2(del) = corr2(del) + (sin(phip(r1,t0))-av11)*(sin(psip(r1,t0+del))-av12)  !phi-psi autocorrelation
        corr3(del) = corr3(del) + (sin(phip(r1,t0))-av11)*(sin(chi1p(r1,t0+del))-av13) !phi-chi1 autocorrelation
        corr7(del) = corr7(del) + (sin(psip(r1,t0))-av12)*(sin(phip(r1,t0+del))-av11)  !psi-phi autocorrelation
        corr8(del) = corr8(del) + (sin(psip(r1,t0))-av12)*(sin(psip(r1,t0+del))-av12)  !psi-psi autocorrelation
        corr9(del) = corr9(del) + (sin(psip(r1,t0))-av12)*(sin(chi1p(r1,t0+del))-av13) !psi-chi1 autocorrelation
	    corr13(del) = corr13(del) + (sin(chi1p(r1,t0))-av13)*(sin(phip(r1,t0+del))-av11) !chi1-phi autocorrelation
        corr14(del) = corr14(del) + (sin(chi1p(r1,t0))-av13)*(sin(psip(r1,t0+del))-av12) !chi1-psi autocorrelation
        corr15(del) = corr15(del) + (sin(chi1p(r1,t0))-av13)*(sin(chi1p(r1,t0+del))-av13) !chi1-chi1 autocorrelation
	  enddo
	  
!       if(del.eq.0) then 
	    norm1 = sqrt(sigma_phi(r1)*sigma_phi(r1))
	    norm2 = sqrt(sigma_phi(r1)*sigma_psi(r1))
	    norm3 = sqrt(sigma_phi(r1)*sigma_chi1(r1))
	    norm7 = sqrt(sigma_psi(r1)*sigma_phi(r1))
	    norm8 = sqrt(sigma_psi(r1)*sigma_psi(r1))
	    norm9 = sqrt(sigma_psi(r1)*sigma_chi1(r1))
        norm13= sqrt(sigma_chi1(r1)*sigma_phi(r1))
	    norm14= sqrt(sigma_chi1(r1)*sigma_psi(r1))
	    norm15= sqrt(sigma_chi1(r1)*sigma_chi1(r1))
	   

	  corr1(del) = corr1(del)/norm1
      corr2(del) = corr2(del)/norm2
	  corr3(del) = corr3(del)/norm3
	  corr7(del) = corr7(del)/norm7
	  corr8(del) = corr8(del)/norm8
	  corr9(del) = corr9(del)/norm9
	  corr13(del)= corr13(del)/norm13
	  corr14(del)= corr14(del)/norm14
      corr15(del)= corr15(del)/norm15
	 

         

       write(20,4)float(del)*delt/1000.0,value,corr1(del),corr2(del),corr3(del),corr7(del),corr8(del),corr9(del)
	   write(97,4) float(del)*delt/1000.0,value,corr13(del),corr14(del),corr15(del)

              
4	  format(8(f10.4,'	'))  
	
      enddo
      close(20)
      close(97)
!        close(22)

      end do
       
      end
  
