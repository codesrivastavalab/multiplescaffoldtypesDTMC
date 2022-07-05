!===========================================================================================================================
!                                  $mcode                        MAIN PROGRAM
!===========================================================================================================================

        PROGRAM triangulation
        USE module_curvcalc ; USE module_rerun ; Use Module_dataformat ; USE module_datastruct ; USE module_mcsmoves
        Use module_compute_analytic_measures ; Use Module_Initialize_System ; Use module_linklist_calc

        IMPLICIT NONE  
        Include 'mpif.h'                                								            ! Header for using MPIF90

        Integer:: i,i1
        Integer :: ier,num_rings,nbuffer,buffsize
        Real(Kind=8),Dimension(:),Allocatable :: pern,Kparstr,Kperstr,cparstr,cperstr
        Real(Kind=8):: bud_area,kappa_unscaled
		Real(Kind=8), Allocatable :: surfaceprop(:,:),energyprop(:,:)
        Character(100) :: nname,fname,start_state,surffile,energyfile,syscmd,nem_profile

        fangle=-0.5 ;  maxtl=SQRT(3.0)
        pi=acos(-1.0000000000000) ; NPhase=6 ; nbuffer=10; buffsize=1

        Call MPI_Init(ier)
        Call MPI_Comm_Size(MPI_Comm_World,Nprocs,ier)
        Call MPI_Comm_Rank(MPI_Comm_World,mynod,ier)
        Nallgrp=MPI_Comm_World

        Write(nname,*),mynod
		Write(*,*)mynod,Nprocs

		Write(syscmd,*),'mkdir -pv rundir-'//Trim(Adjustl(nname))
		Call System(Trim(Adjustl(syscmd)))

        fname='./parameters.in'
        OPEN(01,File=Trim(Adjustl(fname)))                   		! Read the parameters from a  file 
        READ(01,*),kappa_unscaled,pr
        Allocate(pern(Nphase))
        READ(01,*),pern(1:Nphase)               			! Percentage of nematic field of kind A and B  
        Allocate(Kparstr(Nphase))
        READ(01,*),Kparstr(1:Nphase)               			! Kparallel flag for each phase 
        Allocate(Kperstr(Nphase))
        READ(01,*),Kperstr(1:Nphase)                       		! Kperp flag for each phase 
        Allocate(cparstr(Nphase))
        READ(01,*),cparstr(1:Nphase)     				! c0 parallel for each phase
        Allocate(cperstr(Nphase))
        READ(01,*),cperstr(1:Nphase)                   			! c0 perp for each phase
        Allocate(nempar(Nphase,Nphase))
        Allocate(unempar(Nphase,Nphase))
        Do i=1,Nphase,1
        READ(01,*),nempar(i:i,:)
        Enddo 
        Read(01,*),start_state                                   	! starting state 'MEMBRANE' or 'RESTART'                                
		Read(01,*),itime,ftime,dataint,confint
		Read(01,*),nem_profile,bud_area,num_rings	
        Close(01)

	    kappa = kappa_unscaled/2.0                                                         ! the form of kappa/2 used in the Helfrich 
	    Call Initialize_random_number_seed()	                                           ! Start from a random seed 
        
		Allocate(surfaceprop(nbuffer,4))
		Allocate(energyprop(nbuffer,4))


        If(Trim(start_state).Eq. 'MEMBRANE')Then
        Call Read_Memb_Configuration(start_state)            			! Reads only the membrane data
        !************
	mp%nadd=0
	do i=1,NPhase
        mp%nadd=mp%nadd + Nint(nver*pern(i)) 
	end do
  	!************                  		                ! No of additive particles to be present 
        mp%nmem=nver-mp%nadd                           			        ! No of membrane particle  to be present 
		If (Trim(Adjustl(nem_profile)) .Eq. 'CIRCULAR') Then
         Call Initialize_PatchOfProteins(mp%nadd)                      	        ! Arranges all nematics as a single patch  
		Else If (Trim(Adjustl(nem_profile)) .Eq. 'RANDOM') Then
         Call Initialize_Proteins_Randomly(mp%nadd)                     	! Randomly distribute the nematics       
		Else If (Trim(Adjustl(nem_profile)) .Eq. 'BUDSIZE_NUMBER') Then
		 Call Initialize_Annulus_Area_and_Number(bud_area,mp%nadd)
        Else If (Trim(Adjustl(nem_profile)) .Eq. 'BUDSIZE_RINGS') Then
  	  	 Call Initialize_Annulus_Area_and_rings(bud_area,num_rings)
		Else If (Trim(Adjustl(nem_profile)) .Eq. 'BUDSIZE_PATCH') Then
	 	 Call Initialize_Patch_Specified_Area(bud_area)
	 	Endif
        Else If(Trim(start_state).Eq. 'RESTART')Then
        Call Read_Memb_Configuration(start_state)                    								    ! Reads membrane + nematic data

  !************
	mp%nadd=0
	do i=1,NPhase
        mp%nadd=mp%nadd + Nint(nver*pern(i)) 
	end do
  !************                 	 								            ! No of additive particles to be present 
        mp%nmem=nver-mp%nadd
       Endif

        Do i=1,nver,1
        ver(i)%kpar=Kparstr(ver(i)%phase) ; ver(i)%kper=Kperstr(ver(i)%phase)                           ! Allocate Parameters to individual vertices 
        ver(i)%cpar=Cparstr(ver(i)%phase) ; ver(i)%cper=Cperstr(ver(i)%phase)  
        Enddo 

        Allocate(mp%vertex(nver)) ; Allocate(mp%triangle(ntr)) ; Allocate(mp%link(-tlink:tlink))	
		Call Initialize_linklist()

       vertice_Call: DO i=1,nver                               								            ! Normal calculation over each vertex 
       Call normalcalc(i)
       Enddo vertice_Call

       Call vtkformat(0)                                                                                ! Dumps a vtu file that may be viewed in Paraview
       Call Write_Memb_Configuration(0)                                                                 ! Dump a configuration that can be used to restart
       mvinter=Nint(Real(tlink)/nver)                                 								    ! Ratio of link flip to vertex move 

	   fname="./rundir-"//Trim(Adjustl(nname))//'/protein-numbers.dat'
       Open(01,File=fname,FORM='Formatted')
       Write(01,*),'Num proteins,mem ',mp%nadd,mp%nmem
       Write(01,*),Trim(Adjustl(nem_profile)),bud_area,num_rings
       Close(01)
       
       surffile="./rundir-"//Trim(Adjustl(nname))//'/surfaceprops.dat'
       energyfile="./rundir-"//Trim(Adjustl(nname))//'/energystats.dat'
       Open(02,File=Trim(Adjustl(surffile)),FORM='Formatted'); Close(02) 
       Open(03,File=Trim(Adjustl(energyfile)),FORM='Formatted'); Close(03) 

       buffsize=0

       mcs_loop: DO mcs=itime,ftime,1                           								    ! Monte Carlo loop 
       Call Monte_Carlo_Steps()

       If(mod(mcs,confint).Eq.0)Then
        Call vtkformat(mcs/confint)                                            
        Call Write_Memb_configuration(mcs/confint)
       Endif

       If(mod(mcs,dataint).Eq.0) Then      
       Call analqtys()
       buffsize=buffsize+1
       surfaceprop(buffsize,:)=(/ Real(mcs/dataint),Real(mp%rg),Real(mp%area),Real(mp%vol) /)
       energyprop(buffsize,:)=(/ Real(mcs/dataint),Real(mp%enel),Real(mp%ennn),Real(mp%ennc) /)
       
       If (buffsize .Eq. nbuffer) Then
		Print*,'dumping data at MCS : ',mcs
		Open(02,File=Trim(Adjustl(surffile)),POSITION='Append')
		Open(03,File=Trim(Adjustl(energyfile)),POSITION='Append')
         Do i1=1,nbuffer
         Write(02,'(4(F21.6,1X))') surfaceprop(i1,:)                                   ! Surface Quantifiers (MCS, rg, area, vol)
         Write(03,'(4(F21.6,1X))') energyprop(i1,:)                                    ! individual Energy Contributions (elastic, nematic, nematic-curv)
         Enddo
        Close(02); Close(03)
		surfaceprop = 0.0; energyprop=0.0;  buffsize=0
       Endif  
       Endif

      Enddo mcs_loop
   
      Call MPI_Barrier(nallgrp,ier)
      Call MPI_FINALIZE(ier)
      END PROGRAM triangulation

!-------------------------------------------------------------------------------------------
!           Initialize the seed for the random number
!-------------------------------------------------------------------------------------------
	Subroutine Initialize_random_number_seed() 
	Use module_datastruct
	Implicit None
    Integer :: clock
    Call SYSTEM_CLOCK(COUNT=clock); seed = -(abs(int(clock*(mynod+1)/8)))                                                        ! get clock time
    Print*,'node: ',mynod,' ran seed: ',seed
	End Subroutine Initialize_random_number_seed

!------------------------------------------------------------------------------------------------------------------------
!                                  Function to generate a random number
!------------------------------------------------------------------------------------------------------------------------  
      FUNCTION ran2(idum)
      IMPLICIT NONE
      Integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real*8  ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1 &
        ,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791 &
        ,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      Integer idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/

      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue  
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END  FUNCTION ran2
!---------------------------------------------------------------------------------------------------------------------------
