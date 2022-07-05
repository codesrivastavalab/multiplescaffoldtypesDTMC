!==========================================================================================================================
!===========================================================================================================================
!                   $moddata                   MODULE CONTAINING THE DATASTRUCTURE
!===========================================================================================================================
!===========================================================================================================================

      MODULE module_datastruct
        IMPLICIT NONE
        INTEGER,PARAMETER::verno=5000,trino=11000,linno=17000
        INTEGER::itime,ftime,noli,nver,ntr,tlink,mcs,seed,iloop          ! mcs-->mcs step num |iloop --> inner loop num          
        INTEGER :: confint,dataint,Nphase,Ann_Nei,Ann_size,mvinter 
        Integer::mynod,nprocs,nallgrp                                    ! variables for parallelizing          
        REAL(KIND=8) :: kappa,fangle,maxtl,pr                            ! itime,ftime --> start and stop time  
        REAL(KIND=8) :: pscur,nscur,pi,necup1,necup2,Ispar,sppar         ! parameters in the problem            
        REAL(KIND=8), DIMENSION(:,:),ALLOCATABLE::nempar,unempar         ! nematic parameter and unscaled value 
		Integer :: cell_lenx,cell_leny,cell_lenz,Ncellx,Ncelly,Ncellz,Ncell			
        Real(Kind=8),Parameter :: zero=0.0,beta=1.0,blen=1.3
        
         TYPE vertex
           REAL(KIND=8),DIMENSION(3,1)::splo,spgl,vnor,vcoord,t1,t2      ! Orientation of the nematic in local and global sys 
           INTEGER,DIMENSION(10):: vneipt,vneitr                         ! List of all neigh points,triangles,links           
           REAL(KIND=8),DIMENSION(3,3)::L2G
           INTEGER ::nonei,phase,clno,neigh                              ! Number of Neighbours a vertex has                     
           REAL(KIND=8)::mcur,cur1,cur2,totarea,spen,op,WN               ! WN is the winding number of a loop around the vertex  
           REAL(KIND=8):: kap,kpar,kper,cpar,cper,Nemal,Isal             ! The first four are surface  directional parameters    
           REAL(KIND=8):: rkap,rkpar,rkper,rcpar,rcper,rNemal,rIsal      ! Reference values with respect to which TI will be performed
		   REAL(KIND=8):: ukpar,ukper                                    ! unscaled kappa parallel and kappa perp at each vertex 
		   INTEGER :: linkcell,next_vert,prev_vert	   	      
         END TYPE vertex                                                 ! The last two are relative nematic and spin alignments 

         TYPE triangle                                                   ! Triangle details                                 
           REAL(KIND=8)::ar,vol                                          ! Area of each face and their corresponding volume 
           INTEGER,DIMENSION(3) ::li,vert                                ! Links that make up the triangle                  
           REAL(KIND=8),DIMENSION(3,1)::fnor
         END TYPE triangle 
         
         TYPE link
           INTEGER ::tr                                                  !Triangle associated 
           INTEGER,DIMENSION(2) ::sep                                    !Begin and end vertex  
         END TYPE link
    
      TYPE membrane_prop
           REAL(KIND=8) :: enel,ento,ennn,ennc,enIsing,lam_ennn,lam_ennc   ! The different energy terms involved 
	   REAL(KIND=8) :: lambda_ener
           REAL(KIND=8) :: vol,area,rg,op,pdef,ndef                      ! Volume, area, Radius of Gyration    
           Integer      :: phdef,nhdef,p1def,n1def
           REAL(KIND=8) :: imcur,imcursq,lambda
           INTEGER :: nadd,nmem,noclus,Inlen
  	   Type(vertex),Allocatable :: vertex(:)
           Type(triangle),Allocatable :: triangle(:)
           Type(link),Allocatable :: link(:)
	   Character(100) :: Lambda_variable

      END TYPE membrane_prop 

      Type Linklist
         Real(Kind=8) :: start_coord(3,1),end_coord(3,1)                 ! Spatial limit of the given cell
	 Integer :: first_vert,last_vert,num_vert,buffer
	 Integer :: neigh_cell(27)
      End Type Linklist

       
        TYPE(vertex),Allocatable::ver(:)                                               ! Setting the size for type vertex
        TYPE(triangle),Allocatable::tri(:)                                             ! Setting the size for type triangle
        TYPE(link),Allocatable::lin(:)                                                 ! Setting the size for type link
        TYPE(membrane_prop)::mp
	Type(Linklist),Allocatable :: cell(:)                                          ! Size of linklist cells is set to allocatable	
         
      END MODULE module_datastruct

