!===========================================================================================================================
!===========================================================================================================================
!      $$ modmcsmoves                          MODULE TO FLIP THE BONDS AND DISPLACE THE CHOSEN VERTEX POINTS
!===========================================================================================================================
!===========================================================================================================================

      MODULE module_mcsmoves
      IMPLICIT NONE
      INTEGER ::norchk,amov,aflip,tmov,tflip
      REAL(KIND=8)::flen,mven,dmven,prob,rnum,Enecur                     ! flen,mven,nemven are temp var for energy cal 
      REAL(KIND=8)::dfang,Ising,LebLas,SIsing  
      contains

!---------------------------------------------------------------------------------------------------------------------------
!                                        Randomize the Monte Carlo Steps
!--------------------------------------------------------------------------------------------------------------------------- 
      SUBROUTINE Monte_Carlo_Steps()
      USE module_datastruct ;USE module_curvcalc ;USE module_rerun
      IMPLICIT NONE                                                      ! Perform all MC moves in all possible sequences 
       
      REAL(KIND=8)::ran2                                               
      Integer :: rnumber
      
       rnumber=NINT(24.0*ran2(Seed))+1
       If(rnumber.Gt.24) rnumber=24                                            ! Number of combinations = (number of MC moves)!  

       Choose_MC_Seq:Select Case(rnumber)                        
       Case(1)
       DO iloop=1,tlink,1                                                ! Placing loops inside case statement speed up the code
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       ENDDO

       Case(2)
       DO iloop=1,tlink,1                                                ! Inner loop       
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       ENDDO

       Case(3)
       DO iloop=1,tlink,1                                                ! Inner loop       
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       ENDDO

       Case(4)
       DO iloop=1,tlink,1                                                ! Inner loop       
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       ENDDO

       Case(5)
       DO iloop=1,tlink,1                                                ! Inner loop       
       Call Flip_Link()
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       ENDDO

       Case(6)
       DO iloop=1,tlink,1                                                ! Inner loop       
       Call Flip_Link()
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       ENDDO

       Case(7)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       Call Flip_Link()
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       ENDDO

       Case(8)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       ENDDO

       Case(9)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF   
       If(mod(iloop,mvinter).Eq.0)Then;Call Nematic_Flip() ; ENDIF 
       Call Flip_Link()
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       ENDDO
      
       Case(10)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       Call Flip_Link()
       ENDDO

       Case(11)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF   
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then;Call Nematic_Flip() ; ENDIF 
       ENDDO
      
       Case(12)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       Call Flip_Link()
       ENDDO

       Case(13)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       ENDDO

       Case(14)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       ENDDO

       Case(15)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF   
       If(mod(iloop,mvinter).Eq.0)Then;Call Nematic_Flip() ; ENDIF 
       ENDDO
      
       Case(16)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       ENDDO

       Case(17)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       If(mod(iloop,mvinter).Eq.0)Then;Call Nematic_Flip() ; ENDIF 
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF   
       Call Flip_Link()
       ENDDO
      
       Case(18)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       ENDDO

       Case(19)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF     
       ENDDO

       Case(20)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       Call Flip_Link()
       ENDDO

       Case(21)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       Call Flip_Link()
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF      
       ENDDO

       Case(22)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       ENDDO

       Case(23)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF      
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif
       Call Flip_Link()
       ENDDO

       Case(24)
       DO iloop=1,tlink,1                                                ! Inner loop       
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       Call Flip_Link()
       If(mp%nadd.Lt. nver)Then; Call Exchange_Nematic() ; Endif       
       ENDDO
       End Select Choose_MC_Seq
       End Subroutine Monte_Carlo_Steps

!---------------------------------------------------------------------------------------------------------------------------
!                                        Randomize the Monte Carlo Steps
!--------------------------------------------------------------------------------------------------------------------------- 
      SUBROUTINE Monte_Carlo_Steps_fixednematic()
      USE module_datastruct ;USE module_curvcalc ;USE module_rerun
      IMPLICIT NONE                                                                            ! Perform all MC moves in all possible sequences 
       
      REAL(KIND=8)::ran2                                               
      Integer :: rnumber
      
       rnumber=NINT(24.0*ran2(Seed))+1
       If(rnumber.Gt.24) rnumber=24                                                                  ! Number of combinations = (number of MC moves)!  

       Choose_MC_Seq_fixednematic:Select Case(rnumber)                        
       Case(1)
       DO iloop=1,tlink,1                                                                      ! Placing loops inside case statement speed up the code
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       ENDDO

       Case(2)
       DO iloop=1,tlink,1                                                
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       ENDDO

       Case(3)
       DO iloop=1,tlink,1                                                       
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       ENDDO

       Case(4)
       DO iloop=1,tlink,1                                                     
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       ENDDO

       Case(5)
       DO iloop=1,tlink,1                                                      
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       ENDDO

       Case(6)
       DO iloop=1,tlink,1                                                       
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       ENDDO

       Case(7)
       DO iloop=1,tlink,1                                                       
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       ENDDO

       Case(8)
       DO iloop=1,tlink,1                                                      
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       ENDDO

       Case(9)
       DO iloop=1,tlink,1                                                       
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF   
       If(mod(iloop,mvinter).Eq.0)Then;Call Nematic_Flip() ; ENDIF 
       Call Flip_Link()
       ENDDO
      
       Case(10)
       DO iloop=1,tlink,1                                                       
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       Call Flip_Link()
       ENDDO

       Case(11)
       DO iloop=1,tlink,1                                                      
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF   
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then;Call Nematic_Flip() ; ENDIF 
       ENDDO
      
       Case(12)
       DO iloop=1,tlink,1                                                       
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       Call Flip_Link()
       ENDDO

       Case(13)
       DO iloop=1,tlink,1                                                 
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       ENDDO

       Case(14)
       DO iloop=1,tlink,1                                                    
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       ENDDO

       Case(15)
       DO iloop=1,tlink,1                                               
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF   
       If(mod(iloop,mvinter).Eq.0)Then;Call Nematic_Flip() ; ENDIF 
       ENDDO
      
       Case(16)
       DO iloop=1,tlink,1                                                
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       ENDDO

       Case(17)
       DO iloop=1,tlink,1                                                 
       If(mod(iloop,mvinter).Eq.0)Then;Call Nematic_Flip() ; ENDIF 
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF   
       Call Flip_Link()
       ENDDO
      
       Case(18)
       DO iloop=1,tlink,1                                                     
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       ENDDO

       Case(19)
       DO iloop=1,tlink,1                                                    
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF     
       ENDDO

       Case(20)
       DO iloop=1,tlink,1                                                     
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       Call Flip_Link()
       ENDDO

       Case(21)
       DO iloop=1,tlink,1                                            
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF      
       ENDDO

       Case(22)
       DO iloop=1,tlink,1                                                
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       Call Flip_Link()
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       ENDDO

       Case(23)
       DO iloop=1,tlink,1                                                
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF      
       Call Flip_Link()
       ENDDO

       Case(24)
       DO iloop=1,tlink,1                                              
       If(mod(iloop,mvinter).Eq.0)Then; Call Nematic_Flip();ENDIF
       If(mod(iloop,mvinter).Eq.0)Then; Call Vertex_Move() ; ENDIF
       Call Flip_Link()
       ENDDO
       End Select Choose_MC_Seq_fixednematic
       End Subroutine Monte_Carlo_Steps_fixednematic

!---------------------------------------------------------------------------------------------------------------------------
!                                                       SUBROUTINE TO FLIP A LINK
!--------------------------------------------------------------------------------------------------------------------------- 
      SUBROUTINE Flip_Link()
      USE module_datastruct ;USE module_curvcalc ;USE module_rerun ; Use module_energy_calculations
      IMPLICIT NONE 
       
      REAL(KIND=8)::ran2,inen,fien       
      INTEGER ::vt1,vt2,i,j,ep1,ep2,fv1,fv2       							                    ! vt -->vert of triangle; norchk -->nor check  
      INTEGER ::fvp1,fvm1,fvp2,fvm2,tn1,tn2,rnumber,mchk       
      INTEGER :: lp1,lp2,ch1,trn,llist(6)
      INTEGER,DIMENSION(4)::cver                  	    									    ! Temporary variables for case statement        
      INTEGER,DIMENSION(3)::t1,t2,tmp,tmp1                                     			        ! ep -->endpoint of link                        

      rnumber=nint((1-2*ran2(seed))*tlink)                                                      ! Choose a Random link 
      IF (rnumber.Eq.0 .Or. abs(rnumber).Gt.tlink) Return

       t1=tri(lin(rnumber)%tr)%vert ; t2=tri(lin(-rnumber)%tr)%vert 
       tn1=lin(rnumber)%tr ; tn2=lin(-rnumber)%tr                                               ! tn1 and tn2 are the names of the triangles    
       lp1=lin(rnumber)%sep(1) ; lp2=lin(rnumber)%sep(2)
       ep1=0;ep2=0 ;inen=0.0 ;fien=0.0; fv1=0; fv2=0

       llist=(/tri(tn1)%li,tri(tn2)%li/)

        DO i=1,3  
        IF(t1(i).EQ.lp2)THEN ; ep1=i ; ENDIF                                                     ! Position of vertex where the chosen bond ends   
        IF(t2(i).EQ.lp1)THEN ; ep2=i ; ENDIF                                                     ! The position of vert where conj chosen bond ends

        IF((t1(i).NE.lp1).AND.(t1(i).NE.lp2))then
        vt1=t1(i);fv1=i                                                                          ! Position of the free vertex in triang1
        ENDIF                                                                                    ! fv1 and fv2 are the position of the free vert

        IF((t2(i).NE.lp2).AND.(t2(i).NE.lp1))then
        vt2=t2(i);fv2=i                                                                          ! Position of the free vertex in triang2
        ENDIF 
        ENDDO

        cver=(/vt1,vt2,lp1,lp2/)                                                                  ! All relevant vertice num in one array
        
        DO i=1,ver(vt1)%nonei,1                                         
        IF(ver(vt1)%vneipt(i).Eq.vt2) RETURN                                                       ! Do not proceed if already connected
        ENDDO

		If(Neigh_blcheck(vt2,vt1) .Ne.0) Return                                                    ! Check for the new bond length        
         	
        IF((ver(vt1)%nonei.GE.9).OR.(ver(vt2)%nonei.GE.9))  RETURN                                 ! Maximum limit on the number of neighbours
        IF((ver(lp1)%nonei.LE.3).AND.(ver(lp2)%nonei.LE.3)) RETURN                                 ! Minimum limit on the number of neighbours

		mp%triangle(tn1)=tri(tn1) ; mp%triangle(tn2)=tri(tn2)
        Do i=1,6,1
         mp%link(llist(i))=lin(llist(i))
        Enddo
        Do i=1,4,1
         mp%vertex(cver(i))=ver(cver(i))                                                            ! Original state of involved vertices
        Enddo

		Call energy_linkflip(cver,tn1,tn2,inen)

        t1(ep1)=vt2 ; t2(ep2)=vt1                                                                  ! Free vertices are connected to new trian

          j=ver(cver(3))%nonei                                                                     ! lp2 is removed from list of lp1 
          DO i=1,j-1,1
          IF(ver(cver(3))%vneipt(i).EQ.cver(4))THEN
          ver(cver(3))%vneipt(i:j-1)=ver(cver(3))%vneipt(i+1:j)
          ENDIF
          ENDDO
          ver(cver(3))%vneipt(j:10)=0                                                               ! All sites beyond the neighbour size=0 
          ver(cver(3))%nonei=j-1

          j=ver(cver(4))%nonei                                                                      ! lp1 is removed from the list of lp2 
          DO i=1,j-1,1
          IF(ver(cver(4))%vneipt(i).EQ.cver(3))THEN
          ver(cver(4))%vneipt(i:j-1)=ver(cver(4))%vneipt(i+1:j)
          ENDIF
          ENDDO
          ver(cver(4))%vneipt(j:10)=0                                                               ! All sites above neigh size=0 
          ver(cver(4))%nonei=j-1
  
       
         tmp1=0  
         DO i=1,3,1                                                                                                                 ! To find the the links whose triangles will change
         IF(lin(tri(tn1)%li(i))%sep(2).EQ.vt1) tmp1(1)=tri(tn1)%li(i)
         IF(lin(tri(tn2)%li(i))%sep(2).EQ.vt2) tmp1(2)=tri(tn2)%li(i)
         ENDDO


        lin(rnumber)%sep=(/vt2,vt1/)                                                                                                ! Update the start&end of the chosen rnumberom link 
        lin(-rnumber)%sep=(/vt1,vt2/)                                
        tri(tn1)%vert=t1                                                                                                            ! Update the Vertex of the triangle 
        tri(tn2)%vert=t2  

        fvp1=fv1+1 ; IF(fv1.EQ.3)fvp1=1                                                                                             ! Circular boundary conditions
        fvm1=fv1-1 ; IF(fv1.EQ.1)fvm1=3
        fvp2=fv2+1 ; IF(fv2.EQ.3)fvp2=1
        fvm2=fv2-1 ; IF(fv2.EQ.1)fvm2=3

        tmp=0 
        tmp(1)=tri(tn1)%li(fvm1)                                                                                                    ! tmp is a temp array that stores links
        tmp(2)=tri(tn2)%li(fvm2)

        tri(tn1)%li(fvm1)=tri(tn1)%li(fvp1)                                                                                         ! Updating the triangle --link
        tri(tn2)%li(fvm2)=tri(tn2)%li(fvp2)
        tri(tn1)%li(fvp1)=tmp(2)
        tri(tn2)%li(fvp2)=tmp(1)

        lin(tmp1(1))%tr=tn2                                                                                                         ! Updating the link - triangle
        lin(tmp1(2))%tr=tn1

        ver(cver(1))%nonei=ver(cver(1))%nonei+1                                                                                     ! Increment the neigh of ver1 by  1
        ver(cver(2))%nonei=ver(cver(2))%nonei+1                                                                                     ! Increment the neigh of ver2 by  1

        ch1=0 ; i=1
        f1_ver:Do WHILE((i.LE.mp%vertex(cver(1))%nonei).And.(ch1.EQ.0)) 
        If(mp%vertex(cver(1))%vneipt(i).EQ.cver(3))then
        Do j=i+1,mp%vertex(cver(1))%nonei
        ver(cver(1))%vneipt(j+1)=mp%vertex(cver(1))%vneipt(j)                                                                       ! The change in neighbouring  vertex is put in 
        Enddo 
        ver(cver(1))%vneipt(i+1)=cver(2) ; ch1=1
        Endif
        i=i+1
        Enddo f1_ver
	
        ch1=0 ; i=1
        f2_ver:Do WHILE((i.LE.mp%vertex(cver(2))%nonei).And.(ch1.EQ.0))                                                             ! The vertex order is changed
        If(mp%vertex(cver(2))%vneipt(i).EQ.cver(4))then
         Do j=i+1,mp%vertex(cver(2))%nonei                                                                                          ! Chosen link goes from cver(3) to cver(4)        
         ver(cver(2))%vneipt(j+1)=mp%vertex(cver(2))%vneipt(j)                                                                      !      cver(4)               4                    
         Enddo                                                                                                                      !         *                  *  After flip        
         ver(cver(2))%vneipt(i+1)=cver(1) ;ch1=1                                                                                    !       / | \              /   \ chosen link goes 
         Endif                                                                                                                      !      /  |  \            /     \ from 2 to 1     
         i=i+1                                                                                                                      ! vt1 /   |   \  vt2     /  tn2  \                
        Enddo f2_ver                                                                                                                !(or)*tn1 |tn2 *(or)  1 *---------* 2             
                                                                                                                                    !cver \   |   / cver(2)  \  tn1  /                
        ch1=0 ; i=1                                                                                                                 ! (1)  \  |  /            \     /                 
        f1_tr:Do WHILE((i.LE.mp%vertex(cver(1))%nonei).And.(ch1.EQ.0))                                                              !       \ | /              \   /                  
        If(mp%vertex(cver(1))%vneitr(i).EQ.tn1)Then                                                                                 !         *                  *                    
         Do j=i+1,mp%vertex(cver(1))%nonei                                                                                          !       cver(3)              3                    
         ver(cver(1))%vneitr(j+1)=mp%vertex(cver(1))%vneitr(j)
         Enddo 
         ver(cver(1))%vneitr(i+1)=tn2 ;ch1=1                                                                                        ! The first triangle free vertex is linked to tri 2
         Endif
         i=i+1
        Enddo f1_tr

        ch1=0 ; i=1
        f2_tr: Do WHILE((i.LE.mp%vertex(cver(2))%nonei).And.(ch1.EQ.0))
         If(mp%vertex(cver(2))%vneitr(i).EQ.tn2)Then
         Do j=i+1,mp%vertex(cver(2))%nonei
         ver(cver(2))%vneitr(j+1)=mp%vertex(cver(2))%vneitr(j)                                                                      ! The first triangle free vertex is linked to tri 2
         Enddo 
         ver(cver(2))%vneitr(i+1)=tn1 ;ch1=1
         Endif
         i=i+1
        Enddo f2_tr

        ch1=0 ; i=1
        f3_tr:Do WHILE(i.LE.mp%vertex(cver(3))%nonei .And.(ch1.EQ.0))                                                               ! Updating the neig trian list of the vertex where 
         If(ver(cver(3))%vneitr(i).EQ.tn2)Then                                                                                      ! the chosen bond started earlier(either t1 or t2) 
          Do j=i,ver(cver(3))%nonei                                                                                                 ! removed
          ver(cver(3))%vneitr(j)=ver(cver(3))%vneitr(j+1)
          Enddo                
          If(i.EQ.mp%vertex(cver(3))%nonei)Then                                                                                     ! If tn2 is last in the list it is rearranged for
          ver(cver(3))%vneitr(mp%vertex(cver(3))%nonei)=ver(cver(3))%vneitr(1)                                                      ! matching the link list
          ver(cver(3))%vneitr(1:ver(cver(3))%nonei)=ver(cver(3))%vneitr(2:mp%vertex(cver(3))%nonei)
          Endif
          ch1=1
         Endif 
         i=i+1
        Enddo f3_tr
        ver(cver(3))%vneitr(mp%vertex(cver(3))%nonei:10)=0 
     
        ch1=0 ; i=1
        f4_tr:Do WHILE((i.LE.mp%vertex(cver(4))%nonei).And.(ch1.EQ.0))                                                              ! Updating the neig trian list of the vertex 
         If(ver(cver(4))%vneitr(i).EQ.tn1) Then
         Do j=i,ver(cver(4))%nonei                                 
         ver(cver(4))%vneitr(j)=ver(cver(4))%vneitr(j+1)
         Enddo                
         If(i.EQ.mp%vertex(cver(4))%nonei)Then
         ver(cver(4))%vneitr(mp%vertex(cver(4))%nonei)=ver(cver(4))%vneitr(1)
         ver(cver(4))%vneitr(1:ver(cver(4))%nonei)=ver(cver(4))%vneitr(2:mp%vertex(cver(4))%nonei)
         Endif 
         ch1=1
         Endif
         i=i+1
        Enddo f4_tr
        ver(cver(4))%vneitr(mp%vertex(cver(4))%nonei:10)=0 

        Call onlyarea(tn1)                                                                                                          ! Calculate the area of the two new triangles 
        Call onlyarea(tn2)  

        DO i=1,4,1
        ver(cver(i))%totarea=0
         DO j=1,ver(cver(i))%nonei
         trn=ver(cver(i))%vneitr(j)
         ver(cver(i))%totarea=ver(cver(i))%totarea+tri(trn)%ar/3.0
         ENDDO
        ENDDO

        
        flen=0 ; mchk=0
                                                     
        Call faceangchk(tn1)                                        							            ! Check  angle between tri 1 and neigh faces       
        fang_tr1:IF(norchk .EQ.0) THEN                 								                    ! Proceed further only if Yes (Max =150 \degrees)  
        Call faceangchk(tn2)                             									    ! Check for angle between tri2 and its neigh faces 
        fang_tr2 : IF(norchk .EQ.0) THEN                        

        DO i=1,4,1
        Call normalcalc(cver(i))              										            ! Update the curvature at the chosen vertices 
        ENDDO

	Call energy_linkflip(cver,tn1,tn2,fien)
                      
        IF(((fien-inen).GT.0).And.(ran2(seed).Gt.exp(-beta*(fien-inen))))  mchk=1                                                   ! pass the control of the two faces are perpn 

        ELSE
        mchk=1 
        ENDIF fang_tr2
        ELSE
        mchk=1 
        ENDIF fang_tr1

        flip_failed:IF(mchk==1)then
        tri(tn1)=mp%triangle(tn1) ; tri(tn2)=mp%triangle(tn2)
        ForAll(i=1:6) lin(llist(i))=mp%link(llist(i))
        ForAll(i=1:4) ver(cver(i))=mp%vertex(cver(i))
        ENDIF flip_failed

      END SUBROUTINE Flip_Link 
!---------------------------------------------------------------------------------------------------------------------------
!                  $vermov           SUBROUTINE TO MOVE THE VERTEX
!---------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE Vertex_Move()
        USE module_datastruct; USE module_curvcalc ; Use module_rerun
       	USE module_energy_calculations ; Use module_linklist_calc
        IMPLICIT NONE                                                                                                  
        INTEGER::vert,i,j,k,trno,new_cell
        REAL(KIND=8)::ran2,inen,fien                                                                                           	    ! ran2 fn, init and fin ener
        REAL(KIND=8),DIMENSION(3,1)::dr                                                                                             ! Values of the Nematic Vector            

        vert=NINT(ran2(Seed)*nver)+1
		If(vert.Gt.nver) Return	

        inen=0; fien=0
        mp%vertex(vert)=ver(vert)
	
        DO i=1,3
        dr(i,1)=(1.0-2*ran2(SEED))*0.05                                                                                             ! A Small displacement vector
        ENDDO

!        IF(SUM(dr**2).GT.0.0025) RETURN                                                                                             ! Max allowed displacement is sqrt(0.0025)
       	ver(vert)%vcoord=ver(vert)%vcoord+dr                                                                                       !New displaced position 
        
          new_cell=get_cellnumber(vert)
          If(cell(new_cell)%buffer .Eq. 1) Then
          ver(vert)%vcoord=mp%vertex(vert)%vcoord
          Call Initialize_linkcells()
          Return
          Endif

        IF(Neigh_blcheck(vert) .NE.0 ) then                       
        ver(vert)%vcoord=mp%vertex(vert)%vcoord
        RETURN 
        ENDIF                                                                                                                       ! with all neighbouring vertices.

		If(Minimum_bondlength_check(vert,new_cell).Eq.1)Then
		ver(vert)%vcoord=mp%vertex(vert)%vcoord
		Return
		Endif

		nei_ar: DO i=1,ver(vert)%nonei 
		j=ver(vert)%vneipt(i) ; k=ver(vert)%vneitr(i)                                                                               ! Neigh points and triangles
		mp%vertex(j)=ver(j) ; mp%triangle(k)=tri(k)
	ENDDO nei_ar

	Call energy_vertexmove(vert,inen)

    	upd_ntr_area:  DO i=1,ver(vert)%nonei
		Call areacal(ver(vert)%vneitr(i))                                                                                           ! Calculate the new area of each face
        ENDDO upd_ntr_area 

        norchk=0 

        fang_chk:DO trno=1,ver(vert)%nonei                                                                                          ! Call the face angle check for all triangles- 
        Call faceangchk(ver(vert)%vneitr(trno))                                                                                     ! -around the vertex moved                     
        IF(norchk.NE.0)THEN
        ver(vert)%vcoord=mp%vertex(vert)%vcoord                          
        ver(vert)%totarea=mp%vertex(vert)%totarea
        res_ntr: DO i=1,ver(vert)%nonei,1                                                                                           ! Update all neighbours       
        j=ver(vert)%vneipt(i) ; k=ver(vert)%vneitr(i)                                                                               ! Neigh triangles|}variables 
        tri(k)%ar=mp%triangle(k)%ar        
        ver(j)%totarea=mp%vertex(j)%totarea 
        tri(k)%fnor=mp%triangle(k)%fnor
        ENDDO res_ntr
        RETURN
        ENDIF
        ENDDO  fang_chk

        Call normalcalc(vert)                                                                                                       !Calculate the curvature of the moved vertex

        upd_nptcur: DO i=1,ver(vert)%nonei                                                                                          !Calculate the curvature of all neighbours
        Call normalcalc(ver(vert)%vneipt(i))
        ENDDO upd_nptcur   
         
		Call energy_vertexmove(vert,fien)

        mvmetropolis: IF(((fien-inen).GT.0.0).And.(ran2(seed).Gt.exp(-beta*(fien-inen)))) THEN
        ver(vert)=mp%vertex(vert)
        rest_neig: DO i=1,ver(vert)%nonei,1                                 
        j=ver(vert)%vneipt(i) ; k=ver(vert)%vneitr(i)              
        ver(j)=mp%vertex(j) ; tri(k)=mp%triangle(k)
        ENDDO rest_neig
		Else
		!new_cell=get_cellnumber(vert)
        If(new_cell.Ne. ver(vert)%linkcell) Call update_vertex_cell(vert,new_cell)
        ENDIF mvmetropolis

	END SUBROUTINE Vertex_Move

!---------------------------------------------------------------------------------------------------------------------------
!                               SUBROUTINE TO MAKE AN EXCHANGE OPERATION BETWEEN NEIGHBOURING SPINS
!---------------------------------------------------------------------------------------------------------------------------
           SUBROUTINE Exchange_Nematic()
           USE module_datastruct ; USE module_curvcalc ; Use module_energy_calculations
           IMPLICIT NONE
           INTEGER :: v1,v2,Lnk,phase1,phase2
           REAL(KIND=8) :: splo1(3,1),splo2(3,1),spgl1(3,1),spgl2(3,1) 
           REAL(KIND=8) :: ran2,inen,fien,prob,necuen,cpar(2), cper(2)
           REAL(KIND=8) :: en1,en2,parkap(2),perpkap(2)
           REAL(KIND=8) :: temp(4),spenergy,necuenergy


           Lnk=nint((1-2*ran2(seed))*tlink)                               !Kawasaki exchange 
           If (Lnk.Eq.0 .Or. abs(Lnk).Gt.tlink) Return
           v1=lin(Lnk)%sep(1)   ; v2=lin(Lnk)%sep(2)
           phase1=ver(v1)%phase ; phase2=ver(v2)%phase                    ! type of nematics at each vertex
           
           if (phase1 .Eq. phase2) Return                                 ! do not exchage same type of nematics (very useful for low densities)

           splo1=ver(v1)%splo   ; splo2=ver(v2)%splo
           spgl1=ver(v1)%spgl   ; spgl2=ver(v2)%spgl

           parkap=(/ver(v1)%kpar,ver(v2)%kpar/)                                                                   ! \kappa_{\parallel} at v1 and v2
           perpkap=(/ver(v1)%kper,ver(v2)%kper/)                                                                  ! \kappa_{\perp} at v1 and v2
           cpar=(/ver(v1)%cpar,ver(v2)%cpar/) 
           cper=(/ver(v1)%cper,ver(v2)%cper/) 
               
           en1=0.0 ; en2=0.0    ; inen=0.0 
           Call Nematic_Interaction_Potential(spenergy,v1)
           necuenergy=(parkap(1)*((ver(v1)%splo(1,1)**2*ver(v1)%cur1+ver(v1)%splo(2,1)**2*ver(v1)%cur2)-cpar(1))**2+&
           &perpkap(1)*((ver(v1)%splo(2,1)**2*ver(v1)%cur1+ver(v1)%splo(1,1)**2*ver(v1)%cur2)-cper(1))**2)*ver(v1)%totarea
           en1=necuenergy+spenergy                                                                                                        !  Energy (necu + nene) at vertex 1 

           Call Nematic_Interaction_Potential(spenergy,v2)
           necuenergy=(parkap(2)*((ver(v2)%splo(1,1)**2*ver(v2)%cur1+ver(v2)%splo(2,1)**2*ver(v2)%cur2)-cpar(2))**2+&
           &perpkap(2)*((ver(v2)%splo(2,1)**2*ver(v2)%cur1+ver(v2)%splo(1,1)**2*ver(v2)%cur2)-cper(2))**2)*ver(v2)%totarea  	    ! Nematic -curvature part for the neighbours         
     
           en2=necuenergy+spenergy                                                                                                      ! Energy at vertex 2 
           inen=en1+en2                                                                                                             ! Initial energy of the system 
            
           ver(v1)%splo=splo2;ver(v1)%phase=phase2                                                                                  ! Exchange the nematic order( Flag + field).
           ver(v2)%splo=splo1;ver(v2)%phase=phase1                                                                                  ! orientation of nematic in the local frame 
           ver(v1)%spgl=MATMUL(ver(v1)%L2G,ver(v1)%splo)                                                                            ! The orientation in the global cooord system changes
           ver(v2)%spgl=MATMUL(ver(v2)%L2G,ver(v2)%splo)

           en1=0.0 ; en2=0.0 ; fien=0.0 
           Call Nematic_Interaction_Potential(spenergy,v1)
           necuenergy=(parkap(2)*((ver(v1)%splo(1,1)**2*ver(v1)%cur1+ver(v1)%splo(2,1)**2*ver(v1)%cur2)-cpar(2))**2+&
           &perpkap(2)*((ver(v1)%splo(2,1)**2*ver(v1)%cur1+ver(v1)%splo(1,1)**2*ver(v1)%cur2)-cper(2))**2)*ver(v1)%totarea          ! Nematic -curvature part for the neighbours 
           en1=necuenergy+spenergy                                                                                                      !  Energy (necu + nene) at vertex 1 

           necuen=0.0
           Call Nematic_Interaction_Potential(spenergy,v2)
           necuenergy=(parkap(1)*((ver(v2)%splo(1,1)**2*ver(v2)%cur1+ver(v2)%splo(2,1)**2*ver(v2)%cur2)-cpar(1))**2+&
           &perpkap(1)*((ver(v2)%splo(2,1)**2*ver(v2)%cur1+ver(v2)%splo(1,1)**2*ver(v2)%cur2)-cper(1))**2)*ver(v2)%totarea
           en2=necuenergy+spenergy                                                                                                        ! Energy at vertex 2 
           fien=en1+en2                                                                                                               ! Final energy after the exchange of spins 
 
        sfmetro1:IF((fien-inen).GT.0) THEN                           
          prob=(1-tanh(beta*(fien-inen)*0.5))                                                                                     ! Probability is calculated in the Kawasaki scheme   
          sfmetro2: IF(prob.LT.ran2(SEED)) THEN 
          ver(v1)%phase=phase1 ; ver(v2)%phase=phase2                   
          ver(v1)%splo=splo1   ; ver(v1)%spgl=spgl1                                                                                  ! Restore the original spin in local and global frame 
          ver(v2)%splo=splo2   ; ver(v2)%spgl=spgl2
          Else
          temp=(/ver(v2)%kpar,ver(v2)%kper,ver(v2)%cpar,ver(v2)%cper/)                                                               ! Exchange the parameters along with the spins if accepted 
          ver(v2)%kpar=ver(v1)%kpar ; ver(v2)%kper=ver(v1)%kper 
          ver(v2)%cpar=ver(v1)%cpar ; ver(v2)%cper=ver(v1)%cper
          ver(v1)%kpar=temp(1) ; ver(v1)%kper=temp(2) 
          ver(v1)%cpar=temp(3) ; ver(v1)%cper=temp(4) 
          ENDIF sfmetro2
         Else
         temp=(/ver(v2)%kpar,ver(v2)%kper,ver(v2)%cpar,ver(v2)%cper/)
         ver(v2)%kpar=ver(v1)%kpar ; ver(v2)%kper=ver(v1)%kper 
         ver(v2)%cpar=ver(v1)%cpar ; ver(v2)%cper=ver(v1)%cper
         ver(v1)%kpar=temp(1) ; ver(v1)%kper=temp(2) 
         ver(v1)%cpar=temp(3) ; ver(v1)%cper=temp(4)
         ENDIF sfmetro1 
         END SUBROUTINE Exchange_Nematic
!---------------------------------------------------------------------------------------------------------------------------
!                               SUBROUTINE TO FLIP THE CHOSEN SPIN
!---------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE Nematic_Flip()
        USE module_datastruct ; USE module_curvcalc ; Use module_energy_calculations
        IMPLICIT NONE
        INTEGER :: v
        REAL(KIND=8) ::ran2,inen,fien,prob,spenergy,necuenergy
        REAL(KIND=8),DIMENSION(3,1) :: ospl,ospg,rannem              
        
       v=nint(ran2(seed)*nver)+1                                                                                                    ! spin flip  move  
       IF (v.Gt.nver .Or. ver(v)%phase.Eq.2) Return
       
       ospl=ver(v)%splo                                                                                                             ! Store the old global spin  
       ospg=ver(v)%spgl                                                                                                             ! Store the old local spin   
       Call Nematic_Interaction_Potential(spenergy,v)
       Call Nematic_Curvature_Potential(necuenergy,v)
       inen=spenergy+necuenergy

       rannem(1,1)=0.5*(1-2*ran2(SEED))                                                                                             ! Test nematic vector in local frame   
       rannem(2,1)=0.5*(1-2*ran2(SEED))
       rannem(3,1)=0.0! 0.5*(1-2*ran2(SEED))
       rannem=ospl+rannem                                                                                                           ! Rotation of the nematic is carried out in the local frame since  
       rannem=rannem/SQRT(sum(rannem**2))                                                                                           ! it eases the calculation of energy    
       ver(v)%splo=rannem
       ver(v)%spgl=MATMUL(ver(v)%L2G,rannem) 

        Call Nematic_Interaction_Potential(spenergy,v)  
        Call Nematic_Curvature_Potential(necuenergy,v)
        fien=spenergy+necuenergy

         sfmetro1:IF((fien-inen).GT.0) THEN                                                                                         ! Metropolis scheme
         prob=exp(-beta*(fien-inen))
         sfmetro2: IF(prob.LT.ran2(SEED)) THEN
         ver(v)%spgl=ospg
         ver(v)%splo=ospl
         ENDIF sfmetro2
         ENDIF sfmetro1
       END SUBROUTINE Nematic_Flip    
 
!------------------------------------------------------------------------------------------------------------------------
!                                SUBROUTINE TO CHECK FOR THE ANGLE BETWEEN THE GIVEN FACES
!------------------------------------------------------------------------------------------------------------------------
          SUBROUTINE faceangchk(trno)
          USE module_datastruct
          IMPLICIT NONE
          INTEGER :: trno,j,k,c1                                ! trno --> Triangle to be chked for
          REAL(KIND=8) :: npro(1,1)                                      ! Dot product of the Normals       
          INTEGER,DIMENSION(3) :: ll
          norchk=0 ; ll=0 ; c1=0
              
		  ll= tri(trno)%li
          DO k=1,3
           j=lin(-ll(k))%tr 
           npro=MATMUL(TRANSPOSE(tri(trno)%fnor),tri(j)%fnor)  
           IF(npro(1,1).LT.fangle) THEN
            norchk=1                                                     ! norchk --> 0 for the conf to be accepted
            RETURN
           ENDIF
          ENDDO
          RETURN
          END SUBROUTINE faceangchk

!------------------------------------------------------------------------------------------------------------------
!                 SUBROUTINE TO CALCULATE THE AREA OF THE GIVEN TRIANGLE
!------------------------------------------------------------------------------------------------------------------
          SUBROUTINE onlyarea(tr)
          USE module_datastruct  
          IMPLICIT NONE
          INTEGER::i,j,k,tr
          REAL(KIND=8) ::ax,ay,az,area
          REAL(KIND=8),DIMENSION(3,1)::r1,r2,r3,r21,r31

          i=tri(tr)%vert(1) ; j=tri(tr)%vert(2) ; k=tri(tr)%vert(3)      ! The vertices that make up the triangle tr 
          r1=ver(i)%vcoord ; r2=ver(j)%vcoord ; r3=ver(k)%vcoord         ! Their corresponding coordinates           
          r21=r2-r1 ; r31=r3-r1

          ax=(r21(2,1)*r31(3,1)-r21(3,1)*r31(2,1))*0.5
          ay=(r21(3,1)*r31(1,1)-r21(1,1)*r31(3,1))*0.5
          az=(r21(1,1)*r31(2,1)-r21(2,1)*r31(1,1))*0.5
          area=SQRT(ax**2+ay**2+az**2)

          tri(tr)%ar=area ; ax=ax/area ; ay=ay/area ; az=az/area
          tri(tr)%fnor=RESHAPE((/ax,ay,az/),SHAPE(tri(tr)%fnor))

          tri(tr)%vol=(r1(1,1)*(r2(2,1)*r3(3,1)-r2(3,1)*r3(2,1))+ r1(2,1)*(r2(3,1)*r3(1,1)-r2(1,1)*r3(3,1))+ &
		       r1(3,1)*(r2(1,1)*r3(2,1)-r2(2,1)*r3(1,1)))/6.0

          RETURN
          END SUBROUTINE onlyarea

!-------------------------------------------------------------------------------------------------------------------------
!                                   FUNCTION TO CHECK THE MINIMUM BONDLENGTH
!-------------------------------------------------------------------------------------------------------------------------
      FUNCTION Minimum_bondlength_check(vert,vercell) result(fres) 
      USE module_datastruct 
      IMPLICIT NONE
      REAL(KIND=8)::bl
      REAL(KIND=8),DIMENSION(3,1) ::co1,co2     
      INTEGER::vert,ver2,fres,i,cellnumber,vercell,k												     ! fres=1 ===>false ,0===>true
       fres=0 ; k=0
       Do i=1,27,1
       cellnumber=cell(vercell)%neigh_cell(i)
       ver2=cell(cellnumber)%first_vert
       Do while(ver2.Ne.-1)
       If (ver2 .Ne. vert) Then
       co1=ver(vert)%vcoord ; co2=ver(ver2)%vcoord ; bl=Sum((co2-co1)**2)
       If(bl .LE. 1.00)Then 
	      fres=1  ; Return ; Endif   
       Endif	      
       ver2=ver(ver2)%next_vert
       Enddo
       Enddo
       Return
      END FUNCTION Minimum_bondlength_check
!-------------------------------------------------------------------------------------------------------------------------
!               FUNCTION TO CHECK THE MAX AND MIN BONDLENGTH EITHER FOR TWO VERTICES (vert and ver2) OR
!               WITH ALL VERTICES IN THE LINK CELLS AROUND vert
!-------------------------------------------------------------------------------------------------------------------------
      FUNCTION Neigh_blcheck(vert,vert1) result(fres) 
      USE module_datastruct 
      IMPLICIT NONE
      Integer,Optional,Intent(in) ::vert1
      REAL(KIND=8)::bl
      INTEGER::i,ver1,vert,fres                                            !fres=1 ===>false ,0===>true
      fres=0 
      If(Present(vert1))Then
      bl=Sum((ver(vert)%vcoord-ver(vert1)%vcoord)**2)
      IF ((bl .Ge. 3.0) .Or. (bl .Le. 1.0)) fres=1
      Else
      Do i=1,ver(vert)%nonei,1
      ver1=ver(vert)%vneipt(i)
      bl=Sum((ver(vert)%vcoord-ver(ver1)%vcoord)**2)
      IF ((bl .Ge. 3.0) .Or. (bl .Le.1.0)) Then
      fres=1 ; RETURN ; Endif
      Enddo
      Endif
      RETURN
!----------------------------------------------------------------------------------------------------------------------          
      END FUNCTION Neigh_blcheck


      END MODULE module_mcsmoves

