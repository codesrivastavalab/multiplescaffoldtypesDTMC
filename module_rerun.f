!===========================================================================================================================
!===========================================================================================================================
!				MODULE TO START FROM AN ARBITRARY CONFIGURATION
!===========================================================================================================================
!===========================================================================================================================

       MODULE module_rerun
       contains 
!---------------------------------------------------------------------------------------------------------------------------
!                 Subroutine to calculate the area of the given face and update the total area linked to the vertices
!---------------------------------------------------------------------------------------------------------------------------

       SUBROUTINE areacal(tr)
       USE module_datastruct 
       IMPLICIT NONE
       INTEGER::i,j,k,tr
       REAL(KIND=8) ::ax,ay,az,area
       REAL(KIND=8),DIMENSION(3,1)::r1,r2,r3,r21,r31
       ax=0 ; ay=0 ; az=0 ; area=0

       i=tri(tr)%vert(1) ; j=tri(tr)%vert(2) ; k=tri(tr)%vert(3)
       r1=ver(i)%vcoord ; r2=ver(j)%vcoord ; r3=ver(k)%vcoord                                                  ! The coordinates of the 3 vertices      
       r21=r2-r1 ; r31=r3-r1                                                                                   ! Relative position vectors for area calc
       ver(i)%totarea=ver(i)%totarea-tri(tr)%ar/3.0                                                            ! Contrib of exisiting area to totalarea 
       ver(j)%totarea=ver(j)%totarea-tri(tr)%ar/3.0
       ver(k)%totarea=ver(k)%totarea-tri(tr)%ar/3.0

        ax=(r21(2,1)*r31(3,1)-r21(3,1)*r31(2,1))*0.5
        ay=(r21(3,1)*r31(1,1)-r21(1,1)*r31(3,1))*0.5
        az=(r21(1,1)*r31(2,1)-r21(2,1)*r31(1,1))*0.5
        area=SQRT(ax**2+ay**2+az**2)
        tri(tr)%ar=area                                       
        ax=ax/area ; ay=ay/area ;  az=az/area                                                                  ! Normalized area vector 

       tri(tr)%fnor=RESHAPE((/ax,ay,az/),SHAPE(tri(tr)%fnor)) 
       ver(i)%totarea=ver(i)%totarea+tri(tr)%ar/3.0                                                            ! Contrib of new area to the total area of vert
       ver(j)%totarea=ver(j)%totarea+tri(tr)%ar/3.0
       ver(k)%totarea=ver(k)%totarea+tri(tr)%ar/3.0

       tri(tr)%vol=(r1(1,1)*(r2(2,1)*r3(3,1)-r2(3,1)*r3(2,1))+ r1(2,1)*(r2(3,1)*r3(1,1)-r2(1,1)*r3(3,1))+ &
                    r1(3,1)*(r2(1,1)*r3(2,1)-r2(2,1)*r3(1,1)))/6.0


       END SUBROUTINE areacal     
!--------------------------------------------------------------------------------------------------------------------------
!                                     SUBROUTINE TO MAKE A RESTART FILE
!--------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE Write_Memb_Configuration(time)
        USE module_datastruct
        INTEGER:: i,time,Fin
        Character(100) :: nname,fname,tname

318    FORMAT(10(I4,1X))
320    FORMAT(F9.6,1X,F9.6,1X,F9.6,1X,I2)

        Write(nname,*),mynod
        Write(tname,*),time
        fname="./rundir-"//Trim(Adjustl(nname))//'/startdet-'//Trim(Adjustl(tname))//'.dat'
        Fin=11+mynod

	  OPEN(Fin,FILE=Trim(Adjustl(fname)))
          WRITE(Fin,*),nver,"vertexdet"   
          DO i=1,nver,1
          WRITE(Fin,*)ver(i)%vcoord
          WRITE(Fin,*)ver(i)%nonei
          WRITE(Fin,318)ver(i)%vneipt
          WRITE(Fin,318)ver(i)%vneitr
          ENDDO

          WRITE(Fin,*),ntr,"triangledet"
          DO i=1,ntr
          WRITE(Fin,*)tri(i)%vert
          WRITE(Fin,*)tri(i)%li
          ENDDO
   
          WRITE(Fin,*)tlink,"linkdet"
          DO i=-tlink,tlink,1
             IF(i.NE.0)THEN
             WRITE(FIn,*)lin(i)%sep  
             WRITE(Fin,*)lin(i)%tr  
             ENDIF
          ENDDO

          Write(Fin,*),nver,"nematicdet"      
          DO i=1,nver,1
          WRITE(Fin,320)ver(i)%splo,ver(i)%phase
          ENDDO

          Write(Fin,*),nver,"nematicdet-glo"      
          DO i=1,nver,1
          WRITE(Fin,320)ver(i)%spgl,ver(i)%phase
          ENDDO

          Write(Fin,*),nver,"c1  c2 mcur"      
          DO i=1,nver,1
          WRITE(Fin,*)ver(i)%cur1,ver(i)%cur2,ver(i)%mcur
          ENDDO

          Call MPI_Barrier(nallgrp,ier)	  
          CLOSE(Fin)

        END SUBROUTINE Write_Memb_Configuration

!---------------------------------------------------------------------------------------------------------------------------
!                              To start  from one compact startupfile
!---------------------------------------------------------------------------------------------------------------------------

       SUBROUTINE Read_Memb_Configuration(state)
       USE module_datastruct 
       INTEGER::i
       Character(100) :: state,nname,fname

		Write(nname,*),mynod
		fname='./startdet.in'
		OPEN(11,FILE=Trim(Adjustl(fname)),FORM='FORMATTED')
        READ(11,*),nver                                                ! Number of vertex
	  	Allocate(ver(nver))
          DO i=1,nver,1
          READ(11,*)ver(i)%vcoord                                        ! Position of coordinates 
          READ(11,*)ver(i)%nonei                                         ! Position of coordinates 
          READ(11,*)ver(i)%vneipt                                        ! Neighbouring vertex     
          READ(11,*)ver(i)%vneitr                                        ! Neighbouring triangles  
          ENDDO

          READ(11,*),ntr
		  Allocate(tri(ntr)) 
          DO i=1,ntr
          READ(11,*)tri(i)%vert
          READ(11,*)tri(i)%li
          ENDDO

          READ(11,*)tlink
	  Allocate(lin(-tlink:tlink))
          DO i=-tlink,tlink,1
             IF(i.NE.0)THEN
             READ(11,*)lin(i)%sep
             READ(11,*)lin(i)%tr
             ENDIF
          ENDDO

          If(Trim(state).Eq.'RESTART')Then                               ! Read the nematic data only if called for 
          READ(11,*)  
          DO i=1,nver,1
          READ(11,*)ver(i)%splo,ver(i)%phase
          ENDDO
          ElseIf(Trim(state).Eq.'MEMBRANE')Then                          ! Do not read the nematic data otherwise
          ver(1:nver)%splo(1,1)=0.0
          ver(1:nver)%splo(2,1)=0.0
          ver(1:nver)%splo(3,1)=0.0
          Endif

          CLOSE(11)

          DO i=1,ntr,1                                                   ! Calculate the area of the triangles 
          Call areacal(i)
          ENDDO
       
          END SUBROUTINE Read_Memb_Configuration

        END MODULE module_rerun

