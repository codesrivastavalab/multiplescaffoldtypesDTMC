        Module module_dataformat
        Contains
!------------------------------------------------------------------------------------------------------------------------
!                                    SUBROUTINE TO Write THE FILES FOR VTK VIEWERS
!------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE vtkformat(time)
        Use module_datastruct
        Integer::time,i,Fin
        Integer:: foff(ntr),ftype(ntr)
        Character(100) :: nname,fname,tname
        
2000    FORMAT(A23,I5,A1,1X,A15,I5,A2)
     
        Write(nname,*),mynod
        Write(tname,*),time
        fname="./rundir-"//Trim(Adjustl(nname))//'/conf-'//Trim(Adjustl(tname))//'.vtu'
        Fin=11+mynod
        Open(Fin,FILE=Trim(Adjustl(fname)))
        Write(Fin,*)'<VTKFile type="UnstructuredGrid" version="0.1"  byte_order="LittleEndian">'
        Write(Fin,*)'<UnstructuredGrid>'
        Write(Fin,2000)'<Piece NumberOfPoints="',nver,'"','NumberOfCells="',ntr,'">'
        Write(Fin,*)'<PointData Scalars="scalars">'

        Write(Fin,*)'<DataArray type="Float32" Name="H" Format="ascii">'  							    ! Mean curvature
        Do i=1,nver,1
        Write(Fin,*) ver(i)%mcur/2.0
        Enddo
        Write(Fin,*)'</DataArray>'

		Write(Fin,*)'<DataArray type="Int32" Name="linkcell" Format="ascii">'  					 
        Do i=1,nver,1
        Write(Fin,*) ver(i)%linkcell
        Enddo
        Write(Fin,*)'</DataArray>'

       Write(Fin,*)'<DataArray type="Int32" Name="ph" Format="ascii">'  							    ! Winding number 
        Do i=1,nver,1
        Write(Fin,*) ver(i)%phase
        Enddo
        Write(Fin,*)'</DataArray>'

       Write(Fin,*)'<DataArray type="Float32" Name="Vert-field"  NumberOfComponents="3"  Format="ascii">'
        Do i=1,nver,1
        Write(Fin,*) ver(i)%spgl
        Enddo
        Write(Fin,*)'</DataArray>'

        Write(Fin,*)'</PointData>'

        Write(Fin,*)'<Points>'
        Write(Fin,*)'<DataArray type="Float32"  NumberOfComponents="3" Format="ascii">'
        Do i=1,nver,1
        Write(Fin,*),ver(i)%vcoord
        Enddo
        Write(Fin,*)'</DataArray>'
        Write(Fin,*)'</Points>'

        Write(Fin,*)'<Cells>'
        Write(Fin,*)'<DataArray type="Int32"  Name="connectivity" Format="ascii">'
        Do i=1,ntr,1
        Write(Fin,*),tri(i)%vert-1
        Enddo
        Write(Fin,*)'</DataArray>'

        Do i=1,ntr,1
        foff(i)=3*i
        Enddo
        ftype=5

        Write(Fin,*)'<DataArray type="Int32" Name="offsets"  Format="ascii">' 
        Write(Fin,*)foff
        Write(Fin,*)'</DataArray>'

        Write(Fin,*)'<DataArray type="Int32" Name="types"  Format="ascii">'
        Write(Fin,*)ftype
        Write(Fin,*)'</DataArray>'
        Write(Fin,*)'</Cells>'
        Write(Fin,*)'</Piece>'
        Write(Fin,*)'</UnstructuredGrid>'
        Write(Fin,*)'</VTKFile>'
        Call MPI_Barrier(nallgrp,ier)
        Close(Fin)
        End Subroutine vtkformat


!------------------------------------------------------------------------------------------------------------------------
!             SUBROUTINE TO Write THE FILES FOR JAVAVIEW (Yet another viewing format -- use Javaview (www.javaview.de))
!------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE javaviewdat(time)
        USE module_datastruct
        Integer::time,i,j
        Character(100) :: nname,fname,tname
        Integer:: C(Nphase,3)                                            ! Vector field color stored for upto 6 phases

        Do i=1,Nphase,1
        Select Case(i)
        Case (1)
                C(1:1,:)=reshape((/255,0,0/),shape(C(1:1,1:3)))
        Case (2)
                C(2:2,:)=reshape((/0,0,255/),shape(C(1:1,1:3)))
        Case (3)
                C(3:3,:)=reshape((/0,255,0/),shape(C(1:1,1:3)))
        Case (4) 
                C(4:4,:)=reshape((/225,255,0/),shape(C(1:1,1:3)))
        Case (5) 
                C(5:5,:)=reshape((/255,0,255/),shape(C(1:1,1:3)))
        Case (6)         
                C(6:6,:)=reshape((/0,255,255/),shape(C(1:1,1:3)))
        End Select
        Enddo

1111    FORMAT(A3,1X,F19.12,1X,F19.12,1X,F19.12,1X,A4)
1112    FORMAT(A3,1X,I4,1X,I4,1X,I4,1X,A4)

        Write(nname,*),mynod
        Write(tname,*),time
        fname="./rundir-"//Trim(Adjustl(nname)) //'/jv-'//Trim(Adjustl(tname))//'.jvx'

        Open(11,FILE=fname,FORM='FORMATTED')
        Write(11,*)'<?xml version="1.0" encoding="ISO-8859-1" standalone="no"?>'
        Write(11,*)'<!DOCTYPE jvx-model SYSTEM  "http://www.javaview.de/rsrc/jvx.dtd">'
        Write(11,*)'<jvx-model>'
        Write(11,*)'<geometries>'
        Write(11,*)'<geometry name="membrane">'  

         Write(11,*)'<pointSet dim="3" color="show" point="hide">'
         Write(11,*)'<points num="',nver,'">'    
         DO i=1,nver,1
         Write(11,1111)'<p>',ver(i)%vcoord,'</p>'
         ENDDO
         Write(11,*)'<thickness>3.0</thickness>'
         Write(11,*)'</points>'

         Write(11,*)'<colors num="',nver,'">'    
         DO i=1,nver,1
         IF(ver(i)%phase .Eq. 1 )THEN
         Write(11,*)'<c> 255  255 255 </c>'                              ! Positive curvature region is red colored      
         ELSE         
         Write(11,*)'<c> 125  125 125 </c> '                             ! Intermediate curvature region is green colured
         ENDIF   
         ENDDO 
         Write(11,*)'</colors>'


         Write(11,*)'</pointSet>'
     
         Write(11,*)'<faceSet face="show" edge="hide">'
         Write(11,*)'<faces>'
         DO i=1,ntr,1
         Write(11,1112)'<f>',tri(i)%vert-1,'</f>'
         ENDDO
         Write(11,*)'<color> 192  192 192 </color>'
         Write(11,*)'</faces>'
         Write(11,*)'</faceSet>'

         Write(11,*)'<vectorField color="show"arrow="hide"  base="vertex">'                          ! The vector field goes in here
         Write(11,*)'<vectors num="',nver,'">'
         DO i=1,nver,1
         If(ver(i)%phase.Eq.1)Then
         Write(11,1111),'<v>',0.6*ver(i)%spgl,'</v>'
         Else
         Write(11,1111),'<v>',0.0, 0.0, 0.0,'</v>'
         Endif
         ENDDO
         Write(11,*)'<thickness>2.0</thickness>'
         Write(11,*)'</vectors>'

         Write(11,*)'<colors num="',nver,'">'    
         DO i=1,nver,1
         j=ver(i)%phase
         Write(11,*)'<c>',C(j,:),'</c>'                                          ! kappa \parallel supporting spins are red
         ENDDO 
         Write(11,*)'</colors>'
         Write(11,*)'</vectorField>'

         Write(11,*)'<vectorField color="show"arrow="hide"   base="vertex">'                  ! The vector field goes in here
         Write(11,*)'<vectors num="',nver,'">'
         DO i=1,nver,1
         If(ver(i)%phase.Eq.1)Then
         Write(11,1111),'<v>',-0.6*ver(i)%spgl,'</v>'
         Else
         Write(11,1111),'<v>',0.0, 0.0, 0.0,'</v>'
         Endif
         ENDDO
         Write(11,*)'<thickness>2.0</thickness>'
         Write(11,*)'</vectors>'

         Write(11,*)'<colors num="',nver,'">'    
         DO i=1,nver,1
         j=ver(i)%phase
         Write(11,*)'<c>',C(j,:),'</c>'                                   ! kappa \parallel supporting spins are red
         ENDDO 
         Write(11,*)'</colors>'
         Write(11,*)'</vectorField>'

        Write(11,*)'</geometry>'
        Write(11,*)'</geometries>'
        Write(11,*)'</jvx-model>'
        Call MPI_Barrier(Nallgrp,ier)                                    ! To prevent first thread closing the file resulting in incomplete write 
        Close(11)

        END SUBROUTINE javaviewdat

!-----------------------------------------------------------------------------------------------------
!       Write THE VTK FILE FOR THE LINK CELLS (You can display the link cells using the file generated here)
!-----------------------------------------------------------------------------------------------------

        SUBROUTINE Write_linkcells(time)
        Use module_datastruct 
        Integer::time,num,i,j
        Integer,Dimension(:),Allocatable:: foff,ftype 
        Character(1000) :: filename,timename,nname


		Write(timename,*),time
		Write(nname,*),mynod
        filename="./rundir-"//Trim(Adjustl(nname))//'/linklist-'//Trim(Adjustl(timename))//'.vtu'         

        Open(11,FILE=Trim(Adjustl(filename)),FORM='FORMATTED')
        Write(11,*)'<VTKFile type="UnstructuredGrid" version="0.1"', '  byte_order="BigEndian">'
        Write(11,*)'<UnstructuredGrid>'
        Write(11,*)'<Piece NumberOfPoints="',Ncell*8,'"	NumberOfCells="',Ncell*6,'">'

        Write(11,*)'<Points>'
        Write(11,*)'<DataArray type="Float32"','  NumberOfComponents="3" Format="ascii">'
		Do i=1,Ncell,1
        Write(11,*),cell(i)%start_coord(1,1),cell(i)%start_coord(2,1),cell(i)%start_coord(3,1)
		Write(11,*),cell(i)%end_coord(1,1),cell(i)%start_coord(2,1),cell(i)%start_coord(3,1)
		Write(11,*),cell(i)%end_coord(1,1),cell(i)%end_coord(2,1),cell(i)%start_coord(3,1)
		Write(11,*),cell(i)%start_coord(1,1),cell(i)%end_coord(2,1),cell(i)%start_coord(3,1)
		Write(11,*),cell(i)%start_coord(1,1),cell(i)%start_coord(2,1),cell(i)%end_coord(3,1)
		Write(11,*),cell(i)%end_coord(1,1),cell(i)%start_coord(2,1),cell(i)%end_coord(3,1)
		Write(11,*),cell(i)%end_coord(1,1),cell(i)%end_coord(2,1),cell(i)%end_coord(3,1)
		Write(11,*),cell(i)%start_coord(1,1),cell(i)%end_coord(2,1),cell(i)%end_coord(3,1)
		Enddo
        Write(11,*)'</DataArray>'
        Write(11,*)'</Points>'

        Write(11,*)'<Cells>'
        Write(11,*)'<DataArray type="Int32"','  Name="connectivity" Format="ascii">'
		Do j=0,Ncell-1,1
		i=j*8
        Write(11,'(I4,1x,I4,1x,I4,1x,I4)'), i+0, i+1, i+2, i+3
        Write(11,'(I4,1x,I4,1x,I4,1x,I4)'), i+1, i+2, i+6, i+5
        Write(11,'(I4,1x,I4,1x,I4,1x,I4)'), i+2, i+3, i+7, i+6
        Write(11,'(I4,1x,I4,1x,I4,1x,I4)'), i+5, i+6, i+7, i+4
        Write(11,'(I4,1x,I4,1x,I4,1x,I4)'), i+4, i+0, i+3, i+7
        Write(11,'(I4,1x,I4,1x,I4,1x,I4)'), i+4, i+0, i+1, i+5
		Enddo
        Write(11,*)'</DataArray>'

        Allocate(foff(Ncell*6)); Allocate(ftype(Ncell*6))
        num=0
        Do i=1,Ncell*6,1
        num=num+1
        foff(num)=4*num
        Enddo
        ftype=5

        Write(11,*)'<DataArray type="Int32" Name="offsets"  Format="ascii">' 
        Write(11,*)foff
        Write(11,*)'</DataArray>'

        Write(11,*)'<DataArray type="Int32" Name="types"  Format="ascii">' 
        Write(11,*)ftype
        Write(11,*)'</DataArray>'

        Write(11,*)'</Cells>'

        Write(11,*)'</Piece>'
        Write(11,*)'</UnstructuredGrid>'
        Write(11,*)'</VTKFile>'
        Call MPI_Barrier(Nallgrp,ier)                             ! To prevent first thread closing the file resulting in incomplete write 	
        Close(11)
        End Subroutine Write_linkcells


        End Module module_dataformat

