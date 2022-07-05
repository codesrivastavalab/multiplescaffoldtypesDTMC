! the first and last cell in every direction is a buffer. The linklist should be 
! reinitialized everytime a vertex move pushes a membrane vertex into one of the 
! buffer cells. This will work in the present case since the conformations of
! the membrane is entriely dependent upon its relative position and hence
! recentering the membrane at the center of a linkcell should not affect its
! thermodynamic behaviour in any way.

! The linkcell size should be reduced and made equal in all directions if we
! plan to implement MPI parallelization of the code. In the present form most
! of the spatial decomposition would contain no data. Reducing the number of
! cells would lead to high reboxing frequency which is a  trade off to the 
! overhead arising from the number of null data cells.


!============================================================================================
!                              MODULE THAT HANDLES ALL LINKLIST RELATED  CALCULATIONS
!============================================================================================
	Module module_linklist_calc
	Contains
!---------------------------------------------------------------------------------------------
! 			Subroutine to Initialize the link list variables
!---------------------------------------------------------------------------------------------
	Subroutine Initialize_linklist()
	Use Module_datastruct ; Use Module_dataformat
	Implicit None
	Real(Kind=8) :: xmax,xmin,ymax,ymin,zmax,zmin
	Integer :: i
	Real(Kind=8) :: xcoord(nver),ycoord(nver),zcoord(nver)
	Real(Kind=8) :: rcom(3,1),rbox(3,1),shift_vec(3,1)
	cell_lenx=3 ; cell_leny=3; cell_lenz=3
	xcoord=reshape(ver(:)%vcoord(1,1),shape(xcoord))
	ycoord=reshape(ver(:)%vcoord(2,1),shape(ycoord))
	zcoord=reshape(ver(:)%vcoord(3,1),shape(zcoord))
	rcom=reshape((/Sum(xcoord),Sum(ycoord),Sum(zcoord)/),shape(rcom))/nver
	xmax=Maxval(xcoord); xmin=Minval(xcoord)
	ymax=Maxval(ycoord); ymin=Minval(ycoord)
	zmax=Maxval(zcoord); zmin=Minval(zcoord)

	Ncellx=2*Nint((xmax-xmin)/cell_lenx)+2                                                                                      ! Cell number 1 and Ncellx are buffers (reinitialize when membrane reaches here)
	Ncelly=2*Nint((ymax-ymin)/cell_leny)+2
	Ncellz=2*Nint((zmax-zmin)/cell_lenz)+2
	Ncell=Ncellx*Ncelly*Ncellz
	rbox=reshape((/Ncellx*cell_lenx/2,Ncelly*cell_leny/2,Ncellz*cell_lenz/2/),shape(rbox))
	shift_vec=rbox-rcom                                                                                                        ! Shift the membrane by shift_vec

        ver(:)%vcoord(1,1)=ver(:)%vcoord(1,1)+shift_vec(1,1)
        ver(:)%vcoord(2,1)=ver(:)%vcoord(2,1)+shift_vec(2,1) 
        ver(:)%vcoord(3,1)=ver(:)%vcoord(3,1)+shift_vec(3,1)

	If(Allocated(cell)) Deallocate(cell)
	Allocate(cell(Ncell))
	cell(:)%first_vert=-1 ; ver(:)%linkcell=-1 ; ver(:)%next_vert=-1 ; cell(:)%last_vert=-1                                     ! Initialize all linklist related varibles in various structures 
	cell(:)%num_vert=0  ; ver(:)%prev_vert=-1 ; cell(:)%buffer=0
	Call Initialize_buffer()
	Do i=1,nver,1
	Call Compute_cellnumber(i)
	Enddo
	Call Initialize_linkcells()
	Call Write_linkcells(0)
	End Subroutine Initialize_linklist
!---------------------------------------------------------------------------------------------
! 			Subroutine to calculate the linklist for a given vertex
!---------------------------------------------------------------------------------------------
	Subroutine Compute_cellnumber(vert)
	Use Module_datastruct 
	Implicit None
	Integer :: vert,xcell,ycell,zcell,cellnumber
	xcell=INT(ver(vert)%vcoord(1,1)/cell_lenx)+1
	ycell=INT(ver(vert)%vcoord(2,1)/cell_leny)+1
	zcell=INT(ver(vert)%vcoord(3,1)/cell_lenz)+1
	cellnumber=xcell+(ycell-1)*Ncellx+(zcell-1)*Ncellx*Ncelly                                                                              ! The absolute number of the cell
	ver(vert)%linkcell=cellnumber
	if(cell(cellnumber)%first_vert==-1)Then
		cell(cellnumber)%first_vert=vert                                                                                       ! the first vertex is vert 
		cell(cellnumber)%last_vert=vert                                                                             	       ! The last vertex is vert
		ver(vert)%prev_vert=-1
		ver(vert)%next_vert=-1
		cell(cellnumber)%num_vert=cell(cellnumber)%num_vert+1
	Else
		ver(cell(cellnumber)%last_vert)%next_vert=vert                                                                         ! vert is linked to the last vertex in the cell 
		ver(vert)%next_vert=-1
		ver(vert)%prev_vert=cell(cellnumber)%last_vert
		cell(cellnumber)%last_vert=vert                                                                                        ! now vert becomes the last vertex 
		cell(cellnumber)%num_vert=cell(cellnumber)%num_vert+1                                                        	       ! number of vertices in the cell in incremented by 1
	Endif
	End Subroutine Compute_cellnumber


!---------------------------------------------------------------------------------------------
!		Subroutine to assign buffer  to	linkcells
!---------------------------------------------------------------------------------------------
	Subroutine Initialize_buffer()
	Use Module_datastruct 
	Implicit None
	Integer :: i,j,k,cellno
	Do k=1,Ncellz,1
	Do j=1,Ncelly,1
	Do i=1,Ncellx,1
	cellno=i+(j-1)*Ncellx+(k-1)*Ncellx*Ncelly
	If((i.eq.1) .Or. (i.eq.Ncellx)) cell(cellno)%buffer=1
	If((j.eq.1) .Or. (j.eq.Ncelly)) cell(cellno)%buffer=1
	If((k.eq.1) .Or. (k.eq.Ncellz)) cell(cellno)%buffer=1
	Enddo ; Enddo ; Enddo
        End Subroutine Initialize_buffer
!---------------------------------------------------------------------------------------------
!		Subroutine to assign coordinates and neighbour list to	linkcells
!---------------------------------------------------------------------------------------------
	Subroutine Initialize_linkcells()
	Use Module_datastruct 
	Implicit None
	Integer :: i,j,k,cellno,neig_no,ix,iy,iz,nei_cell
	Do k=1,Ncellz,1
	Do j=1,Ncelly,1
	Do i=1,Ncellx,1
	cellno=i+(j-1)*Ncellx+(k-1)*Ncellx*Ncelly
	cell(cellno)%start_coord=Reshape((/(i-1)*cell_lenx,(j-1)*cell_leny,(k-1)*cell_lenz/),(/3,1/)) 
	cell(cellno)%end_coord=Reshape((/i*cell_lenx,j*cell_leny,k*cell_lenz/),(/3,1/)) 

	neig_no=0
	Do ix=-1,1,1
	Do iy=-1,1,1
	Do iz=-1,1,1
	nei_cell=cellno+ix+iy*Ncellx+iz*Ncellx*Ncelly
        neig_no=neig_no+1
	cell(cellno)%neigh_cell(neig_no)=nei_cell
	Enddo ; Enddo ;Enddo

	Enddo ; Enddo ; Enddo
        End Subroutine Initialize_linkcells

!---------------------------------------------------------------------------------------------
! 			Subroutine to calculate the linklist for a given vertex
!---------------------------------------------------------------------------------------------
	Function get_cellnumber(vert) Result(cellnumber)
	Use Module_datastruct
	Implicit None
	Integer :: vert,xcell,ycell,zcell,cellnumber
	xcell=INT(ver(vert)%vcoord(1,1)/cell_lenx)+1
	ycell=INT(ver(vert)%vcoord(2,1)/cell_leny)+1
	zcell=INT(ver(vert)%vcoord(3,1)/cell_lenz)+1
	cellnumber=xcell+(ycell-1)*Ncellx+(zcell-1)*Ncellx*Ncelly                                                                              ! The absolute number of the cell
	Return
	End Function get_cellnumber

!---------------------------------------------------------------------------------------------
! 			Subroutine to calculate the linklist for a given vertex
!---------------------------------------------------------------------------------------------
	Subroutine update_vertex_cell(vert,new_cell) 
	Use Module_datastruct
	Implicit None
	Integer :: vert,cur_cell,new_cell,prev_vert,next_vert
	cur_cell=ver(vert)%linkcell
	prev_vert=ver(vert)%prev_vert
	next_vert=ver(vert)%next_vert 

	if((prev_vert .Ne.-1) .And. (next_vert.Ne.-1)) Then                                                                		    ! in the middle of the link list
		ver(prev_vert)%next_vert=next_vert                                                                              	    ! Relink vertex linked by vert (%next_vert) to its prev_vert(%prev_vert)
		ver(next_vert)%prev_vert=prev_vert                   
	else if((prev_vert .Eq.-1) .And. (next_vert.Ne.-1))Then                                                                             ! first vertex in the link list
	       cell(cur_cell)%first_vert=next_vert
	       ver(next_vert)%prev_vert=-1
        else if((next_vert.Eq.-1) .And. (prev_vert.Ne.-1))Then                                                                              ! Last vertex in the link list 
	       cell(cur_cell)%last_vert=prev_vert                                        	       
	       ver(prev_vert)%next_vert=-1
       else if ((next_vert.Eq.-1) .And. (prev_vert.Eq.-1))Then
	       cell(cur_cell)%last_vert=-1
	       cell(cur_cell)%first_vert=-1
	Endif
        cell(cur_cell)%num_vert=cell(cur_cell)%num_vert-1                                                                                   ! Reduce the current cell occupants by 1


	If(cell(new_cell)%first_vert .Eq. -1) Then
		cell(new_cell)%first_vert=vert
		cell(new_cell)%last_vert=vert
		ver(vert)%prev_vert=-1
		ver(vert)%next_vert=-1
	Else
		ver(cell(new_cell)%last_vert)%next_vert=vert
		ver(vert)%prev_vert=cell(new_cell)%last_vert
		ver(vert)%next_vert=-1
		cell(new_cell)%last_vert=vert
	Endif
	ver(vert)%linkcell=new_cell
	cell(new_cell)%num_vert=cell(new_cell)%num_vert+1  

	End Subroutine update_vertex_cell

!---------------------------------------------------------------------------------------------
! 			Subroutine to add a  vertex(vert) to a cell(new_cell)
!---------------------------------------------------------------------------------------------
	Subroutine addvertex_tocell(vert,new_cell) 
	Use Module_datastruct
	Implicit None
	Integer :: vert,new_cell

	If(cell(new_cell)%first_vert .Eq. -1) Then
		cell(new_cell)%first_vert=vert
		cell(new_cell)%last_vert=vert
		ver(vert)%prev_vert=-1
		ver(vert)%next_vert=-1
	Else
		ver(cell(new_cell)%last_vert)%next_vert=vert
		ver(vert)%prev_vert=cell(new_cell)%last_vert
		ver(vert)%next_vert=-1
		cell(new_cell)%last_vert=vert
	Endif
	ver(vert)%linkcell=new_cell
	cell(new_cell)%num_vert=cell(new_cell)%num_vert+1  
	End Subroutine addvertex_tocell


!---------------------------------------------------------------------------------------------
! 			Subroutine to calculate the linklist for a given vertex
!---------------------------------------------------------------------------------------------
	Subroutine removevertex_fromcell(vert) 
	Use Module_datastruct
	Implicit None
	Integer :: vert,cur_cell,prev_vert,next_vert
	cur_cell=ver(vert)%linkcell
	prev_vert=ver(vert)%prev_vert
	next_vert=ver(vert)%next_vert 

	if(prev_vert .Ne.-1 .And. next_vert.Ne.-1) Then                                                                		    ! in the middle of the link list
	ver(prev_vert)%next_vert=next_vert                                                                              	    ! Relink vertex linked by vert (%next_vert) to its prev_vert(%prev_vert)
	ver(next_vert)%prev_vert=prev_vert                   
	else if(prev_vert .Eq.-1)Then                                                                                               ! first vertex in the link list
	       cell(cur_cell)%first_vert=next_vert
	       ver(next_vert)%prev_vert=-1
        else if(next_vert.Eq.-1)Then                                                                                                ! Last vertex in the link list 
	       cell(cur_cell)%last_vert=prev_vert                                        	       
	       ver(prev_vert)%next_vert=-1
	Endif
        cell(cur_cell)%num_vert=cell(cur_cell)%num_vert-1                                                                           ! Reduce the current cell occupants by 1
	ver(vert)%linkcell=-1
	End Subroutine removevertex_fromcell
	
	End Module module_linklist_calc
