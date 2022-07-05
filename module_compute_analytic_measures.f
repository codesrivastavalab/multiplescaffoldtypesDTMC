! This module contains subroutines to compute various measures  used in the analysis of the membranes

          Module module_compute_analytic_measures
          Contains
!---------------------------------------------------------------------------------------------------------------------------
!                           SUBROUTINE TO CALCULATE ALL ANALYTIC  QUANTITIES
!---------------------------------------------------------------------------------------------------------------------------
         SUBROUTINE analqtys()
         USE module_datastruct ; USE module_curvcalc ; USE module_energy_calculations
         IMPLICIT NONE
         INTEGER :: i,ilen
         REAL(KIND=8) :: xcm,ycm,zcm,elen,neneen,necuen,Isinen        ! COM, Elastic,nematic-nematic,nematic-curvature 
         REAL(KIND=8) :: imcur,imcursq,rcm(3,1)
         REAL(KIND=8) :: rgeif(3,3),ineif(3,3) 
         
         mp%rg=0 ; mp%area=0 ; ilen=0 ; imcur=0. ;imcursq=0.
         rgeif=0.0 ; ineif=0.0 ; Isinen=0.0

          xcm=SUM(ver(1:nver)%vcoord(1,1))/nver                          ! X component of the center of mass
          ycm=SUM(ver(1:nver)%vcoord(2,1))/nver                          ! Y component                      
          zcm=SUM(ver(1:nver)%vcoord(3,1))/nver                          ! Z component                      
          rcm=reshape((/xcm,ycm,zcm/),shape(rcm)) 

         mp%rg=(SUM((ver(1:nver)%vcoord(1,1)-xcm)**2)+SUM((ver(1:nver)%vcoord(2,1)-ycm)**2)+ &
	 &			SUM((ver(1:nver)%vcoord(3,1)-zcm)**2))/nver                                         ! Rg square calculation 

         mp%area=SUM(tri(1:ntr)%ar)                                                                 ! Surface area of the membrane
         elen=0 ; necuen=0 ; neneen=0

         DO i=1,nver
         elen=elen+Helfrich_energy(i)                                                              ! Total elastic bending energy       
	     necuen=necuen+Nematic_Curvature_Energy(i)
         neneen=neneen+0.5*Nematic_Interaction_energy(i)
         ENDDO

         Call volumecalc()                                                                          ! Call volume calculations 
		 elen = elen - pr*Sum(tri(:)%vol)
         mp%enel=elen ;  mp%ennc = necuen ;  mp%ennn = neneen

	 	If ( Trim(Adjustl(mp%Lambda_variable)) .Eq. 'EPSILON') Then
		 mp%lam_ennn = mp%lambda * mp%ennn
		 mp%lam_ennc = mp%ennc
	 	Else 
		 mp%lam_ennc = mp%lambda*mp%ennc
		 mp%lam_ennn = mp%ennn
		Endif

       END SUBROUTINE analqtys
!--------------------------------------------------------------------------------------------------------------------------------------
!                               SUBROUTINE TO CALCULATE THE VOLUME OF THE GEOM
!--------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE volumecalc()
      USE module_datastruct ; USE module_curvcalc 
      USE module_mcsmoves
      IMPLICIT NONE
      INTEGER :: i
      REAL(KIND=8) :: vol,cvol                                       
      REAL(KIND=8),DIMENSION(3,1) :: a,b,c    
      mp%vol=0.0; vol=0.0

      vol_calc: DO i=1,ntr
      a=0; b=0; c=0
      a=ver(tri(i)%vert(1))%vcoord  
      b=ver(tri(i)%vert(2))%vcoord  
      c=ver(tri(i)%vert(3))%vcoord
      cvol=(a(1,1)*(b(2,1)*c(3,1)-b(3,1)*c(2,1))+a(2,1)*(b(3,1)*c(1,1)-b(1,1)*c(3,1))+a(3,1)*(b(1,1)*c(2,1)-b(2,1)*c(1,1)))/6.0	    ! (a.(bxc))/6.0 is the volume 
      vol=vol+cvol 
      ENDDO vol_calc
      mp%vol=vol
      END SUBROUTINE volumecalc

!------------------------------------------------------------------------------------------------------------------------
!                 SUBROUTINE TO CALCULATE THE LOCAL AND AVERAGE ORDER PARAMETER ( The winding number at each vertex)
!------------------------------------------------------------------------------------------------------------------------
          SUBROUTINE orderparameter()
          USE module_datastruct;  USE module_curvcalc
          IMPLICIT NONE
          INTEGER :: i,j,v,ip1,jp1
          REAL(KIND=8) :: tang,theta1,theta2,dtheta,defang           ! The value of cos(\theta) and temporary storage      
          REAL(KIND=8),DIMENSION(3,3):: tapl1,tapl2                      ! tangent planes 1 and 2                              
          REAL(KIND=8),DIMENSION(3,1):: geod,pg1,pg2,cp1,cp2,n1,n2       ! Geodesic, projections 1,2, cross prod 1,2           
          REAL(KIND=8),DIMENSION(3,1):: ln1,ln2
          REAL(KIND=8),DIMENSION(1,1):: c1,c2,s1,s2,fang
  
           DO v=1,nver,1                                                 ! Order parameter calculation at each vertex                               
             tang=0 ; defang=0                                           ! Temporary variable to store the total OP, deficit angle at each vertex   
            DO i=1,ver(v)%nonei,1                                        ! sum over each neighbour                                                  
              ip1=i+1 ; IF(i.EQ.ver(v)%nonei) ip1=1
              j=ver(v)%vneipt(i) ; jp1=ver(v)%vneipt(ip1)                ! Choose two adjacent neighbours in ACW direction
              ln1=ver(j)%vcoord-ver(v)%vcoord                            ! Edge vector joining v to j    
              ln2=ver(jp1)%vcoord-ver(v)%vcoord                          ! Edge vector joining v to jp1  
              ln1=ln1/SQRT(SUM(ln1**2)) ; ln2=ln2/SQRT(SUM(ln2**2))      ! Normalize ln1 and ln2         
              fang=MATMUL(TRANSPOSE(ln1),ln2)                            ! Compute the angle at vertex v for face { v,j,jp1 }
              defang=defang+acos(fang(1,1))                              ! Deficit angle summed over all neighbouring faces  
              
              n1=ver(j)%vnor ; n2=ver(jp1)%vnor
              geod=ver(jp1)%vcoord-ver(j)%vcoord                         ! Geodesic from vertex j to jp1
              geod=geod/SQRT(SUM(geod**2))
              tapl1=unitmat-MATMUL(n1,TRANSPOSE(n1))                     ! Tangent plane at vertex v
              tapl2=unitmat-MATMUL(n2,TRANSPOSE(n2))                     ! Tangent plane at vertex j
              pg1=MATMUL(tapl1,geod)                                     ! Projection of geod on tang plane 
              pg1=pg1/SQRT(SUM(pg1**2))
              pg2=MATMUL(tapl2,geod)                                   
              pg2=pg2/SQRT(SUM(pg2**2))

              c1=MATMUL(TRANSPOSE(pg1),ver(j)%spgl)                      ! cosine of angle of spin1 with the geodesic (x comp. at vertex j)  
              c2=MATMUL(TRANSPOSE(pg2),ver(jp1)%spgl)                    ! cosine of angle of spin2 with the geodesic (x comp. at vertex jp1)

              cp1(1,1)=n1(2,1)*pg1(3,1)-n1(3,1)*pg1(2,1)
              cp1(2,1)=n1(3,1)*pg1(1,1)-n1(1,1)*pg1(3,1) 
              cp1(3,1)=n1(1,1)*pg1(2,1)-n1(2,1)*pg1(1,1) 
              cp1=cp1/SQRT(SUM(cp1**2)) 
              s1=MATMUL(TRANSPOSE(cp1),ver(j)%spgl)                      ! sine of the angle between spin1 and vector perp. to geodesic ( y comp at vertex j) 
              
              cp2(1,1)=n2(2,1)*pg2(3,1)-n2(3,1)*pg2(2,1)
              cp2(2,1)=n2(3,1)*pg2(1,1)-n2(1,1)*pg2(3,1) 
              cp2(3,1)=n2(1,1)*pg2(2,1)-n2(2,1)*pg2(1,1) 
              cp2=cp2/SQRT(SUM(cp2**2))
              s2=MATMUL(TRANSPOSE(cp2),ver(jp1)%spgl)                    ! sine of the angle between spin2 and vector perp. to geodesic ( y comp at vertex jp1) 
              theta1=atan2(s1(1,1),c1(1,1))                             
              theta2=atan2(s2(1,1),c2(1,1))
              IF(theta1.GT. pi/2.)  theta1=theta1-pi                      ! No all negative angles are rotated by 180 degrees
              IF(theta1.LT. -pi/2.) theta1=theta1+pi 
              IF(theta2.GT. pi/2.)  theta2=theta2-pi                      ! No all negative angles are rotated by 180 degrees
              IF(theta2.LT. -pi/2.) theta2=theta2+pi 
              dtheta=theta2-theta1
              IF(dtheta.LT.-pi/2.)dtheta=dtheta+pi 
              IF(dtheta.GT.pi/2.) dtheta=dtheta-pi 
              tang=tang+dtheta
            ENDDO
              ver(v)%op=tang/(2*pi)                                      ! Total winding number at a vertex v 
           ENDDO
          END SUBROUTINE orderparameter
          End Module module_compute_analytic_measures

