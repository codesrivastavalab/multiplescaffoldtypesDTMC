!===========================================================================================================================
!===========================================================================================================================
!        $modcurv                    MODULE TO CALCULATE THE CURVATURE AND OTHER RELATED QUANTITIES
!===========================================================================================================================
!===========================================================================================================================
       MODULE  module_curvcalc
       IMPLICIT NONE
       REAL(KIND=8),DIMENSION(3,3) :: unitmat,cmat,hmat 
       REAL(KIND=8),DIMENSION(10,3):: tanproj 
       
      contains
   
       SUBROUTINE normalcalc(vernu)
       USE module_datastruct ; USE module_rerun ; Use Omp_lib
       IMPLICIT NONE
       INTEGER :: i,ip1,im1,j,vernu,tr1,tr2                            ! Loop integers
       REAL(KIND=8),DIMENSION(3,1)::nor,rtan,amat,smat,zdir    
       REAL(KIND=8),DIMENSION(3,3)::pmat,evec,g2l                        ! g2l is the global to local frame conversion matrix
       REAL(KIND=8),DIMENSION(1,1)::wei,xc,an  
       REAL(KIND=8)::fwei,p,q,r,s,insq,e1,e2,elen,dih,emcu,si 
       REAL(KIND=8),DIMENSION(3,1) :: en,euv,v1,v2,cpm,fn1,fn2       
       unitmat=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),shape(unitmat)) 
       zdir=unitmat(1:3,3:3) ; nor=0 ; cmat=0 ; hmat=0 ; evec=0
       fwei=0

         weightcalc:DO i=1,ver(vernu)%nonei,1
                  j=ver(vernu)%vneitr(i)
                       fwei=fwei+tri(j)%ar                               ! Total area of all neigh triangles
                    ENDDO weightcalc          

         normal : DO i=1,ver(vernu)%nonei
                     j=ver(vernu)%vneitr(i)
                     nor=nor+(tri(j)%ar/fwei)*tri(j)%fnor
                  ENDDO normal

                  nor=nor/SQRT(SUM(nor**2))                              ! Normalize the calculated normal.( Global Frame)
                  ver(vernu)%vnor=nor
                  pmat=(unitmat-MATMUL(nor,TRANSPOSE(nor)))              ! Tangent plane corresponding to normal at vernu 


       !$OMP Parallel Do Reduction(+:cmat)  Private(im1,ip1,v2,elen,tr1,tr2,an,fn1,fn2,v1,si,xc,dih,en,wei,emcu,cpm,rtan)
       DO i=1,ver(vernu)%nonei
       im1=i-1 ; IF(i.EQ.1) im1=ver(vernu)%nonei
       ip1=i+1 ; IF(i.EQ.ver(vernu)%nonei) ip1=1

       v2=ver(ver(vernu)%vneipt(i))%vcoord-ver(vernu)%vcoord

       elen=SQRT(SUM(v2**2))                                       ! Length of the edge from vernu --> i    
       v2=v2/elen ; euv=v2                                         ! Unit vector along the edge and along i 
       tr1=ver(vernu)%vneitr(im1) ; tr2=ver(vernu)%vneitr(i)       ! Two triangles sharing the edge 
       fn1=tri(tr1)%fnor ; fn2=tri(tr2)%fnor   
       an=MATMUL(TRANSPOSE(fn1),fn2)
       IF(an(1,1).GT.1.0000000) an(1,1)=1.0000000000000            ! To compute the angle between the faces sharing e 
       
       v1(1,1)=fn1(2,1)*fn2(3,1)-fn1(3,1)*fn2(2,1)                 ! N1 X N2 
       v1(2,1)=fn1(3,1)*fn2(1,1)-fn1(1,1)*fn2(3,1)
       v1(3,1)=fn1(1,1)*fn2(2,1)-fn1(2,1)*fn2(1,1)

       si=1.000 
       xc=MATMUL(TRANSPOSE(euv),v1)                                ! e.(N1XN2) 
       si=SIGN(si,xc(1,1))                                         ! Convex or concave
       dih=pi+si*acos(an(1,1))                                     ! Dihedral angle
        
       en=tri(tr1)%fnor+tri(tr2)%fnor                              ! Normal along the above edge    
       en=en/SQRT(SUM(en**2))                                      ! Normalized                     
       wei=MATMUL(TRANSPOSE(nor),en)
       emcu=2*elen*cos(dih*0.5)                                    ! Mean curvature along a the edge 
       cpm(1,1)=(euv(2,1)*en(3,1)-euv(3,1)*en(2,1))
       cpm(2,1)=(euv(3,1)*en(1,1)-euv(1,1)*en(3,1))                ! Direction orthogonal to edge and edge normal 
       cpm(3,1)=(euv(1,1)*en(2,1)-euv(2,1)*en(1,1))
       rtan=MATMUL(pmat,cpm) ; rtan=rtan/SQRT(SUM(rtan**2))
       cmat=cmat+0.5*wei(1,1)*emcu*MATMUL(rtan,TRANSPOSE(rtan))
       ENDDO
       !$OMP End Parallel Do

                amat=zdir+nor ; smat=zdir-nor  
                householder:IF(SUM(amat**2) .GT. SUM(smat**2))THEN       ! Use Householder trans to reduce to a 2X2 
                amat=amat/SQRT(SUM(amat**2))
                hmat=-(unitmat-2*MATMUL(amat,TRANSPOSE(amat)))
                ELSE
                smat=smat/SQRT(SUM(smat**2))
                hmat=(unitmat-2*MATMUL(smat,TRANSPOSE(smat)))
                ENDIF householder

                pmat=0
                pmat=MATMUL(TRANSPOSE(hmat),MATMUL(cmat,hmat))           ! Diagonalize the constructed matrix 

                p=pmat(1,1) ; q=pmat(1,2)                                ! Components of the 2X2 minor 
                r=pmat(2,1) ; s=pmat(2,2) 

                non_diagonal:IF(q.NE.0.0 .AND. r.NE.0.0)THEN             ! Non diagonal matrices eigen values 
                insq=(p+s)**2-4*(p*s-q*r)

                   IF((insq .GT. 0.0))THEN                               ! Complex values are avoided                        
                      e1=((p+s)+SQRT(insq))*0.5                          ! Largest eigenvalue is less -ive or more +ive      
                      e2=((p+s)-SQRT(insq))*0.5                          ! The smallest eigenvalue is more -ive or less +ive 
                   ELSE
                     e1=(p+s)*0.5                                        ! Degenerate eigenvalues 
                     e2=e1
                   ENDIF

                  evec(1,1)=q/SQRT(q**2+(e1-p)**2)
                  evec(2,1)=(e1-p)/SQRT(q**2+(e1-p)**2)
                  evec(1,2)=-evec(2,1)
                  evec(2,2)=evec(1,1)
                  evec(3,3)=1 
              ELSE
                 IF(p.GT.s)THEN                                          ! Picking up the largest eigenvalue for diagonal matrix 
                 e1=p ; e2=s
                 evec=unitmat                                            ! The eigenvector is same as the unit matrix 
                 ELSE
                 e1=s ; e2=p
                 evec(:,1:1)=unitmat(:,2:2) 
                 evec(:,2:2)=unitmat(:,1:1) 
                 evec(:,3:3)=unitmat(:,3:3) 
                 ENDIF
               ENDIF non_diagonal

              IF(abs(e1).LT.10E-10) e1=0
              IF(abs(e2).LT.10E-10) e2=0
              
              ver(vernu)%cur1=e1/ver(vernu)%totarea                      ! Principal curvature 1                       
              ver(vernu)%cur2=e2/ver(vernu)%totarea                      ! Principal curvature 2                       
              ver(vernu)%mcur=(ver(vernu)%cur1+ver(vernu)%cur2)          ! Mean curvature (theorema egregium of Gauss) 

              g2l=MATMUL(hmat,evec)                                      ! Calculate the Local to global trans matrix  
              ver(vernu)%L2G=g2l   
              ver(vernu)%spgl=MATMUL(g2l,ver(vernu)%splo)                ! Equivalent to doing MATMUL(hmat,evec(:,1:1 or 2:2))
              ver(vernu)%t1(1:3,1:1)=Matmul(hmat,evec(1:3,1:1))
              ver(vernu)%t2(1:3,1:1)=Matmul(hmat,evec(1:3,2:2))           ! eigenvectors
        END SUBROUTINE normalcalc
        END MODULE module_curvcalc   

