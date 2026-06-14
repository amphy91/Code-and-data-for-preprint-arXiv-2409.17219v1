module constants
    implicit none
    INTEGER,PARAMETER::NQ=40 ! # of q-points in each high symmetry paths along which S(q,\omega) has been calculated
    integer,parameter::dim=2 ! Lattice dimension
    INTEGER,PARAMETER::sub_lat=6 ! # of sites per unit cell
    integer,parameter::L=10 ! system size L\timesL
    integer,parameter::N_SITE=SUB_LAT*L**DIM !Total # of sites
    INTEGER,PARAMETER::PP=N_SITE
    INTEGER,PARAMETER::NS=N_SITE
    INTEGER,parameter::ASF(12)=(/ 1,0,0,1,0,1,0,0,0,0,0,1 /)
    !ASF is an array that specifies which singlet mean-field amplitudes are 
    !allowed in the order \tau^0,\tau^z,\tau^x,\tau^y. 0 and 1 correspond
    !to presence or absence of the mean field amplitudes. The first four, second four
    !and third four elements corresponds to the hexagonal (2NN), tringular and 3NN bonds respectively.
    !The given values correspond to the 8th row of Table I.
    INTEGER,parameter::ATF(12)=(/ 0,1,0,0,1,0,0,1,0,0,0,0 /)
    !ATF has also been defined in the same manner as ASF but for triplet mean-field amplitudes
    REAL(8),parameter::PI=4.0D0*ATAN(1.0D0)  
    COMPLEX(8),PARAMETER::IC=DCMPLX(0.0D0,1.0D0) ! definition of i, i.e, \sqrt{-1} 
    !In the following we define different sign parameters for defining ansatze in table I
    INTEGER,PARAMETER::ET=1 !\eta
    INTEGER,PARAMETER::ET_R=-ET !\eta_{C_6}
    INTEGER,PARAMETER::ET_RS=-1 !\eta_{C_6R}
    INTEGER,PARAMETER::ET_SIG=-1 !\eta_R
    INTEGER,PARAMETER::NW=200 ! # \omega points
    INTEGER,PARAMETER::N_HSP=4 ! # of high symmetry points along the high symmetry paths (\Gamma->M'->K'->\Gamma)
    REAL(8),PARAMETER::deltaNew=0.35D0 ! fixed magnitude of onsite triplet hopping term
end module constants 



 
 
PROGRAM MAIN
  use constants
  IMPLICIT NONE
  INTEGER::I,J,K,M1,M2,NQQ,IR1
  REAL(8)::POS_SITE(NS,DIM)
  COMPLEX(8)::U(4*NS,4*NS)
  INTEGER::SITE_IJS_N(L,L,SUB_LAT),SITE_N_IJS(NS,4)
  REAL(8)::EV(4*NS),H1,H2,H3,qx,qy,q(DIM),szz,HSP(N_HSP,DIM)
  REAL(8)::WH,WL,DW,SIGMA,MF(3),SF(12),TF(12), RATIO(6,3)

  open(40,file='dsf8a.dat',status='unknown')
  
  RATIO(1,:)=(/1.25D0, 1.0D0, 0.0D0/) ! Ratios of mean-field parameters 2NN:1NN:3NN


  DO I=1,3
  MF(I)=RATIO(1,I)/SQRT(1.0D0*SUM(ASF(4*(I-1)+1:4*I))+1.0D0*SUM(ATF(4*(I-1)+1:4*I)))
  END DO


  DO I=1,3
    DO J=1,4
      SF(4*(I-1)+J)=ASF(4*(I-1)+J)*MF(I)
      TF(4*(I-1)+J)=ATF(4*(I-1)+J)*MF(I)
    END DO
  END DO

  TF(1:4)=2.0D0*TF(1:4)
  TF(5:8)=2.0D0*TF(5:8)

  CALL LAT(SF,TF,POS_SITE,SITE_IJS_N,SITE_N_IJS,EV,U)

  WH=EV(1)*2.1 ! Maximum values of \omega has been seit 2.1 times the width of the dispersion energies.
  WL=0.0D0  
  DW=(WH-WL)/NW     
  SIGMA=2.0D0*DW ! delta function has been approxiamted as Gaussian function with standard deviation sigma

  HSP(1,:)=(/0.0D0,0.0D0/) ! \Gamma high symmetry points
  HSP(2,:)=(/2.0D0*PI,0.0D0/) ! M' high symmetry points
  HSP(3,:)=(/2.0D0*PI,2.0D0*PI/SQRT(3.0D0)/) ! K' high symmetry points
  HSP(4,:)=(/0.0D0,0.0D0/) ! \Gamma high symmetry points
  M1=0
  DO I=1,N_HSP-1
  H1=(HSP(I+1,1)-HSP(I,1))/NQ
  H2=(HSP(I+1,2)-HSP(I,2))/NQ
  Q=HSP(I,:)
  IF(I==N_HSP-1) THEN
  NQQ=NQ+1
  ELSE
  NQQ=NQ
  END IF
  DO J=1,NQQ
  M1=M1+1
  CALL DSF(WH,WL,DW,SIGMA,M1,POS_SITE,Q,EV,U)
  Q=Q+(/H1,H2/)
  END DO
  END DO
  
  close(40)
END PROGRAM MAIN

! Construct Hamiltonian and calculate eigenstates and eigenvalues
SUBROUTINE LAT(SF,TF,POS_SITE,SITE_IJS_N,SITE_N_IJS,EV,UV)
USE constants
IMPLICIT NONE
INTEGER,INTENT(OUT)::SITE_IJS_N(L,L,SUB_LAT),SITE_N_IJS(NS,3) 
REAL(8),INTENT(INOUT)::SF(12),TF(12)
REAL(8),INTENT(OUT)::POS_SITE(NS,2),EV(4*NS)
COMPLEX(8),INTENT(OUT)::UV(4*NS,4*NS)
INTEGER::I,J,K,M,U_SITE,D_SITE_IJS_N(0:L+1,0:L+1,SUB_LAT)  
REAL(8)::T1(2),T2(2),SL(SUB_LAT,2),LAT_PARA,DELTA
COMPLEX(8)::HAM(4*NS,4*NS),A(4*NS,4*NS)
COMPLEX(8)::SH(NS,NS),TH(NS,NS),SP(NS,NS),TP(NS,NS)
D_SITE_IJS_N=0
LAT_PARA=1.0D0
T1=(SQRT(3.0D0)+1.0D0/SQRT(3.0D0))*(/SQRT(3.0D0),-1.0D0/)/2.0D0
T2=(/0.0D0,SQRT(3.0D0)+1.0D0/SQRT(3.0D0)/)
SL(1,:)=(/-0.50D0,Sqrt(3.0D0)/2.0D0/)
SL(2,:)=(/-1.0D0,0.0D0/)
SL(3,:)=(/-0.50D0,-Sqrt(3.0D0)/2.0D0/)
SL(4,:)=(/0.50D0,-Sqrt(3.0D0)/2.0D0/)
SL(5,:)=(/1.0D0,0.0D0/)
SL(6,:)=(/0.50D0,Sqrt(3.0D0)/2.0D0/)
M=0
DO I=1,L
DO J=1,L
 DO U_SITE=1,SUB_LAT
    M=M+1
    SITE_IJS_N(I,J,U_SITE)=M
    SITE_N_IJS(M,:)=(/ I, J, U_SITE /)
    POS_SITE(M,:)=I*T1(:)+J*T2(:)+SL(U_SITE,:)
    D_SITE_IJS_N(I,J,U_SITE)=M
 END DO
END DO
END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!PERIODICITY!!!!!!!!!!!!!!!!!!!!!!!!!!
D_SITE_IJS_N(:,0,:)=D_SITE_IJS_N(:,L,:)
D_SITE_IJS_N(:,L+1,:)=D_SITE_IJS_N(:,1,:)
D_SITE_IJS_N(0,:,:)=D_SITE_IJS_N(L,:,:)
D_SITE_IJS_N(L+1,:,:)=D_SITE_IJS_N(1,:,:)
!!!!!!!!!!!!!!!!PERIODICITY!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!HAMILTONIAN!!!!!!!!!!!!!!!!!!!!!!!!!!!
DELTA=deltaNew
A=0.0D0
SH=0.0D0 !contains Hamiltonin matrix elements for siglet hopping
TH=0.0D0 !contains Hamiltonin matrix elements for triplet hopping
SP=0.0D0 !contains Hamiltonin matrix elements for siglet pairing
TP=0.0D0 !contains Hamiltonin matrix elements for triplet pairing
DO I=1,L
DO J=1,L
    SH(D_SITE_IJS_N(I,J,1),D_SITE_IJS_N(I,J,2))=(IC*SF(1)+SF(2))
    SH(D_SITE_IJS_N(I,J,2),D_SITE_IJS_N(I,J,3))=(IC*SF(1)+SF(2))
    SH(D_SITE_IJS_N(I,J,3),D_SITE_IJS_N(I,J,4))=(IC*SF(1)+SF(2))
    SH(D_SITE_IJS_N(I,J,4),D_SITE_IJS_N(I,J,5))=(IC*SF(1)+SF(2))
    SH(D_SITE_IJS_N(I,J,5),D_SITE_IJS_N(I,J,6))=(IC*SF(1)+SF(2))*ET_R
    SH(D_SITE_IJS_N(I,J,6),D_SITE_IJS_N(I,J,1))=(IC*SF(1)+SF(2))
    TH(D_SITE_IJS_N(I,J,1),D_SITE_IJS_N(I,J,2))=(TF(1)+IC*TF(2))
    TH(D_SITE_IJS_N(I,J,2),D_SITE_IJS_N(I,J,3))=(TF(1)+IC*TF(2))
    TH(D_SITE_IJS_N(I,J,3),D_SITE_IJS_N(I,J,4))=(TF(1)+IC*TF(2))
    TH(D_SITE_IJS_N(I,J,4),D_SITE_IJS_N(I,J,5))=(TF(1)+IC*TF(2))
    TH(D_SITE_IJS_N(I,J,5),D_SITE_IJS_N(I,J,6))=(TF(1)+IC*TF(2))*ET_R
    TH(D_SITE_IJS_N(I,J,6),D_SITE_IJS_N(I,J,1))=(TF(1)+IC*TF(2))

    SP(D_SITE_IJS_N(I,J,1),D_SITE_IJS_N(I,J,2))=(SF(3)-IC*SF(4))
    SP(D_SITE_IJS_N(I,J,2),D_SITE_IJS_N(I,J,3))=(SF(3)-IC*SF(4))
    SP(D_SITE_IJS_N(I,J,3),D_SITE_IJS_N(I,J,4))=(SF(3)-IC*SF(4))
    SP(D_SITE_IJS_N(I,J,4),D_SITE_IJS_N(I,J,5))=(SF(3)-IC*SF(4))
    SP(D_SITE_IJS_N(I,J,5),D_SITE_IJS_N(I,J,6))=(SF(3)-IC*SF(4))*ET_R
    SP(D_SITE_IJS_N(I,J,6),D_SITE_IJS_N(I,J,1))=(SF(3)-IC*SF(4))
    TP(D_SITE_IJS_N(I,J,1),D_SITE_IJS_N(I,J,2))=(TF(3)-IC*TF(4))
    TP(D_SITE_IJS_N(I,J,2),D_SITE_IJS_N(I,J,3))=(TF(3)-IC*TF(4))
    TP(D_SITE_IJS_N(I,J,3),D_SITE_IJS_N(I,J,4))=(TF(3)-IC*TF(4))
    TP(D_SITE_IJS_N(I,J,4),D_SITE_IJS_N(I,J,5))=(TF(3)-IC*TF(4))
    TP(D_SITE_IJS_N(I,J,5),D_SITE_IJS_N(I,J,6))=(TF(3)-IC*TF(4))*ET_R
    TP(D_SITE_IJS_N(I,J,6),D_SITE_IJS_N(I,J,1))=(TF(3)-IC*TF(4))


    SH(D_SITE_IJS_N(I,J,1),D_SITE_IJS_N(I,J+1,3))=(IC*SF(5)+SF(6))*(ET)**(I+1)
    SH(D_SITE_IJS_N(I,J,2),D_SITE_IJS_N(I-1,J,4))=(IC*SF(5)+SF(6))
    SH(D_SITE_IJS_N(I,J,3),D_SITE_IJS_N(I-1,J-1,5))=(IC*SF(5)+SF(6))*ET*(ET)**(I+1)
    SH(D_SITE_IJS_N(I,J,4),D_SITE_IJS_N(I,J-1,6))=(IC*SF(5)+SF(6))*ET*ET_R*(ET)**(I+1)
    SH(D_SITE_IJS_N(I,J,5),D_SITE_IJS_N(I+1,J,1))=(IC*SF(5)+SF(6))*ET*ET_R
    SH(D_SITE_IJS_N(I,J,6),D_SITE_IJS_N(I+1,J+1,2))=(IC*SF(5)+SF(6))*ET*(ET)**(I+1)
    TH(D_SITE_IJS_N(I,J,1),D_SITE_IJS_N(I,J+1,3))=(TF(5)+IC*TF(6))*(ET)**(I+1)
    TH(D_SITE_IJS_N(I,J,2),D_SITE_IJS_N(I-1,J,4))=(TF(5)+IC*TF(6))
    TH(D_SITE_IJS_N(I,J,3),D_SITE_IJS_N(I-1,J-1,5))=(TF(5)+IC*TF(6))*ET*(ET)**(I+1)
    TH(D_SITE_IJS_N(I,J,4),D_SITE_IJS_N(I,J-1,6))=(TF(5)+IC*TF(6))*ET*ET_R*(ET)**(I+1)
    TH(D_SITE_IJS_N(I,J,5),D_SITE_IJS_N(I+1,J,1))=(TF(5)+IC*TF(6))*ET*ET_R
    TH(D_SITE_IJS_N(I,J,6),D_SITE_IJS_N(I+1,J+1,2))=(TF(5)+IC*TF(6))*ET*(ET)**(I+1)

    SP(D_SITE_IJS_N(I,J,1),D_SITE_IJS_N(I,J+1,3))=(SF(7)-IC*SF(8))*(ET)**(I+1)
    SP(D_SITE_IJS_N(I,J,2),D_SITE_IJS_N(I-1,J,4))=(SF(7)-IC*SF(8))
    SP(D_SITE_IJS_N(I,J,3),D_SITE_IJS_N(I-1,J-1,5))=(SF(7)-IC*SF(8))*ET*(ET)**(I+1)
    SP(D_SITE_IJS_N(I,J,4),D_SITE_IJS_N(I,J-1,6))=(SF(7)-IC*SF(8))*ET*ET_R*(ET)**(I+1)
    SP(D_SITE_IJS_N(I,J,5),D_SITE_IJS_N(I+1,J,1))=(SF(7)-IC*SF(8))*ET*ET_R
    SP(D_SITE_IJS_N(I,J,6),D_SITE_IJS_N(I+1,J+1,2))=(SF(7)-IC*SF(8))*ET*(ET)**(I+1)
    TP(D_SITE_IJS_N(I,J,1),D_SITE_IJS_N(I,J+1,3))=(TF(7)-IC*TF(8))*(ET)**(I+1)
    TP(D_SITE_IJS_N(I,J,2),D_SITE_IJS_N(I-1,J,4))=(TF(7)-IC*TF(8))
    TP(D_SITE_IJS_N(I,J,3),D_SITE_IJS_N(I-1,J-1,5))=(TF(7)-IC*TF(8))*ET*(ET)**(I+1)
    TP(D_SITE_IJS_N(I,J,4),D_SITE_IJS_N(I,J-1,6))=(TF(7)-IC*TF(8))*ET*ET_R*(ET)**(I+1)
    TP(D_SITE_IJS_N(I,J,5),D_SITE_IJS_N(I+1,J,1))=(TF(7)-IC*TF(8))*ET*ET_R
    TP(D_SITE_IJS_N(I,J,6),D_SITE_IJS_N(I+1,J+1,2))=(TF(7)-IC*TF(8))*ET*(ET)**(I+1)


    SH(D_SITE_IJS_N(I,J,1),D_SITE_IJS_N(I-1,J,4))=(IC*SF(9)+SF(10))
    SH(D_SITE_IJS_N(I,J,2),D_SITE_IJS_N(I-1,J,5))=(IC*SF(9)+SF(10))*(-ET*ET_R*ET_RS)
    SH(D_SITE_IJS_N(I,J,3),D_SITE_IJS_N(I,J-1,6))=(IC*SF(9)+SF(10))*ET*ET_R*(ET)**(I+1)
    SH(D_SITE_IJS_N(I,J,4),D_SITE_IJS_N(I,J-1,1))=(IC*SF(9)+SF(10))*(-ET_RS)*(ET)**(I+1)
    SH(D_SITE_IJS_N(I,J,5),D_SITE_IJS_N(I+1,J+1,2))=(IC*SF(9)+SF(10))*ET*ET_R*(ET)**(I+1)
    SH(D_SITE_IJS_N(I,J,6),D_SITE_IJS_N(I+1,J+1,3))=(IC*SF(9)+SF(10))*(-ET_R*ET_RS)*(ET)**(I+1)
    TH(D_SITE_IJS_N(I,J,1),D_SITE_IJS_N(I-1,J,4))=(TF(9)+IC*TF(10))
    TH(D_SITE_IJS_N(I,J,2),D_SITE_IJS_N(I-1,J,5))=(TF(9)+IC*TF(10))*(ET*ET_R*ET_RS)
    TH(D_SITE_IJS_N(I,J,3),D_SITE_IJS_N(I,J-1,6))=(TF(9)+IC*TF(10))*ET*ET_R*(ET)**(I+1)
    TH(D_SITE_IJS_N(I,J,4),D_SITE_IJS_N(I,J-1,1))=(TF(9)+IC*TF(10))*(ET_RS)*(ET)**(I+1)
    TH(D_SITE_IJS_N(I,J,5),D_SITE_IJS_N(I+1,J+1,2))=(TF(9)+IC*TF(10))*ET*ET_R*(ET)**(I+1)
    TH(D_SITE_IJS_N(I,J,6),D_SITE_IJS_N(I+1,J+1,3))=(TF(9)+IC*TF(10))*(ET_R*ET_RS)*(ET)**(I+1)   
    
    SP(D_SITE_IJS_N(I,J,1),D_SITE_IJS_N(I-1,J,4))=(IC*SF(11)-SF(12))
    SP(D_SITE_IJS_N(I,J,2),D_SITE_IJS_N(I-1,J,5))=(IC*SF(11)-SF(12))*(ET_SIG)*(-ET*ET_R*ET_RS)
    SP(D_SITE_IJS_N(I,J,3),D_SITE_IJS_N(I,J-1,6))=(IC*SF(11)-SF(12))*ET*ET_R*(ET)**(I+1)
    SP(D_SITE_IJS_N(I,J,4),D_SITE_IJS_N(I,J-1,1))=(IC*SF(11)-SF(12))*(ET_SIG)*(-ET_RS)*(ET)**(I+1)
    SP(D_SITE_IJS_N(I,J,5),D_SITE_IJS_N(I+1,J+1,2))=(IC*SF(11)-SF(12))*ET*ET_R*(ET)**(I+1)
    SP(D_SITE_IJS_N(I,J,6),D_SITE_IJS_N(I+1,J+1,3))=(IC*SF(11)-SF(12))*(ET_SIG)*(-ET_R*ET_RS)*(ET)**(I+1)
    TP(D_SITE_IJS_N(I,J,1),D_SITE_IJS_N(I-1,J,4))=(TF(11)-IC*TF(12))
    TP(D_SITE_IJS_N(I,J,2),D_SITE_IJS_N(I-1,J,5))=(TF(11)-IC*TF(12))*(ET_SIG)*(ET*ET_R*ET_RS)
    TP(D_SITE_IJS_N(I,J,3),D_SITE_IJS_N(I,J-1,6))=(TF(11)-IC*TF(12))*ET*ET_R*(ET)**(I+1)
    TP(D_SITE_IJS_N(I,J,4),D_SITE_IJS_N(I,J-1,1))=(TF(11)-IC*TF(12))*(ET_SIG)*(ET_RS)*(ET)**(I+1)
    TP(D_SITE_IJS_N(I,J,5),D_SITE_IJS_N(I+1,J+1,2))=(TF(11)-IC*TF(12))*ET*ET_R*(ET)**(I+1)
    TP(D_SITE_IJS_N(I,J,6),D_SITE_IJS_N(I+1,J+1,3))=(TF(11)-IC*TF(12))*(ET_SIG)*(ET_R*ET_RS)*(ET)**(I+1)
END DO
END DO
A(1:NS,1:NS)=SH+TH
A(NS+1:2*NS,NS+1:2*NS)=SH-TH
A(2*NS+1:3*NS,2*NS+1:3*NS)=-CONJG(SH)-CONJG(TH)
A(3*NS+1:4*NS,3*NS+1:4*NS)=-CONJG(SH)+CONJG(TH)

A(1:NS,3*NS+1:4*NS)=SP+IC*TP
A(NS+1:2*NS,2*NS+1:3*NS)=-SP+IC*TP
A(2*NS+1:3*NS,NS+1:2*NS)=-CONJG(SP)+IC*CONJG(TP)
A(3*NS+1:4*NS,1:NS)=CONJG(SP)+IC*CONJG(TP)

DO I=1,NS
  A(I,I)=DELTA
  A(I+NS,I+NS)=-DELTA
  A(I+2*NS,I+2*NS)=-DELTA
  A(I+3*NS,I+3*NS)=DELTA
END DO

!!!!!!!!!!!!!!!HAMILTONIAN!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
HAM=A
HAM=HAM+CONJG(TRANSPOSE(HAM))

CALL EIGEN(HAM,EV,UV)          
END SUBROUTINE LAT




SUBROUTINE EIGEN(HAM,EV,UV)
USE constants
IMPLICIT NONE
COMPLEX(8),INTENT(IN)::HAM(4*NS,4*NS)
REAL(8),INTENT(OUT)::EV(4*NS)
COMPLEX(8),INTENT(OUT)::UV(4*NS,4*NS)
integer::I,J,K,LWORK,INFO,LDA
real(8)::W(4*NS),RWORK(3*(4*NS)-2)
COMPLEX(8)::WORK(2*(4*NS)-1)
CALL ZHEEV('V','U',4*NS,HAM,4*NS,W,WORK,2*(4*NS)-1,RWORK,INFO)                
EV(2*NS+1:4*NS)=W(1:2*NS)
UV(1:4*NS,2*NS+1:4*NS)=HAM(1:4*NS,1:2*NS)  
DO I=1,2*NS
  EV(I)=W(4*NS-I+1)
  UV(1:4*NS,I)=HAM(1:4*NS,4*NS-I+1)
END DO      
EV=ABS(EV)
!print*,sum(uv(:,1))        
END SUBROUTINE EIGEN

SUBROUTINE DSF(WH,WL,DW,SIGMA1,M2,POS_SITE,Q,EV,U)
  USE constants
  IMPLICIT NONE
  INTEGER,INTENT(IN)::M2
  REAL(8),INTENT(IN)::EV(4*N_SITE),POS_SITE(N_SITE,DIM),Q(DIM)
  REAL(8),INTENT(IN)::WH,WL,DW,SIGMA1
  COMPLEX(8),INTENT(IN)::U(4*N_SITE,4*N_SITE)
  integer::I,J,MU,NU,M1
  COMPLEX(8)::STR1,STR,MAT_MUL(8,N_SITE,N_SITE),YMAT(8,N_SITE,N_SITE),UR(4*N_SITE,4*N_SITE)
  REAL(8)::W,SZZ,SIGMA

  DO M1=1,N_SITE
    UR(M1,:)=U(M1,:)*EXP(IC*DOT_PRODUCT(Q,POS_SITE(M1,:)))
    UR(M1+N_SITE,:)=U(M1,:)*EXP(IC*DOT_PRODUCT(Q,POS_SITE(M1,:)))
    UR(M1+2*N_SITE,:)=U(M1,:)*EXP(IC*DOT_PRODUCT(Q,POS_SITE(M1,:)))
    UR(M1+3*N_SITE,:)=U(M1,:)*EXP(IC*DOT_PRODUCT(Q,POS_SITE(M1,:)))
  END DO
  MAT_MUL(1,:,:)=MATMUL(CONJG(TRANSPOSE(UR(1:NS,2*NS+1:3*NS))),U(1:NS,1:NS))
  MAT_MUL(2,:,:)=MATMUL(CONJG(TRANSPOSE(UR(1:NS,3*NS+1:4*NS))),U(1:NS,NS+1:2*NS))
  MAT_MUL(3,:,:)=MATMUL(CONJG(TRANSPOSE(UR(1:NS,3*NS+1:4*NS))),U(1:NS,1:NS))
  MAT_MUL(4,:,:)=MATMUL(CONJG(TRANSPOSE(UR(1:NS,2*NS+1:3*NS))),U(1:NS,NS+1:2*NS))

  MAT_MUL(5,:,:)=MATMUL(CONJG(TRANSPOSE(UR(NS+1:2*NS,2*NS+1:3*NS))),U(NS+1:2*NS,1:NS))
  MAT_MUL(6,:,:)=MATMUL(CONJG(TRANSPOSE(UR(NS+1:2*NS,3*NS+1:4*NS))),U(NS+1:2*NS,NS+1:2*NS))
  MAT_MUL(7,:,:)=MATMUL(CONJG(TRANSPOSE(UR(NS+1:2*NS,3*NS+1:4*NS))),U(NS+1:2*NS,1:NS))
  MAT_MUL(8,:,:)=MATMUL(CONJG(TRANSPOSE(UR(NS+1:2*NS,2*NS+1:3*NS))),U(NS+1:2*NS,NS+1:2*NS))

  YMAT(1,:,:)=MAT_MUL(1,:,:)-TRANSPOSE( MAT_MUL(1,:,:))-(MAT_MUL(5,:,:)-TRANSPOSE( MAT_MUL(5,:,:)))
  YMAT(2,:,:)=MAT_MUL(2,:,:)-TRANSPOSE( MAT_MUL(2,:,:))-(MAT_MUL(6,:,:)-TRANSPOSE( MAT_MUL(6,:,:)))
  YMAT(3,:,:)=MAT_MUL(3,:,:)-TRANSPOSE( MAT_MUL(4,:,:))-(MAT_MUL(7,:,:)-TRANSPOSE( MAT_MUL(8,:,:)))
  YMAT(4,:,:)=MAT_MUL(4,:,:)-TRANSPOSE( MAT_MUL(3,:,:))-(MAT_MUL(8,:,:)-TRANSPOSE( MAT_MUL(7,:,:)))

  YMAT(5,:,:)=CONJG(TRANSPOSE(MAT_MUL(1,:,:)-MAT_MUL(5,:,:)))
  YMAT(6,:,:)=CONJG(TRANSPOSE(MAT_MUL(2,:,:)-MAT_MUL(6,:,:)))
  YMAT(7,:,:)=CONJG(TRANSPOSE(MAT_MUL(3,:,:)-MAT_MUL(7,:,:)))
  YMAT(8,:,:)=CONJG(TRANSPOSE(MAT_MUL(4,:,:)-MAT_MUL(8,:,:)))

  W=0.0D0
  DO I=1,NW+1
  
    SIGMA=SIGMA1!0.010d0 ! delta function has been approxiamted as Gaussian function with standard deviation sigma

  STR=DCMPLX(0.0D0,0.0D0)
  DO NU=1,NS
  DO MU=1,NS
  STR=STR+&
  YMAT(1,NU,MU)*YMAT(5,MU,NU)*((1.0D0/(SQRT(PI)*SIGMA))*EXP(-(W-EV(MU)-EV(NU+2*NS))**2/SIGMA**2))+&
  YMAT(2,NU,MU)*YMAT(6,MU,NU)*((1.0D0/(SQRT(PI)*SIGMA))*EXP(-(W-EV(MU+NS)-EV(NU+3*NS))**2/SIGMA**2))+&
  YMAT(3,NU,MU)*YMAT(7,MU,NU)*((1.0D0/(SQRT(PI)*SIGMA))*EXP(-(W-EV(MU)-EV(NU+3*NS))**2/SIGMA**2))+&
  YMAT(4,NU,MU)*YMAT(8,MU,NU)*((1.0D0/(SQRT(PI)*SIGMA))*EXP(-(W-EV(MU+NS)-EV(NU+2*NS))**2/SIGMA**2))
  END DO
  END DO
  SZZ=REAL(STR)/(4.0D0*NS)
  write(40,'(3F15.8)')M2*1.0D0,W,SZZ
  W=W+DW
  END DO
END SUBROUTINE DSF
