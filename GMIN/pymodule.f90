
 MODULE PYMODULE

 INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, K1, K2, OFFSET, REALNATOMS
 DOUBLE PRECISION, ALLOCATABLE :: RMIvec(:,:,:), DPI1RMvec(:,:,:), DPI2RMvec(:,:,:), DPI3RMvec(:,:,:)
 DOUBLE PRECISION, ALLOCATABLE :: PSCALEFAC1vec(:),PSCALEFAC2vec(:),epsilon1(:,:,:),AEZR1(:,:,:), AEZR2(:,:,:)


 DOUBLE PRECISION :: angle,angle2,pi,twopi,sigma1(4),cut,ron,ron2,range2inv3,attr(4),vecsbf(3)

 DOUBLE PRECISION :: I3(3,3)


!
! sf344> multisite PY additions
!
       DOUBLE PRECISION, ALLOCATABLE :: PST(:,:),OST(:,:),ELLST1(:,:),ELLST2(:,:),ELLMAT(:,:,:),SITECOORDS(:,:)


END MODULE PYMODULE
 
