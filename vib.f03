INCLUDE "vib_mod.f03"
      Program vib
!
!     This program carries out Harmonic Vibrational analysis.
!
!     -A. J. Bovill, 2023
!
!
!     USE Connections
!
      use vib_mod
!
!     Variable Declarations
!
      integer(kind=int64)::nCommands,iPrint,i,nAtoms, &
        nBasis,nBasisUse,nEl,nElAlpha,nElBeta,l 
      integer(kind=int64)::j,k,n,info,nAt3,iAtom,jAtom
      integer(kind=int64),allocatable,dimension(:)::tmpVectorInt
      real(kind=real64),dimension(:,:),allocatable::hMat,hEVecs,hMatMW,  &
        hMatProjectorVectors,hMatProjector,RotVecs,tmpMat1,tmpMat2,  &
        hEVals,tmpVec
      real(kind=real64),allocatable,dimension(:)::hess_eigenvalues,hess_pack, &
        scratch_array,atomicMasses
      real(kind=real64),allocatable,dimension(:,:)::tmpMatrix1,  &
        tmpMatrix2,hess_eigenvectors
      character(len=512)::matrixFilename1,  &
        matrixFilenameOut
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile
      type(mqc_variable)::hess,grad
      type(mqc_vector)::mqc_hess_eigenvalues,mqc_vector1
      type(mqc_matrix)::mqc_hess_eigenvectors,hess_mw
      real(kind=real64),parameter::conv_factor=48169.11381488435
      real(kind=real64)::tmp
!
!     Format Statements
!

 1000 Format(1x,'Enter Program vib.')
 1010 Format(1x,'Matrix File:      ',A,/)  
 1100 Format(1x,'nAtoms=',I4,3x,'nBasis   =',I4,3x,'nBasisUse=',I4,/,  &
             1x,'nEl   =',I4,3x,'nElAlpha =',I4,3x,'nElBeta  =',I4,/)
      write(IOut,1000)

 8999 format(/,1x,'end of vib program')

!
!     Test if working version of MQC
!
      call mqc_version_print(iOut)
      call get_command_argument(1,matrixFilename1)
      call GMatrixFile%load(matrixFilename1)

      write(IOut,1010) TRIM(matrixFilename1)

!     Do some consistency checks and load the number of atoms, basis functions,
!     and linearly independent basis functions.
!
      nAtoms  = GMatrixFile%getVal('nAtoms')
      nBasis  = GMatrixFile%getVal('nBasis')
      nBasisUse  = GMatrixFile%getVal('nBasisUse')
      nEl1      = GMatrixFile%getVal('nElectrons')
      nElAlpha1 = GMatrixFile%getVal('nAlpha')
      nElBeta1  = GMatrixFile%getVal('nBeta')
!
!     Print out information from the matrix file
!
      write(IOut,1100) nAtoms,nBasis,nBasisUse,nEl1,nElAlpha1,nElBeta1

!
!     Get Hessian and nuclear gradients
!

      call GMatrixFile%getArray('NUCLEAR FORCE CONSTANTS',mqcVarOut=hess)
      call hess%print(header='Hessian Matrix')
      atomicMasses = GMatrixFile%getatomweights()

!
!     Diagonalize Hessian matrix to get eigenvalues and eigenvectors.
!     Before we diagonalize Hessian we massweigh and then put in 
!     compact vector form.
!

      tmpMatrix1 = hess
      n = size(tmpMatrix1,1)
      
      allocate(hess_pack((n*(n+1))/2))
      allocate(hess_eigenvalues(n))
      allocate(hess_eigenvectors(n,n))
      allocate(scratch_array(3*n))
      allocate(tmpMatrix2(n,n))
      allocate(tmpMat2(n,n))

      tmpMatrix2 = 0.0

      Nat3= 3*nAtoms 
      do i = 1,nAt3
        do j = 1,nAt3
         iAtom = int((i-1)/3)+1 
         jAtom = int((j-1)/3)+1 
         tmpMatrix2(i,j) = tmpMatrix1(i,j)/(sqrt(atomicMasses(iAtom)*atomicMasses(jAtom)))
        end do
      end do



      hess_mw = tmpMatrix2
      tmpMatrix2 = 0.0
      !Andrew do scaling
      call hess_mw%print(IOut,header='Hessian matrix massweight')
      !hess_mw = hess_mw*548.5799089940967d0
      call hess_mw%print(IOut,header='Hessian matrix massweight scaled')

      do i = 1,n
        do j = 1,n
          tmp = hess_mw%at(i,j)
          tmpMat2(i,j) = tmp 
        end do
      end do

      k = 0
      do j = 1, n
        do i = 1, j
          k = k + 1
          hess_pack(k) = tmpMat2(j, i)
        end do
      end do 

      write(*,*) "Size of array: ", size(hess_pack)
      mqc_vector1 = hess_pack
      call mqc_vector1%print(IOut,header='Hessian vector massweight packed')

!
!     Call lapack dspev 
!

      call hess_mw%diag(evals=mqc_hess_eigenvalues,evecs=mqc_hess_eigenvectors)

!      call dspev('V', 'L', n, hess_pack,hess_eigenvalues, &
!        hess_eigenvectors,n,scratch_array,info)

!     mqc_hess_eigenvectors = hess_eigenvectors
!     mqc_hess_eigenvalues = hess_eigenvalues

      call mqc_hess_eigenvectors%print(IOut,header='Hessian eigenvectors')
      call mqc_hess_eigenvalues%print(IOut,header='Hessian eigenvalues')

      mqc_hess_eigenvalues = hess_eigenvalues*conv_factor

      call mqc_hess_eigenvalues%print(IOut,header='Hessian eigenvalues after &
        conversion factor')

  999 Continue
      write(iOut,8999)
!      call GMatrixfile1%MQC_Update_Mol_Data(Molecule_Info,label,scalar,vector,matrix)

      End Program vib

