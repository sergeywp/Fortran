	module ModAndDiv

		contains

		!Простое умножение векторов в трехмерный массив
		subroutine MullNorm(A,B,C,D,L,M,N)
			integer, value :: L,M,N
			real*8 :: A(L), B(M), C(N), D(L,M,N)
			do k=1,N
				do j=1,M
					do i=1,L
						D(i,j,k)=A(i)*B(j)*C(k)
					enddo
				enddo
			enddo
		end subroutine MullNorm

		!Заполнение массивов одинаковыми значениями
		subroutine FillMassSim(A,B,C,L,M,N)
			integer, value :: L,M,N
			real*8 :: A(L), B(M), C(N)
			A=2.0
			B=3.0
			C=4.0
		end subroutine FillMassSim
		
		!Заполнение массивов разными значениями
		subroutine FillMassNotSim(A,B,C,L,M,N)
			integer, value :: L,M,N
			real*8 :: A(L), B(M), C(N)	
			do i=1,L
				A(i)=2.0*i
			enddo
			do i=1,M
				B(i)=99/(i+1)
			enddo
			do i=1,N
				C(i)=4.0*i/2.0
			enddo
		end subroutine FillMassNotSim

		!Вывод одномерного массива на экран
		subroutine PrintMass(A,L)
			integer, value :: L
			real*8 :: A(L)			
			do i=1,L
				write (*,*)A(i)
			enddo
		end subroutine PrintMass
		
		!Вывод трехмерного массива на экран
		subroutine PrintMass3(A,L,M,N)
			integer, value :: L,M,N
			real*8 :: A(L,M,N)
			do k=1,N
				do j=1,M
					do i=1,L
						write (*,*)A(i,j,k)
					enddo
				enddo
			enddo
		end subroutine PrintMass3

		!Развертывание одномрного массива в трехмерный
		!Входной массив in одномерный, выходной out трехмерный
		subroutine D1dTo3d(in,out,xdim,ydim,zdim)
			integer, value :: xdim,ydim,zdim
			real*8 ::in(xdim*ydim*zdim),out(xdim,ydim,zdim)
			do k=1,zdim
				do j=1,ydim
					do i=1,xdim
						!Здесь надо проверить правильность задания индексов
						out(i,j,k)=in(k+(j-1)*zdim+(i-1)*zdim*ydim)
					enddo
				enddo
			enddo
		end subroutine D1dTo3d

		!Сворачивание трехмерного массива в одномерный
		!Входной массив in трехмерны, выходной out одномерный
		subroutine D3dTo1d(in,out,xdim,ydim,zdim)
			integer, value :: xdim,ydim,zdim
			real*8 ::in(xdim,ydim,zdim),out(xdim*ydim*zdim)
			integer ind;
			ind=0
			do i=1,xdim
				do j=1,ydim
					do k=1,zdim
						!Здесь надо проверить правильность задания индексов
						ind=ind+1;
						out(ind)=in(i,j,k)						
					enddo
				enddo
			enddo
		end subroutine D3dTo1d

		subroutine GetIJK3D(A,AIJK,I,J,K,L,M,N)
			integer, value :: L,M,N,I,J,K
			real*8 :: A(L*M*N), AIJK
			AIJK=A(K+(J-1)*N+(I-1)*M*N)					
		end subroutine GetIJK3D

		subroutine MullOneThreeDimmMassOne(ResMass,MidMass,A,B,C,L,M,N)
			integer, value :: L,M,N
			real*8 :: A(L),B(M),C(N),ResMass(L*M*N),MidMass(L*M*N)
			do i=1,L*M*N
				intIndA=INT((i-1)/(N*M))+1
				intIndB=MOD(INT((i-1)/N),M)+1
				intIndC=MOD((i-1),N)+1
				ResMass(i)=MidMass(i)*A(intIndA)*B(intIndB)/C(intIndC)
			enddo			
		end subroutine MullOneThreeDimmMassOne

		subroutine MullOneThreeDimmMassThree(ResMass,MidMass,A,B,C,L,M,N)
			integer, value :: L,M,N
			real*8 :: A(L),B(M),C(N),ResMass(L,M,N),MidMass(L,M,N)
			do k=1,N
				do j=1,M
					do i=1,L
						ResMass(i,j,k)=MidMass(i,j,k)*A(i)*B(j)/C(k)
					enddo
				enddo
			enddo			
		end subroutine MullOneThreeDimmMassThree

		subroutine MullOneThreeDimmMassOneMinus(ResMass,MidMass,A,B,C,L,M,N)
			integer, value :: L,M,N
			real*8 :: A(L),B(M),C(N),ResMass(L*M*N),MidMass(L*M*N)
			do i=1,L*M*N
				intIndA=INT((i-1)/(N*M))+1    !стрХстолб
				intIndB=MOD(INT((i-1)/N),M)+1 !столбц
				intIndC=MOD((i-1),N)+1        !строка
				ResMass(i)=MidMass(i)*A(intIndA+1)*B(intIndB)/C(intIndC)
			enddo	
		end subroutine MullOneThreeDimmMassOneMinus

		subroutine MullOneThreeDimmMassThreeMinus(ResMass,MidMass,A,B,C,L,M,N)
			integer, value :: L,M,N
			real*8 :: A(L),B(M),C(N),ResMass(L,M,N),MidMass(L,M,N)
			do k=1,N
				do j=1,M
					do i=1,L
						ResMass(i,j,k)=MidMass(i,j,k)*A(i+1)*B(j)/C(k)
					enddo
				enddo
			enddo
		end subroutine MullOneThreeDimmMassThreeMinus

		subroutine MullOneThreeDimmMassOneCos(ResMass,MidMass,A,B,C,L,M,N)
			integer, value :: L,M,N
			real*8 :: A(L),B(M),C(N),ResMass(L*M*N),MidMass(L*M*N)
			do i=1,L*M*N
				intIndA=INT((i-1)/(N*M))+1    !стрХстолб i
				intIndB=MOD(INT((i-1)/N),M)+1 !столбц    j
				intIndC=MOD((i-1),N)+1        !строка    k
				if((intIndA<L).AND.(intIndB<M).AND.(intIndC<N)) ResMass(i)=3.*(DCOS(A(intIndA))-DCOS(A(intIndA+1)))*&
																			  (DCOS(B(intIndB))-DCOS(B(intIndB+1)))*&
																			  (DCOS(C(intIndC))-DCOS(C(intIndC+1)))
				if((intIndA>=L).OR.(intIndB>=M).OR.(intIndC>=N)) ResMass(i)=0.
			enddo	
		end subroutine MullOneThreeDimmMassOneCos

		subroutine MullOneThreeDimmMassThreeCos(ResMass,MidMass,A,B,C,L,M,N)
			integer, value :: L,M,N
			real*8 :: A(L),B(M),C(N),ResMass(L,M,N),MidMass(L,M,N)
			do k=1,N
				do j=1,M
					do i=1,L
						ResMass(i,j,k)=0.
					enddo
				enddo
			enddo
			do k=1,N-1
				do j=1,M-1
					do i=1,L-1
						ResMass(i,j,k)=3.*(DCOS(A(i))-DCOS(A(i+1)))*(DCOS(B(j))-DCOS(B(j+1)))*(DCOS(C(k))-DCOS(C(k+1)))
					enddo
				enddo
			enddo
		end subroutine MullOneThreeDimmMassThreeCos

		!Одномерное перемножение трехмерного массива
		subroutine MullOneDimm(A,B,C,D,L,M,N)
			integer, value :: L,M,N
			integer :: intIndA,intIndB,intIndC
			real*8 :: A(L), B(M), C(N), D(L,M,N), Buff(L*M*N)
			do i=1,L*N*M
				intIndA=INT((i-1)/(N*M))+1
				intIndB=MOD(INT((i-1)/N),M)+1
				intIndC=MOD((i-1),N)+1
				Buff(i)=A(intIndA)*B(intIndB)*C(intIndC)
			enddo
			!call PrintMass(Buff,L*M*N)
			call D1dTo3d(Buff,D,L,M,N)
		end subroutine MullOneDimm

	end module ModAndDiv

	program prog
	!-------------------------------
	! Интерфейс модуля:
	!-------------------------------
	! MullNorm(A,B,C,D,L,M,N)
	! FillMassSim(A,B,C,L,M,N)
	! FillMassNotSim(A,B,C,L,M,N)
	! PrintMass(A,L)
	! PrintMass3(A,L,M,N)
	! D1dTo3d(in,out,xdim,ydim,zdim)
	! D3dTo1d(in,out,xdim,ydim,zdim)
	! MullOneDimm(A,B,C,D,L,M,N)
	!-------------------------------

		use ModAndDiv
		implicit none
		!Примечание: В этих процедурах быстрее всего изменения идут по k, т е по внешнему циклу
		integer ::  i,j,k
		integer, parameter :: L=3, M=4, N=5
		real*8, allocatable :: D(:,:,:), A(:), B(:), C(:)
		real*8, allocatable :: Buff1(:),Buff3(:,:,:)

		real*8, allocatable :: ResMass3(:,:,:), MidMass3(:,:,:), ResMass1(:), MidMass1(:)

		allocate(D(L,M,N), A(L), B(M), C(N), Buff1(L*M*N), Buff3(L,M,N),&
						ResMass3(L,M,N), MidMass3(L,M,N),ResMass1(L*M*N), MidMass1(L*M*N) )


		call FillMassNotSim(A,B,C,L,M,N)
		call MullOneDimm(A,B,C,D,L,M,N)
		call D3dTo1d(D,Buff1,L,M,N)
		!call PrintMass(Buff1,L*M*N)
		do k=1,N
			do j=1,M
				do i=1,L
					write(*,*) Buff1(K+(J-1)*N+(I-1)*M*N),"   ,   ",D(i,j,k)
				enddo
			enddo
		enddo

		MidMass1=12
		MidMass3=12
		write(*,*)"-------------------------------------------------"
		call MullOneThreeDimmMassThreeCos(ResMass3,MidMass3,A,B,C,L,M,N)
		call MullOneThreeDimmMassOneCos(ResMass1,MidMass1,A,B,C,L,M,N)

		do k=1,N
			do j=1,M
				do i=1,L
					write(*,*) ResMass1(K+(J-1)*N+(I-1)*M*N),"   ,   ",ResMass3(i,j,k)
				enddo
			enddo
		enddo
		
		!call PrintMass3(D,L,M,N)
		read(*,*)
		deallocate(D, A, B, C, Buff1, Buff3, ResMass3, MidMass3, ResMass1, MidMass1)
    end program prog
