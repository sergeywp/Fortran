module a
	contains 
	subroutine ss
				!attributes(global) subroutine gpu_add(A,B,C,N)
		!real(4), device :: A(N), B(N), C(N)
		!integer, value :: N
		!integer :: i
		!i = (blockidx%x-1)*32+threadidx%x
		!if(i<=N)  C(i) = A(i) + B(i)
	!end subroutine
	!do k3=1,NI
	!	do j3=1,NJ
	!		do i3=1,NK
	!			write (*,*)i3," H ", j3,"  "," H ", k3, GAMI_h(i3,j3,k3)
	!		enddo
	!	enddo
	!enddo
	
	!public static void MullOneDim(int []A,int[] B, int[] C, int[][][] D){
	!	int[] Buff = new int[L*N*M];
	!	for(int i=0; i<L*N*M; i++){
	!		Buff[i]=A[i/(N*M)]*B[(i/N)%M]*C[i%N];
	!		//System.out.print(Buff[i]+" ");
	!	}
	!	//PrintMass(Buff);
	!	for(int i=0; i<L; i++){
	!		for(int j=0; j<M; j++){
	!			for(int k=0; k<N; k++){ 
	!				D[i][j][k] = Buff[k+j*N+i*M*N];
	!			}
	!		}	
	!	}
	!
	!}
	!FillMassNotSim(A,B,C);
	!	//PrintMass(A);
	!	//PrintMass(B);
	!	//PrintMass(C);
	!	MullNorm(A, B, C, D);
	!	//PrintMass3(D);
	!	MullOneDim(A, B, C, D);
	!	//System.out.println("\n");
	!	PrintMass3(D);
	!	int [] K1 = new int[L*M*N];
	!	int [][][] K3 =new int[L][M][N];
	!	K1 = d3dTo1d(D, L, M, N);
	!	System.out.println("CHECK:\n");
	!	PrintMass(K1);
	!	System.out.println("\n");
	!	K3 = d1dTo3d(K1, L, M, N);
	!	PrintMass3(K3);

	!public static void MullOneDim(int []A,int[] B, int[] C, int[][][] D){
	!	int[] Buff = new int[L*N*M];
	!	for(int i=0; i<L*N*M; i++){
	!		Buff[i]=A[i/(N*M)]*B[(i/N)%M]*C[i%N];
	!		//System.out.print(Buff[i]+" ");
	!	}
	!	
	!}

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
	! 
	!
	!
	! AJP(I,J,K) = GAMJ(I,J+1,K)*XCV(I)*ZCV(K)/YDIF(J+1)
	!
	!
	!
	end subroutine ss
end module a