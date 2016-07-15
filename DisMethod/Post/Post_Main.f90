PROGRAM POST_MAIN
    IMPLICIT NONE
    
    CHARACTER(LEN = 128) ROOT_PATH, U_INFILE, P_INFILE
    CHARACTER(LEN = 128) MESH_PATH
    CHARACTER(LEN = 128) OUTPUT_PATH, U_OUTFILE
    CHARACTER(LEN = 10) NUM_CHAR
    REAL(8), ALLOCATABLE :: X(:), Y(:), Z(:)
    REAL(8), ALLOCATABLE :: U(:, :, :), V(:, :, :), W(:, :, :), P(:, :, :)
    REAL :: XC, YC, ZC
    INTEGER :: N1, N2, N3
    INTEGER :: START_NUM, END_NUM, PACE
    INTEGER :: FILE_NUM, I, J, K
    REAL(8), PARAMETER :: PI = 3.141592653589397
    PRINT*, 'ROOT_PATH : '
    READ*, ROOT_PATH
    PRINT*, 'MESH_PATH : '
    READ*, MESH_PATH
    PRINT*, 'OUTPUT_PATH : '
    READ*, OUTPUT_PATH
    PRINT*, 'N1, N2, N3 : '
    READ*, N1, N2, N3
    
    ALLOCATE(X(0:N1))
    ALLOCATE(Y(0:N2))
    ALLOCATE(Z(0:N3))
    ALLOCATE(U(N1, N2, N3))
    ALLOCATE(V(N1, N2, N3))
    ALLOCATE(W(N1, N2, N3))
    ALLOCATE(P(N1, N2, N3))
    
    DO I = 0, N1
        X(I) = 2 * PI / N1 * I
    END DO
    
    DO K = 0, N3
        Z(K) = 0.289 * PI / N3 * K
    END DO
    
    OPEN(100, FILE = MESH_PATH, STATUS = 'OLD')
    READ(100, *) (Y(J), J = 0, N2)
    CLOSE(100)
    
    ROOT_PATH = TRIM(ADJUSTL(ROOT_PATH))
    MESH_PATH = TRIM(ADJUSTL(MESH_PATH))
    PRINT*, 'START NUM, END NUM, PACE : '
    READ*, START_NUM, END_NUM, PACE
    
    DO FILE_NUM = START_NUM, END_NUM, PACE
        PRINT*, 'CONVERTING FILE NUMBER ', FILE_NUM
        WRITE(NUM_CHAR, "(I10)") FILE_NUM
        NUM_CHAR = TRIM(ADJUSTL(NUM_CHAR))
        U_INFILE = TRIM(ROOT_PATH)//'U_'//TRIM(NUM_CHAR)//'.DAT'
        P_INFILE = TRIM(ROOT_PATH)//'P_'//TRIM(NUM_CHAR)//'.DAT'
        OPEN(101, FILE = U_INFILE, STATUS = 'OLD', FORM = 'BINARY')
        READ(101) (((U(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        READ(101) (((V(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        READ(101) (((W(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        CLOSE(101)
        OPEN(102, FILE = P_INFILE, STATUS = 'OLD', FORM = 'BINARY')
        READ(102) (((P(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        CLOSE(102)
        
        U_OUTFILE = TRIM(OUTPUT_PATH) // 'UP_' // TRIM(NUM_CHAR) // '.DAT'
        
        OPEN(201, FILE = U_OUTFILE, STATUS = 'REPLACE')
        WRITE(201, *) 'TITLE = "VELOCITY AND PRESSURE"'
        WRITE(201, *) 'VARIABLES = "X" "Y" "Z" "U" "V" "W" "P"'
        WRITE(201, *) 'ZONE'
        WRITE(201, *) 'I = ', N1, ', J = ', N2, ', K = ', N3
        WRITE(201, *) 'DATAPACKING = BLOCK'
        WRITE(201, *) ((((X(I) + X(I - 1)) / 2, I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) ((((Y(J) + Y(J - 1)) / 2, I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) ((((Z(K) + Z(K - 1)) / 2, I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((U(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((V(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((W(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((P(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        CLOSE(201)
    END DO
    
    DEALLOCATE(U, V, W, P)
    DEALLOCATE(X, Y, Z)
END PROGRAM POST_MAIN