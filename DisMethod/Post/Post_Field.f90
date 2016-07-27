MODULE FIELD
    USE PROFILE
    USE MESH
    IMPLICIT NONE
    
    !VELOCITY & PRESSURE FIELD
    REAL, SAVE, ALLOCATABLE :: U(:, :, :), V(:, :, :), W(:, :, :), P(:, :, :), DIV(:, :, :)
    REAL, SAVE, ALLOCATABLE :: UC(:, :, :), VC(:, :, :), WC(:, :, :)
    REAL, SAVE, ALLOCATABLE :: VELABS(:, :, :), TOLPRE(:, :, :)
    
    !VORTEX FIELD
    REAL, SAVE, ALLOCATABLE :: VOR(:, :, :, :, :), VORAVE(:, :, :)
    REAL, SAVE, ALLOCATABLE :: VORX(:, :, :), VORY(:, :, :), VORZ(:, :, :)
    REAL, SAVE, ALLOCATABLE :: VORABS(:, :, :)
    
    !AVERAGE FIELD
    REAL, SAVE, ALLOCATABLE :: VELAVEX(:, :), VELAVEY(:, :), VELAVEZ(:, :)
    REAL, SAVE, ALLOCATABLE :: VORAVEX(:, :), VORAVEY(:, :), VORAVEZ(:, :)
    
    !FLUCTUATE VELOCITY FIELD
    REAL, SAVE, ALLOCATABLE :: DU(:, :, :), DV(:, :, :), DW(:, :, :)
    REAL, SAVE, ALLOCATABLE :: DURMS(:), DVRMS(:), DWRMS(:), DPRMS(:)
    
    !REYNOLDS STRESS
    REAL, SAVE, ALLOCATABLE :: RESAVEY(:, :, :)
    REAL, SAVE, ALLOCATABLE :: RES(:, :, :, :, :)
    
    !STRAIN RATE TENSOR
    REAL, SAVE, ALLOCATABLE :: SR(:, :, :, :, :)
    REAL, SAVE, ALLOCATABLE :: SRAVEY(:, :, :)
    
    !TURBULENT KINECTIC
    REAL, SAVE, ALLOCATABLE :: TK(:, :, :)
    
    !DISSAPATION TERM
    REAL, SAVE, ALLOCATABLE :: EPS(:, :, :)
    
    !Q-CRITIERA
    REAL, SAVE, ALLOCATABLE :: Q(:, :, :)
    
    !TOTAL SHEAR
    REAL, SAVE, ALLOCATABLE :: SHEAR_HIST(:)
    REAL, SAVE, ALLOCATABLE :: UTAU_HIST(:)
    
    !ALLOC IDENTIFIER
    LOGICAL, SAVE :: ALL_ALLOC = .FALSE.
    
    PRIVATE ALL_ALLOC
    
    CONTAINS
    
    SUBROUTINE CALCULATE(NUM)
        IMPLICIT NONE
        INTEGER NUM
        CHARACTER(LEN = 128) U_INFILE, P_INFILE, DIV_INFILE
        CHARACTER(LEN = 10) NUM_CHAR
        LOGICAL U_EXIST, P_EXIST
        INTEGER I, J, K
        
        U = 0
        V = 0
        W = 0
        U_EXIST = .FALSE.
        P_EXIST = .FALSE.
        
        WRITE(NUM_CHAR, "(I10)") NUM
        U_INFILE = TRIM(ADJUSTL(DATA_PATH))//'U_'//TRIM(ADJUSTL(NUM_CHAR))//'.DAT'
        P_INFILE = TRIM(ADJUSTL(DATA_PATH))//'P_'//TRIM(ADJUSTL(NUM_CHAR))//'.DAT'
        !DIV_INFILE = TRIM(ADJUSTL(DATA_PATH))//'DIV_'//TRIM(ADJUSTL(NUM_CHAR))//'.DAT'
        INQUIRE(FILE = U_INFILE, EXIST = U_EXIST)
        INQUIRE(FILE = P_INFILE, EXIST = P_EXIST)

        DO WHILE(.NOT. U_EXIST .OR. .NOT. P_EXIST)
            PRINT*, 'WAITING ...'
            CALL SLEEP(10)
            INQUIRE(FILE = U_INFILE, EXIST = U_EXIST)
            INQUIRE(FILE = P_INFILE, EXIST = P_EXIST)
        END DO
        
        OPEN(101, FILE = U_INFILE, STATUS = 'OLD', FORM = 'BINARY')
        READ(101) (((U(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        READ(101) (((V(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        READ(101) (((W(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        CLOSE(101)
        OPEN(102, FILE = P_INFILE, STATUS = 'OLD', FORM = 'BINARY')
        READ(102) (((P(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        CLOSE(102)
        
        CALL GET_VEL_FIELD
        CALL GET_VOR_FIELD
        CALL GET_AVE_FIELD
        CALL GET_FLUC_FIELD
        CALL GET_RES_FIELD
        CALL GET_SR_FIELD
        CALL GET_TURBE_FIELD
        
    END SUBROUTINE CALCULATE
    
    SUBROUTINE OUTPUT(NUM)
        IMPLICIT NONE
        
        INTEGER NUM
        CHARACTER(LEN = 128) PATH
        CHARACTER(LEN = 10) NUM_CHAR
        INTEGER I, J, K
        WRITE(NUM_CHAR, "(I10)") NUM
        
        PATH = TRIM(ADJUSTL(OUTPUT_PATH))//'FieldData_'//TRIM(ADJUSTL(NUM_CHAR))//'.DAT'
        OPEN(201, FILE = PATH, STATUS = 'REPLACE')
        WRITE(201, *) 'TITLE = "FIELD DATA"'
        WRITE(201, "(A256)") 'VARIABLES = "X" "Y" "Z" "U" "V" "W" "P" "VELABS" "TOLPRE" ' &
            //'"VORX" "VORY" "VORZ" "VORABS" "DU" "DV" "DW" ' &
            //'"TK" "Q"'
        WRITE(201, *) 'ZONE'
        WRITE(201, *) 'I = ', N1, ', J = ', N2, ', K = ', N3
        WRITE(201, *) 'DATAPACKING = BLOCK'
        WRITE(201, *) ((((X(I) + X(I - 1)) / 2, I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) ((((Y(J) + Y(J - 1)) / 2, I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) ((((Z(K) + Z(K - 1)) / 2, I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((UC(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((VC(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((WC(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((P(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((VELABS(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((TOLPRE(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((VORX(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((VORY(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((VORZ(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((VORABS(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((DU(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((DV(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((DW(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((TK(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((Q(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        CLOSE(201)
        CALL SYSTEM('preplot '//TRIM(ADJUSTL(PATH)))
        CALL SYSTEM('del /q/f/s '//TRIM(ADJUSTL(PATH)))
        
        PATH = TRIM(ADJUSTL(OUTPUT_PATH))//'AverageData_'//TRIM(ADJUSTL(NUM_CHAR))//'.DAT'
        OPEN(202, FILE = PATH, STATUS = 'REPLACE')
        WRITE(202, *) 'TITLE = "AVERAGE DATA"'
        WRITE(202, "(A256)") 'VARIABLES = "Y" "W12" "W13" "W23" "U_AVE" "V_AVE" "W_AVE" "VEL_AVE" '&
            //'"WX_AVE" "WY_AVE" "WZ_AVE" "VOR_AVE" "DURMS" "DVRMS" "DWRMS" "DPRMS" ' &
            //'"R11" "R12" "R13" "R22" "R23" "R33" "S11" "S12" "S13" "S22" "S23" "S33"'
        WRITE(202, *) 'ZONE'
        WRITE(202, *) 'I = ', N2
        WRITE(202, *) 'DATAPACKING = BLOCK'
        WRITE(202, *) ((Y(J) + Y(J - 1)) / 2, J = 1, N2)
        WRITE(202, *) (VORAVE(1, 2, J), J = 1, N2)
        WRITE(202, *) (VORAVE(1, 3, J), J = 1, N2)
        WRITE(202, *) (VORAVE(2, 3, J), J = 1, N2)
        WRITE(202, *) (VELAVEY(1, J), J = 1, N2)
        WRITE(202, *) (VELAVEY(2, J), J = 1, N2)
        WRITE(202, *) (VELAVEY(3, J), J = 1, N2)
        WRITE(202, *) (VELAVEY(4, J), J = 1, N2)
        WRITE(202, *) (VORAVEY(1, J), J = 1, N2)
        WRITE(202, *) (VORAVEY(2, J), J = 1, N2)
        WRITE(202, *) (VORAVEY(3, J), J = 1, N2)
        WRITE(202, *) (VORAVEY(4, J), J = 1, N2)
        WRITE(202, *) (DURMS(J), J = 1, N2)
        WRITE(202, *) (DVRMS(J), J = 1, N2)
        WRITE(202, *) (DWRMS(J), J = 1, N2)
        WRITE(202, *) (DPRMS(J), J = 1, N2)
        WRITE(202, *) (RESAVEY(1, 1, J), J = 1, N2)
        WRITE(202, *) (RESAVEY(1, 2, J), J = 1, N2)
        WRITE(202, *) (RESAVEY(1, 3, J), J = 1, N2)
        WRITE(202, *) (RESAVEY(2, 2, J), J = 1, N2)
        WRITE(202, *) (RESAVEY(2, 3, J), J = 1, N2)
        WRITE(202, *) (RESAVEY(3, 3, J), J = 1, N2)
        WRITE(202, *) (SRAVEY(1, 1, J), J = 1, N2)
        WRITE(202, *) (SRAVEY(1, 2, J), J = 1, N2)
        WRITE(202, *) (SRAVEY(1, 3, J), J = 1, N2)
        WRITE(202, *) (SRAVEY(2, 2, J), J = 1, N2)
        WRITE(202, *) (SRAVEY(2, 3, J), J = 1, N2)
        WRITE(202, *) (SRAVEY(3, 3, J), J = 1, N2)
        CLOSE(202)
        CALL SYSTEM('preplot '//TRIM(ADJUSTL(PATH)))
        CALL SYSTEM('del /q/f/s '//TRIM(ADJUSTL(PATH)))
    END SUBROUTINE OUTPUT
    
    SUBROUTINE ALLOC_FIELD
        IMPLICIT NONE
        IF(.NOT. ALL_ALLOC) THEN
            ALLOCATE(U(N1, 0:N2+1, N3))
            ALLOCATE(V(N1, 1:N2+1, N3))
            ALLOCATE(W(N1, 0:N2+1, N3))
            ALLOCATE(P(N1, N2, N3))
            ALLOCATE(DIV(N1, N2, N3))
            ALLOCATE(UC(N1, N2, N3))
            ALLOCATE(VC(N1, N2, N3))
            ALLOCATE(WC(N1, N2, N3))
            ALLOCATE(VELABS(N1, N2, N3))
            ALLOCATE(TOLPRE(N1, N2, N3))
            ALLOCATE(VOR(3, 3, N1, N2, N3))
            ALLOCATE(VORAVE(3, 3, N2))
            ALLOCATE(VORX(N1, N2, N3))
            ALLOCATE(VORY(N1, N2, N3))
            ALLOCATE(VORZ(N1, N2, N3))
            ALLOCATE(VORABS(N1, N2, N3))
            ALLOCATE(VELAVEX(4, N1))
            ALLOCATE(VELAVEY(4, N2))
            ALLOCATE(VELAVEZ(4, N3))
            ALLOCATE(VORAVEX(4, N1))
            ALLOCATE(VORAVEY(4, N2))
            ALLOCATE(VORAVEZ(4, N3))
            ALLOCATE(DU(N1, N2, N3))
            ALLOCATE(DV(N1, N2, N3))
            ALLOCATE(DW(N1, N2, N3))
            ALLOCATE(DURMS(N2))
            ALLOCATE(DVRMS(N2))
            ALLOCATE(DWRMS(N2))
            ALLOCATE(DPRMS(N2))
            ALLOCATE(RES(3, 3, N1, N2, N3))
            ALLOCATE(RESAVEY(3, 3, N2))
            ALLOCATE(SR(3, 3, N1, N2, N3))
            ALLOCATE(SRAVEY(3, 3, N2))
            ALLOCATE(TK(N1, N2, N3))
            ALLOCATE(EPS(N1, N2, N3))
            ALLOCATE(Q(N1, N2, N3))
            ALLOCATE(SHEAR_HIST((END_NUM - START_NUM) / STEP + 1))
            ALLOCATE(UTAU_HIST((END_NUM - START_NUM) / STEP + 1))
        
            ALL_ALLOC = .TRUE.
        END IF
    END SUBROUTINE ALLOC_FIELD
    
    SUBROUTINE DEALLOC_FIELD
        IMPLICIT NONE
        IF(ALL_ALLOC) THEN
            DEALLOCATE(U, V, W, P, DIV)
            DEALLOCATE(UC, VC, WC)
            DEALLOCATE(VELABS, TOLPRE)
            DEALLOCATE(VORX, VORY, VORZ, VORABS)
            DEALLOCATE(VELAVEX, VELAVEY, VELAVEZ)
            DEALLOCATE(VORAVEX, VORAVEY, VORAVEZ)
            DEALLOCATE(DU, DV, DW)
            DEALLOCATE(DURMS, DVRMS, DWRMS, DPRMS)
            DEALLOCATE(RES, RESAVEY)
            DEALLOCATE(SR, SRAVEY)
            DEALLOCATE(TK, EPS, Q, SHEAR_HIST, UTAU_HIST)
            
            ALL_ALLOC = .FALSE.
        END IF
    END SUBROUTINE DEALLOC_FIELD

    SUBROUTINE GET_VEL_FIELD()
        IMPLICIT NONE
        
        INTEGER I, J, K
        
        DO K = 1, N3
            DO J = 1, N2
                DO I = 1, N1
                    UC(I, J, K) = (U(IP(I), J, K) + U(I, J, K)) / 2
                    VC(I, J, K) = (V(I, JP(J), K) + V(I, J, K)) / 2
                    WC(I, J, K) = (W(I, J, KP(K)) + W(I, J, K)) / 2
                END DO
            END DO
        END DO
        
        VELABS = SQRT(UC ** 2 + VC ** 2 + WC ** 2)
        TOLPRE = P + VELABS ** 2 / 2
    END SUBROUTINE GET_VEL_FIELD
    
    SUBROUTINE GET_VOR_FIELD()
        IMPLICIT NONE
        
        REAL DWDY, DVDZ, DUDZ, DWDX, DVDX, DUDY
        REAL W1, W2, W3, U1, U2, U3, V1, V2, V3
        INTEGER I, J, K
        
        VOR = 0
        
        DO K = 1, N3
            DO J = 1, N2
                DO I = 1, N1
                    
                    !VORX = DWDY - DVDZ
                    W1 = (W(I, JP(J), K) + W(I, JP(J), KP(K))) / 2
                    W2 = (W(I, J, K) + W(I, J, KP(K))) / 2
                    W3 = (W(I, JM(J), K) + W(I, JM(J), KP(K))) / 2
                    DWDY = W1 * DYH(1, J) + W2 * DYH(2, J) + W3 * DYH(3, J)
                    V1 = VC(I, J, KM(K))
                    V2 = VC(I, J, K)
                    V3 = VC(I, J, KP(K))
                    DVDZ = (V3 - V1) / 2 / DZ
                    VORX(I, J, K) = DWDY - DVDZ
                    
                    !VORY = DUDZ - DWDX
                    DUDZ = (UC(I, J, KP(K)) - UC(I, J, KM(K))) / 2 / DX
                    DWDX = (WC(IP(I), J, K) - WC(IM(I), J, K)) / 2 / DX
                    VORY(I, J, K) = DUDZ - DWDX
                    
                    !VORZ = DVDX - DUDY
                    DVDX = (V(IP(I), J, K) - V(IM(I), J, K)) / 2 / DX
                    U1 = (U(I, JP(J), K) + U(IP(I), JP(J), K)) / 2
                    U2 = (U(I, J, K) + U(IP(I), J, K)) / 2
                    U3 = (U(I, JM(J), K) + U(IP(I), JM(J), K)) / 2
                    DUDY = U1 * DYH(1, J) + U2 * DYH(2, J) + U3 * DYH(3, J)
                    VORZ(I, J, K) = DVDX - DUDY
                    
                    VOR(1, 2, I, J, K) = DVDX - DUDY
                    VOR(1, 3, I, J, K) = DWDX - DUDZ
                    VOR(2, 1, I, J, K) = DUDY - DVDX
                    VOR(2, 3, I, J, K) = DWDY - DVDZ
                    VOR(3, 1, I, J, K) = DUDZ - DWDX
                    VOR(3, 2, I, J, K) = DVDZ - DWDY
                END DO
            END DO
        END DO
        
        VOR = VOR / 2
        DO J = 1, N2
            DO K = 1, 3
                DO I = 1, 3
                    VORAVE(I, K, J) = SUM(VOR(I, K, :, J, :)) / N1 / N3
                END DO
            END DO
        END DO
        
        VORABS = SQRT(VORX ** 2 + VORY ** 2 + VORZ ** 2)
    END SUBROUTINE GET_VOR_FIELD
    
    SUBROUTINE GET_AVE_FIELD()
        IMPLICIT NONE
        
        !X DIRECTION VELOCITY AVERAGE
        VELAVEX(1, :) = SUM(SUM(UC, 3), 2) / N2 / N3
        VELAVEX(2, :) = SUM(SUM(VC, 3), 2) / N2 / N3
        VELAVEX(3, :) = SUM(SUM(WC, 3), 2) / N2 / N3
        VELAVEX(4, :) = SUM(SUM(VELABS, 3), 2) / N2 / N3
        
        !Y DIRECTION VELOCITY AVERAGE
        VELAVEY(1, :) = SUM(SUM(UC, 1), 2) / N1 / N3
        VELAVEY(2, :) = SUM(SUM(VC, 1), 2) / N1 / N3
        VELAVEY(3, :) = SUM(SUM(WC, 1), 2) / N1 / N3
        VELAVEY(4, :) = SUM(SUM(VELABS, 1), 2) / N1 / N3
        
        !Z DIRECTION VELOCITY AVERAGE
        VELAVEZ(1, :) = SUM(SUM(UC, 1), 1) / N1 / N2
        VELAVEZ(2, :) = SUM(SUM(VC, 1), 1) / N1 / N2
        VELAVEZ(3, :) = SUM(SUM(WC, 1), 1) / N1 / N2
        VELAVEZ(4, :) = SUM(SUM(VELABS, 1), 1) / N1 / N2
        
        !X DIRECTION VORTICITY AVERAGE
        VORAVEX(1, :) = SUM(SUM(VORX, 2), 2) / N2 / N3
        VORAVEX(2, :) = SUM(SUM(VORY, 2), 2) / N2 / N3
        VORAVEX(3, :) = SUM(SUM(VORZ, 2), 2) / N2 / N3
        VORAVEX(4, :) = SUM(SUM(VORABS, 2), 2) / N2 / N3

        !Y DIRECTION VORTICITY AVERAGE
        VORAVEY(1, :) = SUM(SUM(VORX, 1), 2) / N1 / N3
        VORAVEY(2, :) = SUM(SUM(VORY, 1), 2) / N1 / N3
        VORAVEY(3, :) = SUM(SUM(VORZ, 1), 2) / N1 / N3
        VORAVEY(4, :) = SUM(SUM(VORABS, 1), 2) / N1 / N3

        !Z DIRECTION VORTICITY AVERAGE
        VORAVEZ(1, :) = SUM(SUM(VORX, 1), 1) / N1 / N2
        VORAVEZ(2, :) = SUM(SUM(VORY, 1), 1) / N1 / N2
        VORAVEZ(3, :) = SUM(SUM(VORZ, 1), 1) / N1 / N2
        VORAVEZ(4, :) = SUM(SUM(VORABS, 1), 1) / N1 / N2
    END SUBROUTINE GET_AVE_FIELD
    
    SUBROUTINE GET_FLUC_FIELD
        IMPLICIT NONE
        
        INTEGER J
        
        DO J = 1, N2
            DU(:, J, :) = UC(:, J, :) - VELAVEY(1, J)
            DV(:, J, :) = VC(:, J, :) - VELAVEY(2, J)
            DW(:, J, :) = WC(:, J, :) - VELAVEY(3, J)
            DURMS(J) = SQRT(SUM(DU(:, J, :) ** 2) / N1 / N3)
            DVRMS(J) = SQRT(SUM(DV(:, J, :) ** 2) / N1 / N3)
            DWRMS(J) = SQRT(SUM(DW(:, J, :) ** 2) / N1 / N3)
            DPRMS(J) = SQRT(SUM(P(:, J, :) ** 2) / N1 / N3)
        END DO
        
    END SUBROUTINE GET_FLUC_FIELD
    
    SUBROUTINE GET_RES_FIELD()
        IMPLICIT NONE
        
        INTEGER J
        
        RES(1, 1, :, :, :) = DU(:, :, :) * DU(:, :, :)
        RES(1, 2, :, :, :) = DU(:, :, :) * DV(:, :, :)
        RES(1, 3, :, :, :) = DU(:, :, :) * DW(:, :, :)
        RES(2, 1, :, :, :) = DV(:, :, :) * DU(:, :, :)
        RES(2, 2, :, :, :) = DV(:, :, :) * DV(:, :, :)
        RES(2, 3, :, :, :) = DV(:, :, :) * DW(:, :, :)
        RES(3, 1, :, :, :) = DW(:, :, :) * DU(:, :, :)
        RES(3, 2, :, :, :) = DW(:, :, :) * DV(:, :, :)
        RES(3, 3, :, :, :) = DW(:, :, :) * DW(:, :, :)
        
        DO J = 1, N2
            RESAVEY(1, 1, J) = SUM(RES(1, 1, :, J, :)) / N1 / N3
            RESAVEY(1, 2, J) = SUM(RES(1, 2, :, J, :)) / N1 / N3
            RESAVEY(1, 3, J) = SUM(RES(1, 3, :, J, :)) / N1 / N3
            RESAVEY(2, 1, J) = SUM(RES(2, 1, :, J, :)) / N1 / N3
            RESAVEY(2, 2, J) = SUM(RES(2, 2, :, J, :)) / N1 / N3
            RESAVEY(2, 3, J) = SUM(RES(2, 3, :, J, :)) / N1 / N3
            RESAVEY(3, 1, J) = SUM(RES(3, 1, :, J, :)) / N1 / N3
            RESAVEY(3, 2, J) = SUM(RES(3, 2, :, J, :)) / N1 / N3
            RESAVEY(3, 3, J) = SUM(RES(3, 3, :, J, :)) / N1 / N3
        END DO
        
    END SUBROUTINE GET_RES_FIELD
    
    SUBROUTINE GET_SR_FIELD
        IMPLICIT NONE
        
        REAL DUDX, DVDY, DWDZ
        REAL DUDY, DUDZ, DVDX, DVDZ, DWDX, DWDY
        REAL U1, U2, U3, W1, W2, W3, V1, V2, V3
        INTEGER I, J, K
        
        DO K = 1, N3
            DO J = 1, N2
                DO I = 1, N1
                    DUDX = (U(IP(I), J, K) - U(I, J, K)) / DX
                    DVDY = (V(I, JP(J), K) - V(I, J, K)) / DY(J)
                    DWDZ = (W(I, J, KP(K)) - W(I, J, K)) / DZ
                    
                    U1 = (U(I, JP(J), K) + U(IP(I), JP(J), K)) / 2
                    U2 = (U(I, J, K) + U(IP(I), J, K)) / 2
                    U3 = (U(I, JM(J), K) + U(IP(I), JM(J), K)) / 2
                    DUDY = U1 * DYH(1, J) + U2 * DYH(2, J) + U3 * DYH(3, J)
                    DUDZ = (UC(I, J, KP(K)) - UC(I, J, KM(K))) / 2 / DX
                    
                    DVDX = (V(IP(I), J, K) - V(IM(I), J, K)) / 2 / DX
                    DVDZ = (VC(I, J, KP(K)) - VC(I, J, KM(K))) / 2 / DZ
                    
                    DWDX = (WC(IP(I), J, K) - WC(IM(I), J, K)) / 2 / DX
                    W1 = (W(I, JP(J), K) + W(I, JP(J), KP(K))) / 2
                    W2 = (W(I, J, K) + W(I, J, KP(K))) / 2
                    W3 = (W(I, JM(J), K) + W(I, JM(J), KP(K))) / 2
                    DWDY = W1 * DYH(1, J) + W2 * DYH(2, J) + W3 * DYH(3, J)
                    
                    SR(1, 1, I, J, K) = DUDX + DUDX
                    SR(1, 2, I, J, K) = DUDY + DVDX
                    SR(1, 3, I, J, K) = DUDZ + DWDX
                    SR(2, 1, I, J, K) = DVDX + DUDY
                    SR(2, 2, I, J, K) = DVDY + DVDY
                    SR(2, 3, I, J, K) = DVDZ + DWDY
                    SR(3, 1, I, J, K) = DWDX + DUDZ
                    SR(3, 2, I, J, K) = DWDY + DVDZ
                    SR(3, 3, I, J, K) = DWDZ + DWDZ
                END DO
            END DO
        END DO
        
        SR = SR / 2
        
        DO J = 1, N2
            SRAVEY(1, 1, J) = SUM(SR(1, 1, :, J, :)) / N1 / N3
            SRAVEY(1, 2, J) = SUM(SR(1, 2, :, J, :)) / N1 / N3
            SRAVEY(1, 3, J) = SUM(SR(1, 3, :, J, :)) / N1 / N3
            SRAVEY(2, 1, J) = SUM(SR(2, 1, :, J, :)) / N1 / N3
            SRAVEY(2, 2, J) = SUM(SR(2, 2, :, J, :)) / N1 / N3
            SRAVEY(2, 3, J) = SUM(SR(2, 3, :, J, :)) / N1 / N3
            SRAVEY(3, 1, J) = SUM(SR(3, 1, :, J, :)) / N1 / N3
            SRAVEY(3, 2, J) = SUM(SR(3, 2, :, J, :)) / N1 / N3
            SRAVEY(3, 3, J) = SUM(SR(3, 3, :, J, :)) / N1 / N3
        END DO
        
    END SUBROUTINE GET_SR_FIELD
    
    SUBROUTINE GET_TURBE_FIELD()
        IMPLICIT NONE
        
        INTEGER I, J
        
        TK =(RES(1, 1, :, :, :) ** 2 &
            +RES(2, 2, :, :, :) ** 2 &
            +RES(3, 3, :, :, :) ** 2) / 2
        
        Q = 0
        DO J = 1, 3
            DO I = 1, 3
                Q = Q + (VOR(I, J, :, :, :) ** 2 - SR(I, J, :, :, :) ** 2) / 2
            END DO
        END DO
        
    END SUBROUTINE GET_TURBE_FIELD
    
    SUBROUTINE OUTPUT_DEBUG(NUM)
        IMPLICIT NONE
        
        INTEGER NUM
        CHARACTER(LEN = 128) PATH
        CHARACTER(LEN = 10) NUM_CHAR
        INTEGER I, J, K
        WRITE(NUM_CHAR, "(I10)") NUM
        
        !PATH = TRIM(ADJUSTL(OUTPUT_PATH))//'FieldData_'//TRIM(ADJUSTL(NUM_CHAR))//'.DAT'
        !OPEN(201, FILE = PATH, STATUS = 'REPLACE')
        PATH = TRIM(ADJUSTL(OUTPUT_PATH))//'FieldData.DAT'
        IF (NUM == START_NUM) THEN
            OPEN(201, FILE = PATH, STATUS = 'REPLACE')
            WRITE(201, *) 'TITLE = "FIELD DATA"'
            WRITE(201, "(A256)") 'VARIABLES = "X" "Y" "Z" "U" "V" "W" "P" "DIV" ' 
        ELSE
            OPEN(201, FILE = PATH, STATUS = 'OLD', POSITION = 'APPEND')
        ENDIF
        
        WRITE(201, *) 'ZONE'
        WRITE(201, *) 'I = ', N1, ', J = ', N2, ', K = ', N3
        WRITE(201, *) 'DATAPACKING = BLOCK'
        WRITE(201, *) 'SOLUTIONTIME = '//TRIM(ADJUSTL(NUM_CHAR))
        WRITE(201, *) ((((X(I) + X(I - 1)) / 2, I = 1, N1), J = 1, N2), K = 1, N3)
        DO K = 1, N3
            DO J = 1, N2
                DO I = 1, N1
                    WRITE(201, *) GETY(I, J, K, 1.0)
                END DO
            END DO
        END DO
        WRITE(201, *) ((((Z(K) + Z(K - 1)) / 2, I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((UC(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((VC(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((WC(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((P(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(201, *) (((DIV(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        CLOSE(201)
    END SUBROUTINE OUTPUT_DEBUG
    
    END MODULE FIELD
