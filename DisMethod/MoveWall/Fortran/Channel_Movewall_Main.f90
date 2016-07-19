PROGRAM CHANNEL_MAIN
    USE GLOBAL_PARAMETER
    USE MESH
    USE FIELD
    USE MATH
    IMPLICIT NONE
    CHARACTER(LEN = 128) SOURCE
    
    CALL INIT_PARAMETERS()
    
    !INIT MESH
    CALL NEW_MESH()
    CALL INIT_MESH()
    !INIT FIELD
    CALL ALLOC_FIELD()
    CALL INIUP()
    !CALL INIUP_FILE('E:\data\')
    CALL CHECK_DIV
    CALL CHECK_FLOW_RATE
    PRINT*, XFLOW, ZFLOW
    PRINT*, DIVMAX
    
    !OPEN LOG FILE
    !OPEN(200, FILE = TRIM(ADJUSTL(LOG_FILE_PATH))//'LOG.TXT', STATUS = 'REPLACE')
    
    DO WHILE(CURNT_STEP_NUM < TOTAL_STEP_NUM)
        
        CALL SOLVEUP()
        PRINT*, 'ITER : ' , SOLVE_ITER
        IF(MOD(CURNT_STEP_NUM, 1) == 0) THEN
            CALL CHECK_DIV()
            CALL CHECK_FLOW_RATE()
            
            WRITE(*, 100) 'TIME   : ', T
            WRITE(*, 100) 'DIVMAX : ', DIVMAX
            WRITE(*, 100) 'PGX    : ', PGX
            WRITE(*, 100) 'PGZ    : ', PGZ
            WRITE(*, 100) 'XFLOW  : ', XFLOW
            WRITE(*, 100) 'ZFLOW  : ', ZFLOW
            WRITE(*, 100) 'ITER   : ', DBLE(SOLVE_ITER)
            WRITE(*, 110) REPEAT('-', 10), T / END_TIME * 100, '%', REPEAT('-', 10)
            
            WRITE(200, 100) 'TIME   : ', T
            WRITE(200, 100) 'DIVMAX : ', DIVMAX
            WRITE(200, 100) 'PGX    : ', PGX
            WRITE(200, 100) 'PGZ    : ', PGZ
            WRITE(200, 100) 'XFLOW  : ', XFLOW
            WRITE(200, 100) 'ZFLOW  : ', ZFLOW
            WRITE(200, *)
            
100         FORMAT(A10, E15.8)
110         FORMAT(A10, F4.1, A1, A10)            
            PRINT*
        END IF
        
        !CALL OUTPUT
        !IF(MOD(CURNT_STEP_NUM, SAVE_PERIOD) == 0) CALL OUTPUT
        PRINT*, MAXVAL(DIV), MAXLOC(DIV)

        CURNT_STEP_NUM = CURNT_STEP_NUM + 1
        T = T + DT
    END DO
    
    !FINALIZE MESH
    CALL DEL_MESH()
    !FINALIZE FIELD
    CALL DEALLOC_FIELD()
    CALL DEALLOC_FAC()
    CALL DELETEVEC()
    CALL FFT_DEALLOC()
    END PROGRAM CHANNEL_MAIN