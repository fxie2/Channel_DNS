PROGRAM MAIN
    USE PROFILE
    USE MESH
    USE FIELD
    IMPLICIT NONE
    
    INTEGER ITER, I
    CHARACTER(LEN = 128) PATH
    
    CALL READ_PROFILE()
    CALL INIT_MESH()
    CALL ALLOC_FIELD
    
    DO ITER = START_NUM, END_NUM, STEP
        PRINT*
        PRINT*, 'CONVERTING FILE NUMBER ', ITER
        CALL CALCULATE(ITER)
        CALL OUTPUT(ITER)
        SHEAR_HIST((ITER - START_NUM) / STEP + 1) = -RESAVEY(1, 2, 1) + SRAVEY(1, 2, 1) * 2 / RE
        UTAU_HIST((ITER - START_NUM) / STEP + 1) = (SRAVEY(1, 2, 1) + VORAVE(1, 2, 1)) / RE
    END DO
    
    PATH = TRIM(ADJUSTL(OUTPUT_PATH)) // 'SHEAR.DAT'
    OPEN(300, FILE = PATH, STATUS = 'REPLACE')
    WRITE(300, *) 'TITLE = "SHEAR HISTORY"'
    WRITE(300, *) 'VARIABLES = "T" "SHEAR" "U_TAU"'
    WRITE(300, *) 'ZONE'
    WRITE(300, *) 'I = ', SIZE(SHEAR_HIST)
    WRITE(300, *) 'DATAPACKING = BLOCK'
    WRITE(300, *) ((I), I = 1, SIZE(SHEAR_HIST))
    WRITE(300, *) (SHEAR_HIST(I), I = 1, SIZE(SHEAR_HIST))
    WRITE(300, *) (UTAU_HIST(I), I = 1, SIZE(SHEAR_HIST))
    CLOSE(300)
    CALL SYSTEM('preplot '//TRIM(ADJUSTL(PATH)))
    
    PRINT*, 'FINISHED'
    READ*
    END PROGRAM MAIN