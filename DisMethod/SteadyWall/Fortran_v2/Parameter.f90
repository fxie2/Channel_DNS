MODULE GLOBAL_PARAMETER
    IMPLICIT NONE
    
    REAL, PARAMETER :: PI = 3.141592653589793
    LOGICAL, SAVE :: INITIALIZED = .FALSE.
    
    !GEOMETRY PARAMETER
    INTEGER, SAVE :: N1, N2, N3
    REAL, SAVE :: LX, LY, LZ
    LOGICAL, SAVE :: USE_DEFAULT_LENGTH
    REAL, PARAMETER :: DEFAULT_LX = 3.14159265358
    REAL, PARAMETER :: DEFAULT_LY = 2.0
    REAL, PARAMETER :: DEFAULT_LZ = 0.289 * PI
    
    !TIME PARAMETER
    REAL, SAVE :: DT
    REAL, SAVE :: START_TIME, END_TIME
    REAL, SAVE :: T
    INTEGER, SAVE :: TOTAL_STEP_NUM
    INTEGER, SAVE :: CURNT_STEP_NUM
    LOGICAL, SAVE :: ADAPTED_TIME_STEP
    
    !IO PARAMETER
    LOGICAL, SAVE :: INIT_FROM_FILE
    LOGICAL, SAVE :: SAVE_STEP_DATA
    LOGICAL, SAVE :: PRNT_STEP_INFO
    CHARACTER(LEN = 128), SAVE :: INIT_FILE_PATH
    CHARACTER(LEN = 128), SAVE :: GRID_FILE_PATH
    CHARACTER(LEN = 128), SAVE :: SAVE_FILE_PATH
    CHARACTER(LEN = 128), SAVE :: LOG_FILE_PATH
    INTEGER, SAVE :: PRNT_PERIOD
    INTEGER, SAVE :: SAVE_PERIOD
    INTEGER, SAVE :: INSF_SAVE_PERIOD
    
    !FIELD PARAMETER
    LOGICAL, SAVE :: AVERAGE_FIELD
    LOGICAL, SAVE :: INST_FLOW      !INSTANTANEOUS FLOW FIELD
    REAL, SAVE :: RE
    REAL, SAVE :: MAXCFL
    REAL, SAVE :: INIT_TURB_INTENSITY
    REAL, PARAMETER :: DEFAULT_TURB_INTENSITY = 0.5
    
    !SOLVE PARAMETER
    INTEGER, SAVE :: MAX_SOLVE_ITER
    INTEGER, PARAMETER :: DEFAULT_MAX_ITER = 50
    REAL, SAVE :: MAX_SOLVE_ERR
    REAL, PARAMETER :: DEFAULT_MAX_ERR = 1E-10
    
    !BOUNDARY PARAMETER
    LOGICAL, SAVE :: UP_WAVE_WALL
    LOGICAL, SAVE :: DN_WAVE_WALL
    REAL, SAVE    :: DEVELOP_TIME
    INTEGER, SAVE :: UP_WAVE_NUMX, UP_WAVE_NUMZ
    INTEGER, SAVE :: DN_WAVE_NUMX, DN_WAVE_NUMZ
    REAL, SAVE    :: UP_WAVE_PSDX, UP_WAVE_PSDZ !PHASE SPEED
    REAL, SAVE    :: DN_WAVE_PSDX, DN_WAVE_PSDZ
    REAL, SAVE    :: MAX_UP_AMPX, MAX_UP_AMPZ   !AMPITUDE
    REAL, SAVE    :: MAX_DN_AMPX, MAX_DN_AMPZ
    
CONTAINS
    SUBROUTINE INIT_PARAMETERS(SOURCE)
        IMPLICIT NONE
        CHARACTER(LEN = *), OPTIONAL :: SOURCE
        IF(PRESENT(SOURCE)) THEN
            CALL READ_FROM_FILE(SOURCE)
        ELSE
            CALL USE_DEFAULT_PARAMETER()
        ENDIF
        INITIALIZED = .TRUE.
    END SUBROUTINE INIT_PARAMETERS
    
    SUBROUTINE READ_FROM_FILE(SOURCE)
        IMPLICIT NONE
        CHARACTER(LEN = *) :: SOURCE
        LOGICAL FILE_EXIST
        INQUIRE(FILE = SOURCE, EXIST = FILE_EXIST)
        IF(FILE_EXIST == .FALSE.) THEN
            PRINT*, 'PARAMETER FILE NOT EXIST IN PATH : ', SOURCE
        ELSE
            OPEN(100, FILE = SOURCE)
            READ(100, 110)
            READ(100, 110)
            READ(100, 110)
            READ(100, 120) N1
            READ(100, 120) N2
            READ(100, 120) N3
            READ(100, 130) RE
            READ(100, 130) LX
            READ(100, 130) LZ
            READ(100, 140) USE_DEFAULT_LENGTH
            READ(100, 130) START_TIME
            READ(100, 130) END_TIME
            READ(100, 130) DT
            READ(100, 120) TOTAL_STEP_NUM
            READ(100, 140) ADAPTED_TIME_STEP
            READ(100, 140) INIT_FROM_FILE
            READ(100, 140) SAVE_STEP_DATA
            READ(100, 140) PRNT_STEP_INFO
            READ(100, 120) PRNT_PERIOD
            READ(100, 120) SAVE_PERIOD
            READ(100, 130) INIT_TURB_INTENSITY
            READ(100, 130) MAXCFL
            READ(100, 140) AVERAGE_FIELD
            READ(100, 140) INST_FLOW
            READ(100, 120) INSF_SAVE_PERIOD
            READ(100, 120) MAX_SOLVE_ITER
            READ(100, 130) MAX_SOLVE_ERR
            READ(100, 140) UP_WAVE_WALL
            READ(100, 120) UP_WAVE_NUMX
            READ(100, 120) UP_WAVE_NUMZ
            READ(100, 130) UP_WAVE_PSDX
            READ(100, 130) UP_WAVE_PSDZ
            READ(100, 130) MAX_UP_AMPX
            READ(100, 130) MAX_UP_AMPZ
            READ(100, 140) DN_WAVE_WALL
            READ(100, 120) DN_WAVE_NUMX
            READ(100, 120) DN_WAVE_NUMZ
            READ(100, 130) DN_WAVE_PSDX
            READ(100, 130) DN_WAVE_PSDZ
            READ(100, 130) MAX_DN_AMPX
            READ(100, 130) MAX_DN_AMPZ
            READ(100, 130) DEVELOP_TIME
            READ(100, 110)
            READ(100, 150) INIT_FILE_PATH
            READ(100, 150) GRID_FILE_PATH
            READ(100, 150) SAVE_FILE_PATH
            READ(100, 150) LOG_FILE_PATH
            CLOSE(100)
        ENDIF
    
110         FORMAT(65X)
120         FORMAT(45X,I15)
130         FORMAT(45X,E15.7)
140         FORMAT(45X,L5)
150         FORMAT(40X,A128)    
    END SUBROUTINE READ_FROM_FILE
    
    SUBROUTINE USE_DEFAULT_PARAMETER()
        IMPLICIT NONE
        N1 = 16
        N2 = 128
        N3 = 32
        LX = DEFAULT_LX
        LY = DEFAULT_LY
        LZ = DEFAULT_LZ
        
        DT = 0.01
        START_TIME = 0
        TOTAL_STEP_NUM = 100000
        END_TIME = START_TIME + DT * TOTAL_STEP_NUM
        CURNT_STEP_NUM = 0
        T = DT * CURNT_STEP_NUM
        ADAPTED_TIME_STEP = .TRUE.
        
        INIT_FROM_FILE = .FALSE.
        SAVE_STEP_DATA = .TRUE.
        PRNT_STEP_INFO = .TRUE.
        GRID_FILE_PATH = 'channel.grd'
        SAVE_FILE_PATH = 'F:\DATA\'
        LOG_FILE_PATH = SAVE_FILE_PATH
        PRNT_PERIOD = 10
        SAVE_PERIOD = 100
        INSF_SAVE_PERIOD = 10
        
        AVERAGE_FIELD = .TRUE.
        INST_FLOW = .TRUE.
        MAXCFL = 2.0
        RE = 1000.0
        INIT_TURB_INTENSITY = DEFAULT_TURB_INTENSITY
        
        MAX_SOLVE_ITER = DEFAULT_MAX_ITER
        MAX_SOLVE_ERR  = DEFAULT_MAX_ERR
        
        UP_WAVE_WALL = .FALSE.
        DN_WAVE_WALL = .TRUE.
        UP_WAVE_NUMX = 0
        UP_WAVE_NUMZ = 0
        DN_WAVE_NUMX = 2
        DN_WAVE_NUMZ = 0
        UP_WAVE_PSDX = 0
        UP_WAVE_PSDZ = 0
        DN_WAVE_PSDX = 0
        DN_WAVE_PSDZ = 0
        MAX_UP_AMPX = 0
        MAX_UP_AMPZ = 0
        MAX_DN_AMPX = 0
        MAX_DN_AMPZ = 0
        DEVELOP_TIME = 1.0
    END SUBROUTINE USE_DEFAULT_PARAMETER
        
END MODULE GLOBAL_PARAMETER