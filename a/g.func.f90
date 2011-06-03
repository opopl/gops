
MODULE FUNC

IMPLICIT NONE
SAVE

INTEGER :: NUMBER_OF_ATOMS

CONTAINS

SUBROUTINE OPENF(FILEHANDLE,MODE,FILENAME)

INTEGER FILEHANDLE
CHARACTER (LEN=*) MODE,FILENAME

SELECTCASE(MODE)
        CASE(">")
                OPEN(FILEHANDLE,FILENAME,STATUS="UNKNOWN",FORM="FORMATTED")
        CASE("<")
        CASE(">>")
                OPEN(FILEHANDLE,FILENAME,STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
ENDSELECT

ENDSUBROUTINE OPENF

ENDMODULE