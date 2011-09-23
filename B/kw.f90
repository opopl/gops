
      SUBROUTINE KEYWORD

      include "kwde.i.f90"
      include "kwi.i.f90"

      OPEN (5,FILE=D_FILE,STATUS='OLD')

190   CALL INPUT(END)
      IF (.NOT. END) CALL READU(WORD)

      include "kwes.i.f90"
      include "kwif.i.f90"
      
      CALL FLUSH(LFH)
      GOTO 190

      RETURN
      END
