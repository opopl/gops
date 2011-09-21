
      SUBROUTINE KEYWORD

      include "keyword.dec.i.f"
      include "keyword.init.i.f"

      OPEN (5,FILE='data',STATUS='OLD')

190   CALL INPUT(END)
      IF (.NOT. END) CALL READU(WORD)

      include "keyword.endstop.i.f"
      !include "keyword.r.chmd.i.f"
      include "keyword.if.i.f"
      
      CALL FLUSH(MYUNIT)
      GOTO 190

      RETURN
      END
