!op226 {{{
      IF (ADDPATH) THEN
         CALL MYSYSTEM(STATUS,DEBUG,'cp ' // TRIM(ADJUSTL(PATHNAME)) // ' path.info')
         IF (ADDTRIPLES) THEN
            CALL GETALLPATHS
         ELSE
            CALL GETNEWPATH(0,0)
         ENDIF
      ENDIF
!op226 }}
