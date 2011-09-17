
gbln: $(BLN_OBJS)
	$(FC) $(FFLAGS) $(SEARCH_PATH) -o $@ $(BLN_OBJS) $(LDFLAGS) $(LIBS)
	cp gbln $(BINPATH)
