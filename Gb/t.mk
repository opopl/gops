
gbln: $(BASE_OBJS) $(BLN_OBJS)
	$(FC) $(FFLAGS) $(SEARCH_PATH) -o $@ $(BASE_OBJS) $(BLN_OBJS) $(LDFLAGS) $(LIBS)
	cp gbln $(BINPATH)

