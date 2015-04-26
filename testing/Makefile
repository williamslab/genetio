GPP=g++
INBENCH=benchmark.cc
OUTBENCH=benchmark


LIBDIR=../
LIBNAME=genetio
INTEST=testing.cc
OUTTEST=test

benchmark:
	$(GPP) -O2 $(INBENCH) -o $(OUTBENCH)
	chmod +x $(OUTBENCH)

test:
	$(GPP) -O2 $(INTEST) -o $(OUTTEST) -L $(LIBDIR) -l$(LIBNAME)

clean:
	rm -f $(OUTBENCH)
	rm -f $(OUTTEST)
