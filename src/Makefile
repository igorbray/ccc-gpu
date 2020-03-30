system = UNKNOWN

all:
	-mkdir -p $(system)
	-mkdir -p $(system)/obj
	-mkdir -p $(system)/obj/double
	-mkdir -p $(system)/obj/single
	-mkdir -p $(system)/bin
	+make PRECISION=single -f Makefile.sd system=$(system)
	+make PRECISION=double -f Makefile.sd system=$(system)

clean:
	make PRECISION=single -f Makefile.sd clean system=$(system)
	make PRECISION=double -f Makefile.sd clean system=$(system)

single:
	-mkdir -p $(system)
	-mkdir -p $(system)/obj
	-mkdir -p $(system)/obj/single
	-mkdir -p $(system)/bin
	+make PRECISION=single -f Makefile.sd system=$(system)

double:
	-mkdir -p $(system)
	-mkdir -p $(system)/obj
	-mkdir -p $(system)/obj/double
	-mkdir -p $(system)/bin
	+make PRECISION=double -f Makefile.sd system=$(system)
