system = UNKNOWN

all:
	+make PRECISION=single -f Makefile.sd system=$(system)
	+make PRECISION=double -f Makefile.sd system=$(system)

clean:
	make PRECISION=single -f Makefile.sd clean system=$(system)
	make PRECISION=double -f Makefile.sd clean system=$(system)

single:
	+make PRECISION=single -f Makefile.sd system=$(system)

double:
	+make PRECISION=double -f Makefile.sd system=$(system)
