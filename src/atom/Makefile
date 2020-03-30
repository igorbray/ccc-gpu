#a = accobm.o akvint.o amod1.o coef.o cqd1.o cqd11.o cqd2.o dh1.o direct.o exchen.o f1.o f3.o f4.o fchf.o fcross.o g1.o g3.o g5.o ic.o ifact.o int.o inver.o invert.o iot3.o iot6.o ipari.o komon.o komon1.o komon2.o mael1.o mael2.o mael3.o mael4.o mael5.o maelcn.o main1.o main2.o matr.o mezero.o mxmlt.o offdia.o phase.o phoabs.o pot.o printm.o rk.o scatwr.o scfhf.o selfen.o simps.o simps1.o tape.o var1.o wrico.o yk.o zero.o

SOURCE = accobm.f akvint.f amod1.f coef.f cqd1.f cqd11.f cqd2.f dh1.f direct.f exchen.f f1.f f3.f f4.f fchf.f fcross.f g1.f g3.f g5.f ic.f ifact.f int.f inver.f invert.f iot3.f iot6.f ipari.f komon.f komon1.f komon2.f mael1.f mael2.f mael3.f mael4.f mael5.f maelcn.f mm12.f matr.f mezero.f mxmlt.f offdia.f phase.f phoabs.f pot.f printm.f rk.f scatwr.f scfhf.f selfen.f simps.f simps1.f tape.f var1.f yk.f zero.f

OBJS=$(SOURCE:%.f=$(OBJDIR)/%.o)


#e = -O2 -openmp -traceback -mcmodel=medium

$(OBJDIR)/hflib.a: $(OBJS)
	ar crs $(OBJDIR)/hflib.a $(OBJS)
#	mv hflib.a $(OBJDIR)
#	ranlib hflib.a

clean:
	rm -f $(OBJDIR)/*
	rm -f *.o
	rm -f *.mod
	rm -f *.a

$(OBJDIR):
	if ( test ! -d $@ ) then ( mkdir -p $@ ) fi

$(OBJDIR)/coef.o : coef.f paratom.f 
	$(FC) -c $a coef.f  -o $@

$(OBJDIR)/exchen.o : exchen.f paratom.f 
	$(FC) -c $a exchen.f  -o $@

$(OBJDIR)/fchf.o : fchf.f paratom.f 
	$(FC) -c $a fchf.f  -o $@


$(OBJDIR)/int.o : int.f paratom.f 
	$(FC) -c $a int.f  -o $@


$(OBJDIR)/%.o : %.f 
	$(FC) $a -c $< -o $@
