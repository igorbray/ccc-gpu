SOURCE = modules.f mainhe.f hevmat.f setmax.f setmax_dif.f hel-3.f eigv.f ort.f osc.f rearrange.f r2.f scat.f structure.f rsg.f hforbnum.f e2eovp.f cfhe.f gsort.f nonuniq.f so.f so_1el.f posHeVmat.f
#SOURCE = CI_MODULE.f VDCORE_MODULE.f PO_MODULE.f DM_MODULE.f main.f hevmat.f setmax.f hel-3.f eigv.f ort.f osc.f rearrange.f r2.f scat.f structure.f rsg.f hforbnum.f e2e.f cfhe.f gsort.f state.f nonuniq.f so.f so_1el.f

OBJS=$(SOURCE:%.f=$(OBJDIR)/%.o)

#a = -O3 -fpp -openmp -mcmodel=medium
#a = -C -g -fpp -openmp -traceback

$(OBJDIR)/helib.a: $(OBJS)
	ar crs $(OBJDIR)/helib.a $(OBJS)
#	mv helib.a $(OBJDIR)
#	ranlib helib.a


clean:
	rm -f $(OBJDIR)/*
	rm -f *.o
	rm -f *.mod
	rm -f helib.a

.PHONY: $(OBJDIR)
$(OBJDIR):
	if ( test ! -d $@ ) then ( mkdir -p $@ ) fi

$(OBJDIR)/hevmat.o : hevmat.f $(OBJDIR)/modules.o
	$(FC) $(aGPU) -c hevmat.f -o $@

#$(OBJDIR)/modules.o:  modules.f
#	$(FC) $a -c modules.f -o $@

#$(OBJDIR)/so.o:  so.f par.f state.f modules.f
#	$(FC) $a -c so.f -o $@

$(OBJDIR)/nonuniq.o:  nonuniq.f par.f modules.f
	$(FC) $a -c  nonuniq.f -o $@

$(OBJDIR)/posHeVmat.o:  posHeVmat.f par.f par.pos
	$(FC) $a -c  posHeVmat.f -o $@

$(OBJDIR)/scat.o:  scat.f par.f par.pos
	$(FC) $(aGPU) -c  scat.f -o $@

$(OBJDIR)/%.o : %.f par.f
	$(FC) $a -c $< -o $@
