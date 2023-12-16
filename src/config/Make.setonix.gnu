FC := ftn -fbacktrace #-ffpe-trap=invalid -C -g 
FCFLAGS_OPT2 := -dynamic  -cpp -O2 -fopenmp -fno-signed-zeros -fno-trapping-math  -mcmodel=large -std=legacy -fcray-pointer -fno-range-check #-Wno-missing-include-dirs #-fallow-argument-mismatch
FCFLAGS_OPT1 := -dynamic -cpp -O2 -fopenmp -fno-signed-zeros -fno-trapping-math -mcmodel=large -std=legacy -fcray-pointer -fno-range-check #-Wno-missing-include-dirs #-fallow-argument-mismatch
FCFLAGS_FAST := -ffast-math 
FCFLAGS_GPU := -dynamic -cpp -O3 -fopenmp -ffast-math -mcmodel=large #-fallow-argument-mismatch #-fno-signed-zeros -fno-trapping-math -mcmodel=large #-fallow-argument-mismatch
FCFLAGS_SAVE := -fno-automatic
FCFLAGS_REAL8 := -fdefault-real-8 -fdefault-double-8
FCFLAGS_MODFLAG := -J
LDFLAGS := #-lopenblas -llapack  -lscalapack 


#FCFLAGS_OPT2 := -dynamic  -cpp -O0 -fopenmp -mcmodel=large -std=legacy -Wno-missing-include-dirs -fcray-pointer
#FCFLAGS_OPT1 := -dynamic -cpp -O0 -fopenmp -mcmodel=large -std=legacy -Wno-missing-include-dirs -fcray-pointer
