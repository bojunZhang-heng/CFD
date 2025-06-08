MYDIR = ./
FCOMPILER = gfortran -c -O3

OBJS = incomp2d.o fftcos.o fftpackt.o

incomp2d:	$(OBJS)
		gfortran -O3 -o incomp2d $(OBJS)

incomp2d.o:	$(MYDIR)incomp2d.f90
		$(FCOMPILER) $(MYDIR)incomp2d.f90
		
fftcos.o:	$(MYDIR)fftcos.f90
		$(FCOMPILER) $(MYDIR)fftcos.f90
		
fftpackt.o:     $(MYDIR)fftpackt.f90
		$(FCOMPILER) $(MYDIR)fftpackt.f90
