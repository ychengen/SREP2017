case = BGI_GLC_modelname_biome_min

FC = pgf90
INCLUDES = -I/sw/rhea/netcdf/4.3.3.1/rhel6.6_pgi14.4/include
LFLAGS = -L/sw/rhea/netcdf/4.3.3.1/rhel6.6_pgi14.4/lib -lnetcdff -lnetcdf

TARGET = regrid-$(case)
EXE = regrid-$(case).exe

all: $(TARGET)

$(TARGET): $(TARGET).f90
	$(FC) $(INCLUDES) $(LFLAGS) -o $(EXE) $(TARGET).f90

clean:
	rm -f $(EXE)
