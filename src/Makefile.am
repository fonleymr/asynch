noinst_LIBRARIES = libasynch.a
libasynch_a_SOURCES = rkmethods.c problems.c blas.c vector.c riversys.c sort.c comm.c system.c processdata.c partition.c builtin.c misc.c rainfall.c solvers.c io.c forcings.c compression.c asynch_interface.c modeloutputs.c data_types.c consistency.c \
  rksteppers.c \
  steppers/explicit.c
#  steppers/explicit_discont.c \
#  steppers/explicit_index1.c \
#  steppers/explicit_index1_dam.c

bin_PROGRAMS = asynch
asynch_SOURCES = optparse.c asynchdist.c
asynch_LDADD = libasynch.a $(HDF5_LIBS) $(POSTGRESQL_LIBS)
asynch_LDFLAGS = $(HDF5_LDFLAGS) $(POSTGRESQL_LDFLAGS)

include_HEADERS = structs.h structs_fwd.h asynch_interface.h

AM_CFLAGS = $(HDF5_CPPFLAGS) $(POSTGRESQL_CPPFLAGS)
