all:
	cc -O2 -o gcell                                                         \
	cmain.c coper.c catom.c ccops.c ccout.c cimch.c cimol.c cipar.c cipdb.c \
	cpack.c cpacc.c cresi.c cmove.c cstat.c cclus.c cmvol.c 				\
	-march=native

linux:
	cc -O2 -o gcell                                                         \
	cmain.c coper.c catom.c ccops.c ccout.c cimch.c cimol.c cipar.c cipdb.c \
	cpack.c cpacc.c cresi.c cmove.c cstat.c cclus.c cmvol.c 				\
	-fPIC -mcmodel=medium -lm -march=native -g -ggdb

libgcell.so:
	cc -shared -o libgcell.so -O2 -fPIC -mcmodel=medium -march=native -g -ggdb \
	coper.c catom.c ccops.c ccout.c cimch.c cimol.c cipar.c cipdb.c \
	cpack.c cpacc.c cresi.c cmove.c cstat.c cclus.c cmvol.c c_interface.c \
	-lm

clean:
	rm -f gcell libgcell.so *.o
