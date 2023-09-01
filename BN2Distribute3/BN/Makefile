# 
# ------ MAKEFILE FILE for ALPHA ----------------
# 

XMLROOT = /home/jnr3hh/BN2Distribute/xerces-c-src_2_8_0

#INCLUDES = -I. -I $(XMLROOT)/include -I /usr/local/include/g++-v3/
INCLUDES = -I. -I $(XMLROOT)/include 

CC	= g++

AR      = ar
 
CFLAGS =  $(INCLUDES) -O -DLINUX -c 

LFLAGS =  -L$(XMLROOT)/lib -lpthread -lxerces-c -L/usr/lib64 -O3 -ffast-math -m64 
 
PGM    = testBN

OBJS    = BN.o Node.o Relationship.o BNHandler.o eQTL.o Gaussian.o Graph.o Matrix.o nrutil.o common.o 

.PRECIOUS :   $(LIB) $(PGM)

.SUFFIXES: .C

.C.o:
	$(CC) $(CFLAGS) $<

include Makefile.inc

$(PGM) :  $(PGM).o $(OBJS)
	@echo "linking $(PGM) ..."
	$(CC) $(PGM).o $(OBJS) $(LFLAGS)-o $@

SAX2Test: SAX2Test.o BNHandler.o
	@echo "linking $< ...."
	$(CC) SAX2Test.o BNHandler.o -o $@ $(LFLAGS) 

cmi: cmi.o $(OBJS)
	@echo "linking $< ...."
	$(CC) cmi.o $(OBJS) $(LFLAGS) -o $@

printBICs: printBICs.o $(OBJS)
	@echo "linking $< ...."
	$(CC) printBICs.o $(OBJS) $(LFLAGS) -o $@

testBN_tmp: testBN_tmp.o $(OBJS)
	@echo "linking $< ...."
	$(CC) testBN_tmp.o $(OBJS) $(LFLAGS) -o $@
