# ===============================
# Makefile for LAPPDreco
# ===============================


ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

# Linux with egcs
INC := -I include
CC := g++
CFLAGS := -g -Wall -Werror -fPIC $(ROOTCFLAGS)
LIBS := $(ROOTLIBS) 

SRCDIR := src
BUILDDIR := build

SOURCES := $(shell find $(SRCDIR) -type f -name *.cpp)
OBJS := $(patsubst $(SRCDIR), $(BUILDDIR), $(SOURCES:.CPP=.O))

NAME := lappd
#------------------------------------------------------------------------------


all: $(NAME)

$(NAME): $(OBJS)
	@echo "Linking..."
	@$(CC) $(CFLAGS) $(INC) -o $@ $^ $(LIBS)
	@echo "Done"

$(BUILDDIR)%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	@$(CC) $(CFLAGS) $(INC) -o $@ -c $<

clean:
	@rm -rf $(BUILDDIR) $(NAME)

.PHONY: clean

# clean:
# 		@rm -f $(OBJS) lappd.$(ObjSuf) *\~ core SplineFitDict*
# 		@rm -f $(PROGRAM) $(PROGRAMSO)
# 		@rm -f $(ROOTDEV)/bin/$(PROGRAM)
# 		@rm -f $(ROOTDEV)/lib/$(LIBNAME).$(DllSuf)
# 		@rm -f $(ROOTDEV)/include/TPoly3.h
# 		@rm -f $(ROOTDEV)/include/TOnePadDisplay.h
# 		@rm -f $(ROOTDEV)/include/TBandedLE.h
# 		@rm -f $(ROOTDEV)/include/TZigZag.h
# 		@rm -f $(ROOTDEV)/include/TSplineFit.h
# 		@rm -f $(ROOTDEV)/include/UVevent.h
# 		@rm -f $(ROOTDEV)/include/TemplateFit.h
# 		@rm -f $(ROOTDEV)/include/vector.h #added
# 		@rm -f $(ROOTDEV)/include/matrix.h #added
# 		@rm -f $(ROOTDEV)/include/denseMatrix.h #added
# 		@rm -f $(ROOTDEV)/include/sparseMatrix.h #added
# 		@rm -f $(ROOTDEV)/include/nnls.h #added
# 		@rm -f UVevent.o UVeventDict.o UVeventDict.cxx UVeventDict.h libUVevent.so nnls.o nnlsDriver.o sparseMatrix.o denseMatrix.o #added
# 		@rm -f TemplateFit.o TemplateFitDict.o TemplateFitDict.cxx TemplateFitDict.h libTemplateFit.so
#
# install:
# 		@rm -f $(ROOTDEV)/include/UVevent.h
# 		@rm -f $(ROOTDEV)/include/TemplateFit.h
# 		@rm -f $(ROOTDEV)/bin/$(PROGRAM)
# 		@rm -f $(ROOTDEV)/lib/$(LIBNAME).$(DllSuf)
# 		@rm -f $(ROOTDEV)/include/TPoly3.h
# 		@rm -f $(ROOTDEV)/include/TOnePadDisplay.h
# 		@rm -f $(ROOTDEV)/include/TBandedLE.h
# 		@rm -f $(ROOTDEV)/include/TZigZag.h
# 		@rm -f $(ROOTDEV)/include/TSplineFit.h
# 		@rm -f $(ROOTDEV)/include/vector.h #added
# 		@rm -f $(ROOTDEV)/include/matrix.h #added
# 		@rm -f $(ROOTDEV)/include/denseMatrix.h #added
# 		@rm -f $(ROOTDEV)/include/sparseMatrix.h #added
# 		@rm -f $(ROOTDEV)/include/nnls.h #added
# 		@cp UVevent.h $(ROOTDEV)/include/UVevent.h
# 		@cp TemplateFit.h $(ROOTDEV)/include/TemplateFit.h
# 		@cp $(PROGRAM) $(ROOTDEV)/bin/$(PROGRAM)
# 		@cp $(LIBNAME).$(DllSuf) $(ROOTDEV)/lib/$(LIBNAME).$(DllSuf)
# 		@cp TPoly3.h $(ROOTDEV)/include/TPoly3.h
# 		@cp TOnePadDisplay.h $(ROOTDEV)/include/TOnePadDisplay.h
# 		@cp TBandedLE.h $(ROOTDEV)/include/TBandedLE.h
# 		@cp TZigZag.h $(ROOTDEV)/include/TZigZag.h
# 		@cp TSplineFit.h $(ROOTDEV)/include/TSplineFit.h
# 		@cp UVevent.h $(ROOTDEV)/include/UVevent.h
# 		@cp TemplateFit.h $(ROOTDEV)/include/TemplateFit.h
# 		@cp vector.h $(ROOTDEV)/include/vector.h #added
# 		@cp vector.h $(ROOTDEV)/include/matrix.h #added
# 		@cp vector.h $(ROOTDEV)/include/denseMatrix.h #added
# 		@cp vector.h $(ROOTDEV)/include/sparseMatrix.h #added
# 		@cp vector.h $(ROOTDEV)/include/nnls.h #added
# 		@echo "libraries, shared libraries and includes copied to $(ROOTDEV)"
