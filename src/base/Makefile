# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Base directory.
HYPERIONCLIMATEDIR= ../..

# Load Makefile framework. 
include $(HYPERIONCLIMATEDIR)/mk/framework.make

CXXFLAGS+=-I$(HYPERIONCLIMATEDIR)/src/netcdf-cxx-4.2

FILES= Announce.cpp \
	   Exception.cpp \
       NetCDFUtilities.cpp \
       TimeObj.cpp \
	   MultiTypeDataArray2D.cpp \
	   MultiTypeDataHeaders.cpp \
	   Object.cpp \
	   RecapConfigObject.cpp \
	   GridObject.cpp \
	   FileListObject.cpp \
	   VariableLookupObject.cpp \
	   GlobalFunction.cpp \
       Variable.cpp \
       kdtree.cpp \
	   LegendrePolynomial.cpp \
	   GaussQuadrature.cpp \
	   GaussLobattoQuadrature.cpp \
	   TriangularQuadrature.cpp \
	   PolynomialInterp.cpp \
	   FiniteElementTools.cpp \
	   MeshUtilities.cpp \
	   MeshUtilitiesFuzzy.cpp \
	   GridElements.cpp \
	   ShpLoader.cpp

LIB_TARGET= libhyperionbase.a

.PHONY: all clean

# Build rules. 
all: $(LIB_TARGET)

$(LIB_TARGET): $(FILES:%.cpp=$(BUILDDIR)/%.o)
	rm -f $(LIB_TARGET)
	ar -cqs $(LIB_TARGET) build/*.o

# Clean rules.
clean:
	rm -f $(LIB_TARGET) 
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

# Include dependencies.
-include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE
