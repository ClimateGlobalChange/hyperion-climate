# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

BUILD_TARGETS= base pointwise
CLEAN_TARGETS= $(addsuffix .clean,$(BUILD_TARGETS))

# Base directory.
HYPERIONCLIMATEDIR= ..

# Load Makefile framework. 
include $(HYPERIONCLIMATEDIR)/mk/framework.make

HYPERIONCLIMATEBASEDIR= $(HYPERIONCLIMATEDIR)/src/base

HYPERIONCLIMATEBASELIB= $(HYPERIONCLIMATEBASEDIR)/libhyperionbase.a

HYPERIONCLIMATEPOINTWISELIB= $(HYPERIONCLIMATEDIR)/src/pointwise/libhyperionpointwise.a

HYPERIONCLIMATELIBS= $(HYPERIONCLIMATEBASELIB) $(HYPERIONCLIMATEPOINTWISELIB)

EXEC_FILES= recap.cpp

EXEC_TARGETS= $(EXEC_FILES:%.cpp=%)

FILES= $(UTIL_FILES) $(EXEC_FILES)

.PHONY: all clean $(BUILD_TARGETS)

# Build rules. 
all: $(BUILD_TARGETS) $(EXEC_TARGETS)

$(BUILD_TARGETS): %:
	cd $*; $(MAKE)

$(EXEC_TARGETS): %: $(UTIL_FILES:%.cpp=$(BUILDDIR)/%.o) $(BUILDDIR)/%.o $(HYPERIONCLIMATELIBS)
	$(CXX) $(LDFLAGS) -o $@ $(UTIL_FILES:%.cpp=$(BUILDDIR)/%.o) $(BUILDDIR)/$*.o $(HYPERIONCLIMATELIBS) $(LIBRARIES)
	mv $@ $(HYPERIONCLIMATEDIR)/bin

# Clean rules.
clean: $(CLEAN_TARGETS)
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

$(CLEAN_TARGETS): %.clean:
	cd $*; $(MAKE) clean

# Include dependencies.
-include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE
