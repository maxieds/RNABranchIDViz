# Makefile to copy external library dependencies and include files 

# The text interface command line parser.
CMD_LINE_PARSER = \
	src/ParseCommandLine.o \
	src/common_utils.o
	
# The utility that creates a structure image.
STRUCTURE_IMAGER = \
	src/StructureImageHandler.o \
	src/loop_utils.o
	
# Common files for the RNA library.
RNA_FILES = \
	RNA_class/RNA.o \
	RNA_class/thermodynamics.o \
	RNA_class/design.o \
	RNA_class/RsampleData.o \
	src/algorithm.o \
	src/alltrace.o \
	src/DynProgArray.o \
	src/dotarray.o \
	src/draw.o \
	src/extended_double.o \
	src/forceclass.o \
	src/MaxExpect.o \
	src/MaxExpectStack.o \
	src/outputconstraints.o \
	src/pfunction.o \
	src/probknot.o \
	src/random.o \
	src/rna_library.o \
	src/stackclass.o \
	src/stackstruct.o \
	src/stochastic.o \
	src/structure.o \
	src/substructure.o \
	src/TProgressDialog.o
	
# Common files for ProbScan library
PROBSCAN_FILES = \
    RNA_class/ProbScan.o

INCLUDEDIR=~/RNABranchIDViz/src/include-external/

copyexternal:
	make all
	gcc draw/DrawStructure.o -shared -o RNAStructureDraw.so \
		$(CMD_LINE_PARSER) $(STRUCTURE_IMAGER) $(RNA_FILES) $(PROBSCAN_FILES)
	mkdir -p ~/RNABranchIDViz/src/include-external
	cp RNAStructureDraw.so draw/DrawStructure.h src/*.h $(INCLUDEDIR)
	mkdir -p $(INCLUDEDIR)/RNA_class
	cp RNA_class/*.h $(INCLUDEDIR)/RNA_class/
	cd $(INCLUDEDIR) && sed -i 's/\.\.\/src\///g' *.h \
		&& sed -i 's/\.\.\/src\//\.\.\//g' RNA_class/*.h\
		&& sed -i 's/\.\.\/RNA_class/RNA_class/g' *.h RNA_class/*.h

