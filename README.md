# RNABranchIDViz

## Application Description

## Installation Notes

### Cloning the source

We begin by cloning the (currently private) repo into the user's home directory:
```
$ cd ~
$ git clone https://github.com/maxieds/RNABranchIDViz.git
```

### Dependencies 

This utility is based on the source for [RNAStructViz](https://github.com/gtDMMB/RNAStructViz). 
The complications that can come with the latter utility's dependency on ``FLTK`` have been 
removed from this standalone commandline utility. In it's place we require ``cairo2`` which is 
easily installed on modern Linux systems:
```
$ sudo apt-get install libcairo2-dev
``` 
We will also require the use of the [RNAstructure library](https://rna.urmc.rochester.edu/RNAstructure.html), or 
at least its ``draw`` package and shared libraries. We first download the library from the official site and 
extract it to our working directory:
```
$ cd ~
$ export RNASTRUCTZIP=/path/to/RNAstructureSource.tgz
$ tar xvzf $RNASTRUCTZIP
$ cd RNAstructure
``` 
Next we will need to build the ``draw`` utility and compile it as a shared C++ library so we 
can link against it when we compile our program. We do this and copy the corresponding header files 
into our local source repo directory as follows: 
```
$ make -f ~/RNABranchIDViz/CopyExternal.Makefile copyexternal
$ sudo mkdir -p /usr/local/RNAstructure
$ sudo cp ~/RNABranchIDViz/src/include-external/RNAStructureDraw.so /usr/local/RNAstructure/
```
Now we are all set to begin compiling the local utility!

### Compiling from source on Linux (Debian-Variants)

First, change back into the working source directory from above:
```
$ cd ~/RNABranchIDViz
```
On most recent Debian and Ubuntu desktops it should not be necessary for users to have to 
re-run the build scripts to run the ``aclocal`` and ``auto*`` family of configuration 
utilities. If you are running a version of ``aclocal<1.15`` however, the next step may be 
necessary for you to proceed with the install:
```
$ ./buildconfig.sh
```
If this step succeeds you can build the program with the following sequence of commands:
```
$ ./configure
$ make
$ echo "alias RNABranchIDViz=\'~/RNABranchIDViz/src/RNABranchIDViz\'" >> ~/.bashrc
$ source ~/.bashrc
$ RNABranchIDViz
``` 
If all goes well in the above compilation steps, you now have a working binary for the 
utility built and ready to use at your terminal the next time you login!
