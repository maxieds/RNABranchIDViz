# RNABranchIDViz

<center>
<img src="https://github.com/maxieds/RNABranchIDViz/blob/master/images/d.16.e.E.cuniculi_nop-RNAStructViz.png" width="500" />
</center>

## Application Description

This utility is based on the [RNAStructViz](https://github.com/gtDMMB/RNAStructViz) source developed by 
researchers at GA Tech. It is a standalone application that is used to identify the four biological 
branch structures of an organism specified by an input 
[CT file](https://rna.urmc.rochester.edu/Text/File_Formats.html#CT). 
The usage of the program is as follows: 
```
Usage: RNABranchIDViz CTFileName [--quiet] [--debug] [--no-renumber-CT] [--no-images]
```
Most of the options specified above are self explanatory. By default, the utility identifies the four 
branches in the input CT file (first argument above) and then outputs re-numbered CT files for each 
branch which can be passed into other utilities, plots a color-coded circular image of each of the 
identified branches in the same format as output in the viewer for 
[RNAStructViz](https://github.com/gtDMMB/RNAStructViz), and plots each of the branches using a modification of 
the source code for the [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) project (*this last feature is currently 
under development*). To keep the original indices in the output branch CT files, you can specify the 
**--no-renumber-CT** option at the commandline. Similarly, to prevent the program from generating the image files 
for the identified branches, append the **--no-images** option to the command at runtime. The 
**--quiet** command can be used in batch or shell scripts to keep the program from producing any output at runtime. 

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
We will require the [ViennaRNA]https://www.tbi.univie.ac.at/RNA/() package to generate some of the diagrams in our utility. 
Here is how we configure it:
```
$ cd ~
$ wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.7.tar.gz
$ tar xvzf ViennaRNA-2.4.7.tar.gz
$ cd ViennaRNA-2.4.7
$ ./configure
$ make all
$ sudo make install
```
Now test that the ViennaRNA library has been installed correctly:
```
$ pkg-config --cflags --lib RNAlib2
-pthread -I/usr/local/include -I/usr/local/include/ViennaRNA -L/usr/local/lib -lRNA -fopenmp
```
If you obtain close to the output above, we are all set to begin compiling the local utility!

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
