
# newhybrids

This is the source code repository for the program *newhybrids* written by Eric C.
Anderson.  I don't have all the documentation up here yet.  So, just find the old
distribution at http://ib.berkeley.edu/labs/slatkin/eriq/software/software.htm and
use stuff from there, like the test data and the PreDefined Views files, etc,
but swap out the old executable, which is built for the Mac Classic
environment (which is really *quite* classic at this point in 2014)).

I have only tested this on Macs.

## Building the program
Here are the basic steps to build this on a Mac.  Note that you have to have the
developer tools installed (install XCode for free from the app store).  From the terminal
do the following:
```sh
# first clone it from GitHub
git clone https://github.com/eriqande/newhybrids.git
cd newhybrids
git submodule init
git submodule update

# then configure it without-x for the Mac and make it
./configure --without-x
make

# the above makes two binaries:
# newhybs    : newhybrids with the GLUT interface
# newhybsng  : newhybrids with no GLUT interface


# at this point, if you want to install these into
# /usr/local/bin you can do this:
sudo make install 
```

## Command line interface
Newhybrids has a nicer command line interface now than it did before.  To read about the 
available options you can do:
```sh
newhybs --help

# or

newhybs --help-full

# or if you are feeling fun you can do

newhybs --help-nroff | nroff -man | less
```