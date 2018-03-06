
# newhybrids

This is the source code repository for the program *newhybrids* written by Eric C.
Anderson.  The old documentation which I have
put into `new_hybs_doc1_1Beta3.pdf` here is mostly relevant.

I have only tested this on Macs.

## Binaries

In `bin/PC` you will find `newhybrids.exe` which I compiled on a virtual PC 
(Windows 7) on my Mac, using `Compile-with-no-gui.sh`.  You can download it and then
do:

```sh
newhybrids.exe --help-full
```
for a full explanation of the command line syntax.

In `bin/OSX` you will find Mac OSX compatible binaries,
`newhybs` and `newhybsng`, which are the versions 
with and without the GUI, respectively.  You can download the whole repository
from GitHub, then if you are on the Terminal in the top level of the repository you can do
like this:
```sh
./bin/OSX/newhybs --help-full
```
for a full listing of options.  If this seems mysterious, find a friend who knows their
way around the Unix/Linux command line.


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

## About the version with the graphical user interface

If you are using the Mac OS X version with the GUI (or you have 
compiled it up that way on Linux), after the program launches and
reads the data, it opens up an information window and then it waits for
you to hit the space bar to start the program running.  At that point,
if you have installed the `NewHybrids_PreDefdViews.txt` file as described 
below, you should be able to hit the "1" key on your keyboard and have a bunch of
different windows open up to watch the progression of the MCMC.

## Getting predefined view to work
Until  I  write  a little installer package and change a few more
things this is how you are going  to  have  to  manually  install
NewHybrids Predefined Views for OS X:

After  cloning and building the repository you can move the file
`NewHybrids_PreDefdViews.txt` to the directory `~/Library/gfmcmc`,  
creating the directory if it is not already there.  In the above
the `~` stands for your home directory.


## Running a test data set
To see if you have it up and running on the test data set try issuing this 
command from the `newhybrids` directory:
```sh
./newhybs -d test_data/TestDat.txt
```
Then try hitting `1` key and then the `space-bar`.
