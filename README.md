
# newhybrids

This is the source code repository for the program *newhybrids* written by Eric C.
Anderson.  The old documentation which I have
put into `new_hybs_doc1_1Beta3.pdf` here is mostly relevant.

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
