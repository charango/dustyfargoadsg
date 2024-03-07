# Installation, compilation and first run

The code’s github repository contains 4 sub-directories : `docs` (for the documentation files), `src` (where the source files are), `in` (where the input parameter files are), and `lib` with an archive file to install the FFTW 2.1.5 library required for simulations with gas self-gravity. This library needs not be installed if you don’t plan on using gas self-gravity. If you do, enter the `lib` sub-directory and install the library as follows, which assumes that MPI is already installed on your environment :

```bash
tar zxvf fftw-2.1.5.tar.gz
cd fftw-2.1.5
./configure --enable-mpi --prefix=’/home/cbaruteau/dustyfargoadsg/lib/fftw2_1_5’
make
make install
```

where /home/cbaruteau/dustyfargoadsg/lib should be changed to the full path of the dustyfargoadsg/lib directory on your machine. Also, if you use a bash environment, you will need to add the following line to your .bashrc file in your home directory :

```bash
export FFTW_PREFIX=/home/cbaruteau/dustyfargoadsg/lib/fftw2_1_5
```

where again the full path to the dustyfargoadsg/lib/fftw2_1_5 directory should be indicated accordingly. Type `bash` in the command line of your terminal to account for the above changes.

To compile the code’s source files, go to the `src` directory. There are three compilation options depending on whether MPI and/or FFTW are installed on your machine :

```bash
make BUILD=parallelfftw # with MPI and FFTW
make BUILD=parallel     # with MPI but without FFTW
make BUILD=sequential   # without MPI nor FFTW
```

These compilation options (`BUILD=...`) need only be specified the first time source files are compiled. Any time after you may simply type `make`. Other compilation options can be set by editing the makefile in the `src` directory (options set by default are `-O3 -m64 -ffast-math`).

Next, go back to the main directory where the executable `fargo` has been written. Before running the code, you need to create the output directory where the simulations results are going to be written (contrary to FARGO3D, the present code does not automatically create the directory upon execution !). The output directory is set by `Outputdir` in the parameter file, so if you’re using one of the provided template input files in the `in` directory, simply type `mkdir out1`.

Now, to run your first simulation with Dusty Fargo-ADSG, simply type `./fargo -v in/template.par`. If you have compiled with MPI, you may type instead `mpirun -np 4 ./fargo -vm in/template.par` or any equivalent MPI execution command specific to your MPI environment. Over 4 cpus, it takes about 10 minutes to complete the run with the in/template.par parameter file (simply reduce the values of `Nrad` and `Nsec` to speed up the run). Outputs may be visualized with the public python program `fargo2python` located [here](https://github.com/charango/fargo2python).
