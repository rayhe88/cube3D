![GitHub tag (release)](https://img.shields.io/github/v/release/rayhe88/cube3D?label=version)
![GitHub last commit](https://img.shields.io/github/last-commit/rayhe88/cube3D?label=last%20modified)
![GitHub top language](https://img.shields.io/github/languages/top/rayhe88/cube3D?color=green)

---

             ______      __        _____ ____
            / ____/_  __/ /_  ___ |__  // __ \
           / /   / / / / __ \/ _ \ /_ </ / / /
          / /___/ /_/ / /_/ /  __/__/ / /_/ /
          \____/\__,_/_.___/\___/____/_____/

                 README file for Cube3D
                     February, 2024

---

## What is Cube3D?

Cube3d is a code part of the GPUAM project designed to search for critical
points of scalar fields that are used in computational/quantum chemistry
such as electron density. This search is carried out completely numerically
and takes as input a file in the Gaussian "cube" format.

This code has been developed at the UAM-I (Metropolitan Autonomous University

- Iztapalapa) in the area of Theoretical Physical Chemistry and updated during
  my postdoctoral stays at McMaster University and Argonne National Lab.

Perform a completely numerical search using 3-dimensional Lagrange interpolators.

## How to install?

To install it, download the latest version from the GitHub repository.
Go to the `src` subdirectory, inside you will find a makefile file,
which has all the information necessary to compile it.
This is a code that is parallelized at the OMP level, so it is beneficial to have the OMP library.

Inside `src` type make

```bash
$ make
gcc -c -Wunused  -Wno-unused-result -O2 main.c
gcc -c -Wunused  -Wno-unused-result -O2 -DD file.c
gcc -c -Wunused  -Wno-unused-result -O2 array.c
gcc -c -Wunused  -Wno-unused-result -O2 graph.c
gcc -c -Wunused  -Wno-unused-result -O2 fields.c -fopenmp
gcc -c -Wunused  -Wno-unused-result -O2 transU.c
gcc -c -Wunused  -Wno-unused-result -O2 tableP.c
gcc -c -Wunused  -Wno-unused-result -O2 timing.c
gcc -c -Wunused  -Wno-unused-result -O2 jacobi.c
gcc -c -Wunused  -Wno-unused-result -O2 -DGRAPH  kernels.c
gcc -c -Wunused  -Wno-unused-result -O2 lagForm.c
gcc -c -Wunused  -Wno-unused-result -O2 lebedev.c
gcc -c -Wunused  -Wno-unused-result -O2 lectura.c
gcc -c -Wunused  -Wno-unused-result -O2 findCrit.c -fopenmp
gcc -c -Wunused  -Wno-unused-result -O2 geomData.c
gcc -c -Wunused  -Wno-unused-result -O2 rotation.c
gcc -c -Wunused  -Wno-unused-result -O2 cubeIndex.c
gcc -c -Wunused  -Wno-unused-result -O2 lagrange2.c
gcc -c -Wunused  -Wno-unused-result -O2 replicate.c
gcc -c -Wunused  -Wno-unused-result -O2 mathTools.c
gcc -c -Wunused  -Wno-unused-result -O2 numBondPath.c -fopenmp
gcc -c -Wunused  -Wno-unused-result -O2 runCommands.c -fopenmp
gcc -c -Wunused  -Wno-unused-result -O2 utils.c
gcc -o cube3D.x main.o file.o array.o graph.o fields.o transU.o tableP.o timing.o jacobi.o kernels.o lagForm.o lebedev.o lectura.o findCrit.o geomData.o rotation.o cubeIndex.o lagrange2.o replicate.o mathTools.o refinement.o numBondPath.o runCommands.o utils.o -lm  -fopenmp
```

and finally, install by the next command:

```bash
$ make install
```

The binary is located in `/home/USER/bin/cube3d.x`.

Make sure to export this subdirectory in your `PATH` so you can invoke it from any location.

```bash
$ export PATH= $PATH:/home/USER/bin/cube3d.x
```

## First steps

Once the binary is installed and the location of the binary is exported in the path, you can go to `proof` subdir. There is an example cube file, called `ureaII.cube` which corresponds to a urea crystal. To generate an input file type the following.

```bash
cube3d.x -f ureaII.cube -i input.in
```

This will create a cube3d input file, the input file contains the information to run cube3d, the name of the cube file, the tasks it can evaluate, and various extra options. Take a look at that file.
