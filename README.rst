*****************************
 2d_hydro Build Instructions
*****************************

1) Clone using ssh:

    $ git clone git@github.com:zbyerly/2d_hydro.git

   or using https:

    $ git clone https://github.com/zbyerly/2d_hydro.git

   or download zipped version from https://github.com/zbyerly/2d_hydro and unzip

2) Create a separate build directory (prefered):
    
    $ mkdir my_2d_hydro_build
    $ cd my_2d_hydro_build

3) Compile with ifort:

    $ ifort -openmp -o 2d_hydro_executable ../2d_hydro/*.f90

   openmp flag can be left out to disable openmp support

4) Set openmp environment variables:

    $ export OMP_NUM_THREADS=8

5) create data directory inside run directory

    $ mkdir data

6) 2d_hydro_executable takes 3 arguments, the first determines which set of momenta is evolved (cart or cyl), the second determines the reconstruction scheme(ppm or minmod), and the third sets the rotation of the grid (rotate norotate). For example, to run the code with cartesian momenta, using minmod reconstruction, on a rotating grid:

    $ ./2d_hydro_executable cart minmod rotate

7) ???

8) Profit!

