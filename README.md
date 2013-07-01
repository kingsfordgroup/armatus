armatus
=======

Installation
------------

Dependencies:

* C++11
* Boost (Graph, Test, ublas, program options, system)

Note: Although this should be conceptually easy, properly installing these dependencies can be tricky since on a system like OS X you can choose to use different compilers and standard libraries.  Please make sure that boost is installed with the same standard libary and compiler that you are using to compile Armatus. For this reason, binaries of Armatus will also be released.

Compilation:

    $ mkdir build
    $ cmake [options] ..

Here is an example of how to compile boost with Clang++ and use Clang++ to compile Armatus on OS X:

Boost-specific details:

    ./bootstrap --prefix=$HOME/boost
    ./b2 clean
    ./b2 install toolset=clang cxxflags="-stdlib=libc++" linkflags="-stdlib=libc++"

CMake:

    cmake -DCMAKE_CXX_COMPILER=clang++ -DBOOST_ROOT=$HOME/boost -DBoost_NO_SYSTEM_PATHS=true ..
    make
    export DYLD_FALLBACK_LIBRARY_PATH=$HOME/boost/lib

Components
----------

Core:

* armatus executable: (Boost Program Options)
* Parse3CMatrix(File): 3C matrix parser (Dixon et al. format)
* BuildDAG(Matrix): Build the DAG for the dynamic program (Boost Graph)
* Viterbi(DAG): Run the Viterbi algorithm on DAG to solve the dynamic program
* MultipleSolutions(DAG): Contains algorithms for multiple solutions

Data structures:

* SolutionDAG
* 

Testing:

* Boost test

Pipeline:

* Code and data links to generate results in application note