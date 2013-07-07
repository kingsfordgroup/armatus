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
    $ cd build
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

* armatus executable: (Boost Program Options): Armatus.cpp

ArmatusUtil.{cpp,hpp}:

* parseMatrix(File): 3C matrix parser (Dixon et al. format)
* outputDomains(set of domains)
* consensus(set of domains)

Data structures {cpp,hpp}:

* ArmatusParams(Matrix, gamma): parameters to DP and pre-computable quantities such as the quality function
* ArmatusDAG(Params): encodes structure of dynamic program
    * BackPointers and SubProblem classes: see below
    * build(): Build the DAG for the dynamic program (Boost Graph)
    * computeTopK(k): At every node, store 'SubProblem': k 3-tuples: (edge, child solution, score of kth best solution for this subproblem)
    * extractTopK(k): Returns a set of domains
* WeightedIntervalScheduling
    * WeighedInterval
    * WeightedIntervalScheduler
        * previousDisjointInterval()
        * computeScheduling()
        * extractIntervals()

Testing:

* Boost test

Pipeline:

* Code and data links to generate results in application note