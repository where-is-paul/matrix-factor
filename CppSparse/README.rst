===========================================================
CppSparse: A C++ sparse matrix library with Python bindings
===========================================================

About
-----
I have been fascinated with sparse direct solvers ever since I
encountered them whilst learning finite element analysis back during
my undergraduate days. One particularly clever trick that stood out
was the *multi-frontal* method developed by Bruce M. Irons
[Irons1970]_ back in the 1960s that performed Gaussian elimination on
the stiffness matrix *while it was being assembled* in the usual
element-by-element fashion.

As a graduate student, I was fortunate enough to re-learn finite
element methods from Prof. K-J Bathe. One of the intriguing aspects of
sparse solvers he mentioned in his class was the use of methods based
on graph theory that had surpassed classical *skyline* solvers such as
COLSOL described in his book [Bathe1995]_. This was the first time I
had heard about the connection between finite element methods and
graph theory - two very distinct worlds that are deeply connected at a
rather profound level. Although my curiosity was piqued, I did not dig
deeper into sparse solvers until 2006 when I noticed an announcement
in SIAM News about an upcoming book by Tim Davis on sparse direct
solvers [Davis2006]_ and immediately preordered it.

I started reading Davis' book in earnest during Fall 2006 and quickly
realized that the only way I was going to benefit from the book was to
digest the (rather terse) algorithms and reimplement them myself along
with tackling all the exercises. Unfortunately, this has been a rather
painstaking process and progress was glacially slow until 2011 when I
resolved to finally dedicate a part of my weekends to understand
sparse solvers, matrix reordering algorithms and the graph theoretical
foundations of sparse linear solvers.

*CppSparse* has been the result of my work so far. It is a sparse matrix
library based on the algorithms and exercises in Tim Davis' book. It
is implemented in C++ and extensively uses the containers and
algorithms in the standard library. This allows me to elide details
such as checking for out-of-memory conditions or managing the
lifetime of my temporary buffers. The templated nature of *CppSparse*
allows it to be used in both 32- and 64-bit libraries and the same
code base supports both real and complex-valued sparse matrices.

*CppSparse* is distributed as a self-contained header-only library. It
has very few dependencies besides a half-decent C++ compiler. I
primarily use *CppSparse* through its Python bindings and aim for it
to be a fully-featured sparse matrix library for Python.

Building CppSparse
------------------
*CppSparse* currently reliably builds on Windows using Microsoft
Visual Studio 2010 and MingW GCC 4.5 (64-bit) and is tested with
64-bit Python 2.7 and NumPy. I also occassionally build it on my Mac
OS Snow Leopard and back-port the changes to Windows.

Some of the routines in *CppSparse* use C++11 features such as lambda
expressions that are not supported out-of-the-box on the older GCC
installed on Snow Leopard; I will eventually remove these features
once I have completed a few more chapters in Davis' book. I also had
move constructors on the matrix classes that I disabled because SWIG
seemed to choke on them.

Prequisites
~~~~~~~~~~~
* `SWIG`_ version 2.0.1 or better for generating the Python bindings
* `Microsoft Windows 7 SDK (64-bit)`_ that contains the 64-bit version of the Visual C++ compiler, or, 
* `MingW64`_ with GCC 4.5 or better for compiling the generated bindings
* `Python`_ 64-bit version 2.7 or better (I have not tested *CppSparse* with Python 3.x)
* `NumPy`_ compatible with your Python installation. You can get a 64-bit installer for Windows from `Ch. Golke`_ 's website.
* `GraphViz`_ (optional). Some of the functions generate graphs in the GraphViz format; you can use tools like ``dot`` to post-process these graphs and create pretty pictures.
.. _SWIG: http://www.swig.org 
.. _MingW64: http://tdm-gcc.tdragon.net/
.. _Microsoft Windows 7 SDK (64-bit): http://www.microsoft.com/download/en/details.aspx?id=8279
.. _Python: http://www.python.org
.. _NumPy: http://www.numpy.org
.. _Ch. Golke: http://www.lfd.uci.edu/~gohlke/pythonlibs
.. _GraphViz: http://www.graphviz.org

Building with the Microsoft C++ compiler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Ensure that SWIG is in your path.
* Start up a 64-bit Visual C++ command prompt, navigate to the source
  directory and run ``buildswig.bat``. 

Building with MingW GCC
~~~~~~~~~~~~~~~~~~~~~~~
* Ensure that SWIG and GCC are in your path.
* Navigate to the source directory and run ``make -f Makefile.mingw``

Future plans
------------
* Finish Cholesky factorization
* Work on QR factorization
* Work on LU factorization
* Work on reordering techniques (phew!)
* Make all return types be ``std::shared_ptr`` so they don't get
  copied constantly in SWIG.
* Add more PyUnit tests
* Finish matrix arithmetic operations
* ...


License
-------
*CppSparse* is licensed under the `MIT/X11 license`_:

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

.. _`MIT/X11 license`: http://www.opensource.org/licenses/mit-license.php


References
----------
.. [Bathe1995] Bathe, K-J. *Finite Element Procedures*, Prentice Hall, 1995. 
.. [Davis2006] Davis, T. *Direct Methods for Sparse Linear Systems*, SIAM, Philadelphia, PA, 2006.
.. [Irons1970] Irons, B.M. *A frontal solution scheme for finite element analysis*, IJNME, 2(5--32), 1970.

