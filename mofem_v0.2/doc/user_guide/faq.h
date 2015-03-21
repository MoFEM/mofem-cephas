/*! \page faq_and_bugs FAQ and Bugs

\section bugs Bugs

- BasicMoFEMEntity in should be linked to directly to MoAB data structures such
  that connectivity and nodal coordinates could be quickly accessed, without
  need of using native MoAB functions.

\section faq Frequently Asked Questions

\subsection update_on_memory_stick Ho to update MoFEM on Live USB Stick?

MoFEM update on Live USB stick:
\code
$ mofem_update.sh
$ mofem_build.sh
\endcode

Following command run test verifying updated code:
\code
$ mofem_fast_check.sh
\endcode

If you run MoFEM update at University of Glasgow behind proxy server, set proxy
servers as follows:
\code
$ export http_proxy=http://wwwcache.gla.ac.uk:8080
$ export https_proxy=http://wwwcache.gla.ac.uk:8080
\endcode

\subsection ctest How to run ctest?

You can run tests and report results to MoFEM CDash web page. Form mofem user
modules build directory executing  run script
\code./bin/mofem_fast_check.sh\endcode
Results of test can be seen on <http://cdash.eng.gla.ac.uk/cdash/>. 

Note that test tests for MoFEM library and \em User \em Modules are run independently
and can be seen as a two different projects.

If you run test behind proxy server you can have to set \em http_proxy and \em
http_proxy environmental variables. For example if you run mofem at Glasgow
University, please do:
\code
$ export http_proxy=http://wwwcache.gla.ac.uk:8080
$ export https_proxy=http://wwwcache.gla.ac.uk:8080
\endcode

You can as well run ctest directly by simply executing command line:
\code 
ctest -V -D Experimental
\endcode
where option -V sets verbose version and all test output is printed on screen
and -D Experimental tels ctest to submit results to Experimental build on CDash
MoFEM server.

*/


