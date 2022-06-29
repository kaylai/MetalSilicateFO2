About the Tests Folder and Files
================================

Unit testing
------------
python unittest library is implemented for tests of the code. This is used to
ensure code is behaving as expected. Benchmarks were made early on by comparing
values computed with this python code to those reported in the literature and
computer for the same data using other tools (see "Benchmarking" section).

Run unit tests as:
python3 -m unittest

Benchmarking
------------
Activity coefficient values were computed for samples published in Corgne et al.
(2008) "Metal–silicate partitioning and constraints on core composition and
oxygen fugacity during Earth accretion". Values computed with this python code
were compared to those reported in Corgne et al. (2008) and those computed
using the Norris Scientific Metal Activity calculator available online
http://norris.org.au/expet/metalact/.

