==========================
MPAS-Seaice testing system
==========================

The MPAS-Seaice testing system consists of a series of python scripts that are
used to test the MPAS-Seaice model.

Usage
=====

.. code::

   > test_mpas-seaice.py [-h] -d MPASDEVELOPMENTDIR [-b MPASBASEDIR]
                         [-t TESTSUITE] [-o DOMAINSDIR] [-a] [-c]

Options
=======

+-----------------+----------------------------------------------------------------+
| Name            | Description                                                    |
+=================+================================================================+
| -h              | Display the help screen which describes the available options. |
+-----------------+----------------------------------------------------------------+
| -d, --dev       | Specify a path to the MPAS checkout to be tested.              |
+-----------------+----------------------------------------------------------------+
| -b, --base      | [optional]: Specify a path to the MPAS checkout to be.         |
|                 | tested against, in the case of regression tests. If this       |
|                 | option is not specified, only tests that do not require        |
|                 | a base MPAS checkout to check against will be performed.       |
+-----------------+----------------------------------------------------------------+
| -t, --testsuite | [optional, default: /testsuites/testsuite.standard.xml]:       |
|                 | Specify the testsuite for the testing system to test           |
|                 | with.                                                          |
+-----------------+----------------------------------------------------------------+
| -o, --domaindir | [optional, default: env variable                               |
|                 | MPAS_SEAICE_DOMAINS_DIR]: This specifies the domains           |
|                 | directory for the system to use to get domain files to         |
|                 | build testing cases with.                                      |
+-----------------+----------------------------------------------------------------+
| -a, --avail     | [optional]: List the tests that have been implemented.         |
+-----------------+----------------------------------------------------------------+
| -c, --check     | [optional]: This option set namelist options which will        |
|                 | cause MPAS-Seaice to fail all the tests. This is for           |
|                 | testing the testing system.                                    |
+-----------------+----------------------------------------------------------------+

Testsuite .xml files
====================

The tests that will be performed are defined in a testsuite .xml file, stored in
the testsuites directory in the testing system. Below is an example testsuite
file that specifies the configurations (namelist and streams file), domains
(grid, forcing and graph files) and tests that should be performed:

.. code::

   <?xml version="1.0" encoding="UTF-8"?>
   <testsuite name="standard">
       <configuration name="standard_physics">
           <domain name="domain_QU120km">
               <test name="regression"/>
               <test name="parallelism"/>
               <test name="restartability"/>
           </domain>
       </configuration>
   </testsuite>

This example specifies a series of tests (**regression**, **parallelism** and
**restartability**) to be performed with the **domain_QU120km** domain and the
**standard_physics** configuration. Test suites are stored in
``MPAS/testing_and_setup/seaice/testing`` within the MPAS repo.

Configurations
==============

The testing system uses configurations (i.e. namelist and stream file pairs)
listed in the testsuite xml file in the test runs that it performs. The
files are stored in the repository at ``MPAS/testing_and_setup/seaice/configurations``.

Domains
=======

The testing system uses domains listed in the testsuite xml file. The domain
provides grid, forcing and graph files for the testing runs. The system being
used must have a directory containg various domain directories. These
directories must contain ``get_domain.py`` script files to allow the testing
system to set up runs.

Available Tests
===============

1. **Regression**: Bit reproducibility is tested between the dev and base MPAS checkouts.

2. **Parallelism**: Bit reproducibility is tested between different processor counts for the same
dev MPAS checkout.

3. **Restartability**: Bit reproducibility is tested between a standard run and a run with a restart
half way through.


Required python packages
========================

Required python packages can be installed with Conda:

.. code::

   > conda create -n mpas_seaice python=3.6
   > conda install -n mpas_seaice netCDF4
   > conda install -n mpas_seaice -c conda-forge f90nml
   > conda install -n mpas_seaice -c anaconda colorama
