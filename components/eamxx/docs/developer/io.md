# Input-Output

In EAMxx, I/O is handled through the SCORPIO library, currently a submodule of E3SM.
The `eamxx_io` library within eamxx allows to peform read/write operations directly
into/from eamxx data structures. If read/write need to be performed with raw data types
the `eamxx_scorpio_interface` library provides more low-level interfaces to the
SCORPIO library.
