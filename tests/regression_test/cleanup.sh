#!/bin/bash

# This script cleans up output and log files created 
# by regression tests, and restores the status before
# the tests.
#
# Usage:
#    $ ./cleanup.sh

rm test_spdyn/*/*/error*     >& /dev/null
rm test_spdyn/*/*/test*      >& /dev/null
rm -r test_spdyn/*/*/output.?*      >& /dev/null

