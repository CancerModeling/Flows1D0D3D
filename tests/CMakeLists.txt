#////////////////////////////////////////////////////////////////////////////////
#//  Copyright (c) 2020-2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
#//
#//  Distributed under the Boost Software License, Version 1.0. (See accompanying
#//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#////////////////////////////////////////////////////////////////////////////////

# ****************************************************************************
# Model tests
# ****************************************************************************
add_test(NAME two_vessels
        COMMAND ${BASH_PROGRAM} ./run.sh
        WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/../tests/test_data/NetFVFE/two_vessels
        )
