#!/bin/bash

# if [ ! -f ".clang-format" ]; then

# 	echo "Generate clang format configuration file"
# 	clang-format -style=Google -dump-config > .clang-format 
#     sed -i 's/SortIncludes:    true/SortIncludes:    false/g' .clang-format
# fi

echo "Formating header files"
find ./models/ -name "*.hpp" -exec clang-format  -i "{}" ";"

echo "Formating source files"
find ./models/ -name "*.cpp" -exec clang-format -i "{}" ";"
