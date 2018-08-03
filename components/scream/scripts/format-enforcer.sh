#!/usr/bin/env sh

formatted=1

for f in $(find ./ -iname "*.cpp"); do
    clang-format ${f} > ${f}.fmt
    diff ${f} ${f}.fmt > /dev/null
    if [ $? -ne 0 ]; then
	echo "${f} is not formatted"
	formatted=0
    fi
done

if [ $formatted -eq 0 ]; then
    echo "Some files were not formatted, run the version of clang-format below"
    clang-format --version
    exit 1
fi
