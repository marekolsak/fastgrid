#!/bin/bash
find \( -name "*.cpp" -o -name "*.h" \) -exec grep -q -e "$1" {} \; -exec echo "In file {}:" \; -exec grep -e "$1" {} \; -exec echo "" \;

