#!/bin/sh

find . -name "*.h" -exec sed -i .org "s/$1/$2/g" {} \;
find . -name "*.c" -exec sed -i .org "s/$1/$2/g" {} \;
find . -name "*.hpp" -exec sed -i .org "s/$1/$2/g" {} \;
find . -name "*.cpp" -exec sed -i .org "s/$1/$2/g" {} \;
find . -name "*.txt" -exec sed -i .org "s/$1/$2/g" {} \;
find . -name "*.jou" -exec sed -i .org "s/$1/$2/g" {} \;
find . -name "README" -exec sed -i .org "s/$1/$2/g" {} \;


