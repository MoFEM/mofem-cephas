#!/bin/sh

#How to use:
#./name_change.sh "current" "new"
#This will search the files in this directory upwards for the word "current" and it will replace it with the word "new". Note: this can also be used for several words.
#In case these is a mistake, a copy of the original file is kept with the extension .org added to the end.
#If you are satisfied that the changes haven't broken anything and you no longer need the backup files then remove them using this command: 
#find . -name "*org" -exec rm -rf {} \;

find . -name "*.h" -exec sed -i .org "s/$1/$2/g" {} \;
find . -name "*.c" -exec sed -i .org "s/$1/$2/g" {} \;
find . -name "*.hpp" -exec sed -i .org "s/$1/$2/g" {} \;
find . -name "*.cpp" -exec sed -i .org "s/$1/$2/g" {} \;
find . -name "*.txt" -exec sed -i .org "s/$1/$2/g" {} \;
find . -name "*.jou" -exec sed -i .org "s/$1/$2/g" {} \;
find . -name "*.sh" -exec sed -i .org "s/$1/$2/g" {} \;
find . -name "*.cmake" -exec sed -i .org "s/$1/$2/g" {} \;
find . -name "README" -exec sed -i .org "s/$1/$2/g" {} \;


