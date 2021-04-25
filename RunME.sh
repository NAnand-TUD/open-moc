#!/usr/bin/env sh

HOME=$(pwd)

MOC=$HOME'/bin/'

echo "\n\n******************************************************\n"
echo "*********** SUPERSONIC CASCADE DESIGN TOOL ***********\n"
echo "******************************************************\n\n"
echo "Add the following in your ~/.bashrc\n"
echo "\texport MOC_HOME="$HOME
echo "\texport BIN="$MOC
echo "\tPATH=\$PATH:\$BIN\n\n"
echo "******************************************************\n\n"
