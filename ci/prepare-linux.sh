#!/bin/bash

if [ $(uname -m) == 'x86_64' ];
then
    TAG='x86_64';
else
    TAG='i386';
fi

curl -O https://cmake.org/files/v3.8/cmake-3.8.2-Linux-$TAG.sh
bash cmake-3.8.2-Linux-$TAG.sh --skip-license --prefix=/

hash -r
which cmake
cmake --version
