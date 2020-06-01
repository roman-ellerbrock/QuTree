#!/bin/bash
source scl_source enable devtoolset-7
source scl_source enable rh-git218
exec "$@"