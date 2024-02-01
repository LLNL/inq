#!/bin/bash
# Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

set -e 
set -x

tmpdir=`mktemp -d`
cd $tmpdir
sed "s|^inq|$INQ_EXEC_ENV inq|g" $1 > script.sh
chmod +x script.sh
shift
./script.sh $*
ret=$?
cd -
rm -fr $tmpdir
exit $ret
