#### buildconfig.sh : 
#### Local script to rebuild the config and make files using the aclocal and auto* 
#### family of configuration tools. Ordinary users should not need to run this script. 
#### Author:  Maxie D. Schmidt (maxieds@gmail.com)
#### Created: 2018.06.17

#!/bin/bash

aclocal
autoheader
automake
automake --add-missing
autoconf

