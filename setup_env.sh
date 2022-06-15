#!/bin/bash


#voms-proxy-init -voms cms
voms-proxy-init --rfc --voms cms -valid 192:00
cp /tmp/x509up_u124616  /afs/cern.ch/user/c/cbasile/
export X509_USER_PROXY=/afs/cern.ch/user/c/cbasile/x509up_u124616
#export PATH=/afs/cern.ch/work/c/cbasile/bin:$PATH
