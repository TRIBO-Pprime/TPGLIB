#!/bin/sh

cd    STATS ; make all
cd ../FILTR ; make all
cd ../DERIV ; make all
cd ../ASFC2 ; make all
cd ../MORPH ; make all
cd ../ABBOT ; make all
cd ../ANISO ; make all

echo "OK"
read wait

