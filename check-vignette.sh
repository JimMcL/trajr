#!/bin/sh
#
# Searches for, and reports, any exported functions which aren't mentioned in trajr-vignette.Rmd

grep export NAMESPACE | sed 's/export(//; s/)//' | tr -d "\r" | while read fn
do
    n=$(grep "$fn" vignettes/trajr-vignette.Rmd | wc -l)
    [ "$n" -eq 0 ] && echo $fn
done
