#!/usr/bin/env bash


for ID in `cat $1`; do
  IDSTR=${ID}":"
  OUTCOME=$(cat $@ | grep ^TR2: | grep -w $IDSTR | cut -d':' -f4)
  FILTER=
  if [[ "${OUTCOME}" =~ "intron did not pass the filter" ]]; then
    FILTER=$(cat $@ | grep -B3 ^TR2: | grep -B3 -w $IDSTR | head -1)
  elif [[ "${OUTCOME}" =~ "as intron chain exists" ]]; then
    OUTCOME2=" "$(cat $@ | grep -B1 ^TR: | grep -B1 -w $IDSTR | head -1 | cut -d'-' -f2)
    if [[ "${OUTCOME2}" =~ "Skipping transcript as partially redundant" ]]; then
      OUTCOME=$OUTCOME2
    fi
  fi    
  echo -e $ID"\t"$OUTCOME"\t"$FILTER
done
