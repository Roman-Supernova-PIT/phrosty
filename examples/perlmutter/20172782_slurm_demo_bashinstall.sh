#!/bin/bash 

cd /home/phrosty 
pip install -e . 
cd /home/sfft 
pip install -e . 
cd /home/phrosty 
SNPIT_CONFIG=/home/phrosty/examples/perlmutter/phrosty_config.yaml \
python phrosty/pipeline.py \
      -oc ou2024 \
      -ic ou2024 \
      --oid 20172782 \
      -r 7.551093401915147 \
      -d -44.80718106491529 \
      -b R062 \
      -t /home/phrosty/examples/perlmutter/20172782_instances_templates_1.csv \
      -s /home/phrosty/examples/perlmutter/20172782_instances_science.csv \
      -p 9 \
      -w 5 \
      -v