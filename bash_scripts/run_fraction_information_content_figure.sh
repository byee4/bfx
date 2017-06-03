#!/bin/bash

python /home/bay001/projects/codebase/bfx/pyscripts/clipseq/fraction_information_content_figure.py \
--broad-table /home/bay001/projects/encode/permanent_data/fraction_information_content_figure/ALLsubmitted.output20160908.broadtable.csv \
--color-list /projects/ps-yeolab3/bay001/maps/current/fraction_information_content_figure/color_list_269.lines \
--order-file /home/bay001/projects/encode/permanent_data/fraction_information_content_figure/correlation_order.csv \
--low-cutoff 0.01 \
--output fraction.svg
