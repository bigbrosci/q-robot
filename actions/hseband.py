#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from eigenval import *

## Read the EIGENCAL file from HSE calculations and save the weight 0 parts.
save_eigenval_hse()

write_band_energy('EIGENVAL-HSE')
