#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 21:53:43 2019
@author: qli
"""
import sys, re
import pdftotext

# Step1: Get the DOI from PDF file 
pdf_in = sys.argv[1]

def get_doi(page):
    doi_re = re.compile("10.(\d)+/([^(\s\>\"\<)])+")
    doi = doi_re.search(page).group(0)
    return doi
#doi = get_doi(pdf[0])

with open(pdf_in, "rb") as f:
    pdf = pdftotext.PDF(f)
print(pdf[0])


