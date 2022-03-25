#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 08:42:52 2020

@author: qli
"""

import urllib, http.cookiejar
import sys, re

smiles = sys.argv[1]

# Step2: get the 302 redirect link from http://doi.org/XXX link 
cj = http.cookiejar.CookieJar()
starturl = 'https://cactus.nci.nih.gov/chemical/structure/' + smiles + '/smiles'
opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj))
#f = urllib.request.urlopen(starturl)
f = opener.open(starturl)
url2 = f.geturl()

hdr = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.11 (KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11',
       'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
       'Accept-Charset': 'ISO-8859-1,utf-8;q=0.7,*;q=0.3',
       'Accept-Encoding': 'none',
       'Accept-Language': 'en-US,en;q=0.8',
       'Connection': 'keep-alive'}
    
request = urllib.request.Request(url2, headers=hdr)
page = urllib.request.urlopen(request)
bib_infor = page.read().decode("utf8")
print(bib_infor)