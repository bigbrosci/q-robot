#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 21:53:43 2019
@author: qli
"""
import urllib, http.cookiejar
from bs4 import BeautifulSoup

import scholarly
import sys, re
import pdftotext
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By

# Step1: Get the DOI from PDF file 
pdf_in = sys.argv[1]
#pdf_in = 'li2017_acscatal.pdf'
doi_re = re.compile("10.(\d)+/([^(\s\>\"\<)])+")
with open(pdf_in, "rb") as f:
    pdf = pdftotext.PDF(f)
doi = doi_re.search(pdf[0]).group(0)
print(doi)

# Step2: get the 302 redirect link from http://doi.org/XXX link 
cj = http.cookiejar.CookieJar()
starturl = 'https://doi.org/' + doi
opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj))
#f = urllib.request.urlopen(starturl)
f = opener.open(starturl)
url2 = f.geturl()
print(url2)
#
## Step3: Get the final link of the paper in the journal and read the webpage 
wd = webdriver.Firefox(executable_path=r'$HOME/bin/Q_robot/actions/geckodriver')
wd.get(url2)
if 'elsevier' in url2:
    wait = WebDriverWait(wd, 10).until(EC.url_changes(url2))
new_url = wd.current_url
html_page = wd.page_source
wd.quit()
print(new_url)
#
# Step4: get the title of this paper 
soup = BeautifulSoup(html_page, 'html.parser')
title = soup.title.string
paper_title = re.sub(' - .*$', '', title)
key_words = doi + ' ' + paper_title
print(key_words)

# Step5 USe the DOI and Title as key words to search the google scholar and  get the bib file link 
try: 
    search_query = scholarly.search_pubs_query(key_words)
    search_out = next(search_query)
    site = search_out.url_scholarbib
    print(site)
    
    ## Step6 get the bib file 
    hdr = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.11 (KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11',
           'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
           'Accept-Charset': 'ISO-8859-1,utf-8;q=0.7,*;q=0.3',
           'Accept-Encoding': 'none',
           'Accept-Language': 'en-US,en;q=0.8',
           'Connection': 'keep-alive'}
    
    request = urllib.request.Request(site, headers=hdr)
    page = urllib.request.urlopen(request)
    bib_infor = page.read().decode("utf8")
    print(bib_infor)
    
    # Step7 Save the bib file 
    out_name = pdf_in.replace('.pdf', '') + '.bib'
    file_out = open(out_name, 'w')
    file_out.writelines(bib_infor)
    file_out.close()
    
except Exception as e:
    '''If the paper is verrrry new, for example, publish today, we can not find the bib from google scholar '''
    print('An error occus: %s' %e)    
