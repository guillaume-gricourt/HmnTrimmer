#!/usr/bin/env python
# coding: utf8

import argparse
import base64
import codecs
import datetime
import django
import io
import json
import logging
import os
import sys
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import seaborn as sns

from collections import namedtuple
from datetime import date
from django.conf import settings
from django.template.loader import render_to_string
from packaging import version

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%d-%m-%Y %H:%M')

##############
## Variable ##
##############

VERSIONUSED = {
    "min" : "0.6.0",
    "max" : "0.7.0"
}

##############
## Function ##
##############

def valid_files(choices, fname):
    isValid = True
    ext = os.path.splitext(fname)
    if len(ext) > 0:
        ext = ext[1][1:]
        if not ext in choices:
            isValid = False
    else:
        isValid = False
    if not isValid:
        parser.error("Filename must be finished by one of {}".format(choices))
    return fname

def readJson(fpath):
    with open(fpath) as fid:
        dico = json.load(fid)
    return dico

def formatNb(nb):
    return '{:,}'.format(round(float(nb))).replace(","," ").split(".")[0]

def plt2base64(plt):
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek( 0 )
    base = base64.b64encode( buf.read() )
    buf.close()
    return base

def graph_statistics_number(total, kept):
    '''
    Arg : total reads, kept reads
    '''
    plt.close()
    colors, labels = ['#ff7f0e', '#9467bd'], ['total', 'discarded']
    total, kept = int(total), int(kept)
    values = [total, total-kept]

    fig, ax = plt.subplots(figsize=(4,4))
    ax.pie(values, wedgeprops=dict(width=0.5, edgecolor='w'), colors=colors)# startangle=-30, colors=colors)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), labels=labels, frameon=False, borderpad=0, borderaxespad=0)

    fig.subplots_adjust(left=0,right=0.7,bottom=0,top=1)
    return plt2base64(plt).decode()

def graph_statistics_distribution(data):
    '''
    Arg : dict, key : length of read, value, number of read
    '''
    plt.close()

    x, y = [], []
    for read_length, number in sorted(data.items(), key=lambda x: int(x[0])):
        read_length, number = int(read_length), int(number)
        x.append(read_length)
        y.append(number)
    #Graph
    fig, ax = plt.subplots()
    ax.bar(x, y)
    #Trim border
    sns.despine(ax=ax)
    #Annotate
    ax.set_xlabel('Taille (pb)', {'weight':'bold'})
    ax.set_xlim(left=0, right=x[-1]+10)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), labels=['Mean length : %.0f'%(round(sum(x)/len(x)),)], frameon=False)
    plt.tight_layout()
    return plt2base64(plt).decode()

def write(foutput, template):
    with codecs.open(foutput, "w", "utf-8") as fod:
        fod.write(template)

##############
## Argparse ##
##############
description = 'Rendering report file created by HmnTrimmer'
parser = argparse.ArgumentParser(description=description)
#Input
parser.add_argument('-i', '--input-file', type=lambda x: valid_files(['json'], x), required=True, help='Report file json')
parser.add_argument('-t', '--template-file', type=lambda x: valid_files(['html'], x), required=True, help='Django template file used for rendering')
parser.add_argument('-o', '--output-file', type=lambda x: valid_files(['html'], x), required=True, help='Output file rendered')

logging.info('Checking arguments')
if len(sys.argv) <= 1:
    parser.print_help()
    sys.exit(1)
else:
    args = parser.parse_args()

finput = args.input_file
ftemplate_file = args.template_file
foutput = args.output_file

#################
## Load Django ##
#################
settings.configure(DEBUG=True,
    TEMPLATES = [{
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [os.path.dirname(ftemplate_file)]
    }],
)
django.setup()

######################
##  Read file input ##
######################
logging.info('Read input file')
data = readJson(finput)
context = {}

if version.parse(data.get('software',{}).get('version', '0')) < version.parse(VERSIONUSED['min']) or \
    version.parse(data.get('software',{}).get('version', '0')) > version.parse(VERSIONUSED['max']):

    parser.error('Version of report file supported between %s and %s'%(VERSIONUSED['min'], VERSIONUSED['max']))

######################################
##  Construct template information  ##
######################################

logging.info('Construct information')
# Software informations.
SoftwareTuple = namedtuple('SoftwareTuple', ['name', 'version'])
software = SoftwareTuple(
    name=data.get('software',{}).get('name', 'undefined'),
    version=data.get('software',{}).get('version', 'undefined'))
context['software'] = software

# Args informations
def format_filename(filenames):
	filenames = [ os.path.basename(x) for x in filenames ]
	return ', '.join(filenames)

FileTuple = namedtuple('FileTuple', ['finput', 'foutput'])
ffile = FileTuple(
    finput = format_filename(data.get('analyze', {}).get('file',{}).get('input',[])),
    foutput = format_filename(data.get('analyze', {}).get('file',{}).get('output',[])))
context['file'] = ffile

TrimmerTuple = namedtuple('TrimmerTuple', ['name', 'value'])
trimmers = []
for trimmer, value in data.get('analyze',{}).get('trimmers',{}).items():
    trimmers.append(TrimmerTuple(
        name = trimmer,
        value = value))
context['trimmers'] = trimmers

# Statistics number.
StatisticsNumberTuple = namedtuple('StatisticsNumberTuple', ['total', 'kept', 'discarded', 'graph'])
total = data.get('statistics', {}).get('total', 0)
kept = data.get('statistics', {}).get('kept', 0)
discarded = data.get('statistics', {}).get('discarded', 0)
graph = graph_statistics_number(total, kept)
statistics_number = StatisticsNumberTuple(
    total = formatNb(total),
    kept = formatNb(kept),
    discarded = formatNb(discarded),
    graph = graph)
context['statistics_number'] = statistics_number

# Statistics distribution.
StatisticsDistributionTuple = namedtuple('StatisticsDistributionTuple', ['before', 'after'])
before = graph_statistics_distribution(data.get('statistics', {}).get('length_reads_before', {}))
after = graph_statistics_distribution(data.get('statistics', {}).get('length_reads_after', {}))
statistics_distribution = StatisticsDistributionTuple(
    before = before,
    after = after)
context['statistics_distribution'] = statistics_distribution

##############################
##  Rendering information   ##
##############################

logging.info('Rendering information')
template = render_to_string(os.path.basename(ftemplate_file), context)

##########################
##  Write output file   ##
##########################
logging.info('Write output file')
write(foutput, template)
