#!/usr/bin/env python

# for tgi cluster:
#/gapp/x64linux/opt/pythonbrew/venvs/Python-2.7.6/gemini/bin/python
# for uva cluster:

import pysam
import json
import sys
import argparse
from argparse import RawTextHelpFormatter
import string
from string import *

__author__ = "Colby Chiang (colbychiang@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2017-04-12 01:36 $"


# get the number of entries in the set
def countRecords(myCounter):
    numRecords = sum(myCounter.values())
    return numRecords

# sum of the entries
def sumRecords(myCounter):
    mySum = 0.0
    for c in myCounter:
        mySum += c * float(myCounter[c])
    return mySum

# calculate percentiles from the depth counter.
# uses the Nearest Rank definition of percentile
def percentile(myCounter, p):
    p = float(p)
    if p <= 0 or p >= 100:
        return
    
    #length is the number of bases we're looking at
    numEntries = countRecords(myCounter)
    # the ordinal value of the output element
    limit = int(round(p/100 * numEntries + 0.5))

    # a list of the values, sorted smallest to largest
    # note that this list contains unique elements only
    valueList = list(myCounter)
    valueList.sort()
    numValues = len(valueList)

    # the number of entries through the set we've gone
    runEntries = 0
    # iterator i
    i = 0
    # initiate v, in case one element in input
    v = valueList [i]

    # move through the value list, iterating by number of
    # entries for each value
    while runEntries < limit:
        v = valueList[i]
        runEntries += myCounter[v]
        i += 1
        # if i is greater than numValues, just return the largest
        if i == (numValues):
            v = valueList[i - 1]
            break
    return v

# calculate the arithmetic mean, given a counter and the
# length of the feature (chromosome or genome)
# for x percentile, x% of the elements in the set are
# <= the output value
def mean(myCounter):
    # the number of total entries in the set is the
    # sum of the occurrences for each value
    numRecords = countRecords(myCounter)

    # u holds the mean
    u = float()

    u = sumRecords(myCounter) / numRecords
    return u

def stdev(myCounter):
    # the number of total entries in the set is the
    # sum of the occurrences for each value
    numRecords = countRecords(myCounter)

    # u holds the mean
    u = mean(myCounter)
    sumVar = 0.0

    # stdev is sqrt(sum((x-u)^2)/#elements)
    for c in myCounter:
        sumVar += myCounter[c] * (c - u)**2
    myVariance = float(sumVar) / numRecords
    stdev = myVariance**(0.5)
    return stdev


# # calculate the density curve for and insert size histogram
# def calc_insert_density(self):
#     dens = Counter()
#     for i in list(self.hist):
#         dens[i] = float(self.hist[i])/countRecords(self.hist)
#     self.dens = dens
        
def bamfixflags(bamfile,
                lib_info_file,
                limit,
                is_sam,
                bam_out,
                uncompressed_out,
                outlier_bound = 2.0,
                mapping_bound = 3.0,
                max_stddev = 4.0):

    proper = {}

    lib_info = json.load(lib_info_file)
    for sample in lib_info:
        # print lib_info[sample]['libraryArray'][0]
        for lib in lib_info[sample]['libraryArray']:
            lib_hist = {int(k):int(v) for k,v in lib['histogram'].items()}
            p25 = percentile(lib_hist, 25)
            p75 = percentile(lib_hist, 75)
            lib_mean = mean(lib_hist)
            lib_sd = stdev(lib_hist)

            print 'p25', p25
            print 'p75', p75
            print 'mean', lib_mean
            print 'sd', lib_sd

            low = int(p25 - mapping_bound * (p75 - p25) + .499)
            high = int(p75 + mapping_bound * (p75 - p25) + .499)

            if (low > high):
                sys.stderr.write("Error: insert size upper bound is smaller than read length\n");
                exit(1)
            
            # ensure z-score limits do not exceed IQR limits
            low = min(low, int(lib_mean - max_stddev * lib_sd))
            high = max(high, int(lib_mean + max_stddev * lib_sd))
        
            # minimum of zero
            if (low < 0):
                low = 0

            # assign the readgroups
            for rg in lib['readgroups']:
                proper[rg] = [low, high]

    # set input file
    if bamfile == None: 
        if is_sam:
            in_bam = pysam.Samfile("-", "r")
        else:
            in_bam = pysam.Samfile('-', 'rb')
    else:
        if is_sam:
            in_bam = pysam.Samfile(bamfile, 'r')
        else:
            in_bam = pysam.Samfile(bamfile, "rb")

    # set output file
    if uncompressed_out:
        out_bam = pysam.Samfile('-', 'wbu', template=in_bam)
    elif bam_out:
        out_bam = pysam.Samfile('-', 'wb', template=in_bam)
    else:
        out_bam = pysam.Samfile('-', 'wh', template=in_bam)
        
    print proper
    for al in in_bam:
        # out_bam.write(al)
        print al

        if al.is_supplementary:
            pass

        elif al.is_unmapped or al.mate_is_unmapped:
            if al.is_proper_pair:
                print 'mismarked proper (unmapped)'
            al.is_proper_pair = False

        elif al.reference_id != al.next_reference_id:
            if al.is_proper_pair:
                print 'mismarked proper (chrom)'
            al.is_proper_pair = False

        elif (al.reference_start < al.next_reference_start
            and (al.is_reverse or not al.mate_is_reverse)):
            if al.is_proper_pair:
                print 'mismarked proper (orient +)'
            al.is_proper_pair = False

        elif (al.reference_start > al.next_reference_start
            and (not al.is_reverse or al.mate_is_reverse)):
            if al.is_proper_pair:
                print 'mismarked proper (orient -)'
            al.is_proper_pair = False

        # if al.supp
        elif (al.template_length >= proper[al.opt('RG')][0]
              and al.template_length <= proper[al.opt('RG')][1]):
            if not al.is_proper_pair:
                print 'mismarked improper (insert size)'
            al.is_proper_pair = True
        else:
            if al.is_proper_pair:
                print 'mismarked proper (insert size)'
            al.is_proper_pair = False

        # out_bam.write(al)
        print al
        # print proper[al.opt('RG')], al.template_length
        # print al
        # # must be in a user specified readgroup
        # if al.opt('RG') not in rg_list:
        #     continue

        # # write out alignment
        # out_bam.write(al)
        # counter += 1
        
        # # bail if reached limit
        # if (limit != None
        #     and counter >= limit):
        #     break


# ============================================
# functions
# ============================================

# class that holds reads from a sequence fragment
class Namegroup():
    def __init__(self, al):
        self.alignments = list()
        self.name = al.qname
        self.sa = 0
        self.num_prim = 0
        self.add_alignment(al)

    def add_alignment(self, al):
        self.alignments.append(al)
        if not al.is_secondary:
            self.num_prim += 1
            try:
                self.sa += len(al.opt('SA').rstrip(';').split(';'))
                # print self.sa
            except KeyError:
                pass

    def is_complete(self):
        return self.num_prim == 2 and len(self.alignments) == self.sa + 2

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
bamfixflags.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: filter readgroup(s) from a BAM file")
    parser.add_argument('-i', '--input', metavar='BAM', required=False, help='Input BAM file')
    parser.add_argument('-l', '--lib_info', metavar='FILE', dest='lib_info_file', type=argparse.FileType('r'), required=True, default=None, help='SVTyper JSON file of library information')
    parser.add_argument('-n', metavar='INT', type=int, default=None, required=False, help='Output first n alignments and quit')
    parser.add_argument('-S', required=False, action='store_true', help='Input is SAM format')
    parser.add_argument('-b', required=False, action='store_true', help='Output BAM format')
    parser.add_argument('-u', required=False, action='store_true', help='Output uncompressed BAM format (implies -b)')

    # parse the arguments
    args = parser.parse_args()

    # bail if no BAM file
    if args.input is None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
    
    # send back the user input
    return args

# ============================================
# driver
# ============================================

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():
    args = get_args()
    bamfixflags(args.input,
                args.lib_info_file,
                args.n,
                args.S,
                args.b,
                args.u)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
    
