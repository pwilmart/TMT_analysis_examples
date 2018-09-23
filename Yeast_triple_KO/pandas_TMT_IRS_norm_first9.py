"""program "pandas_TMT_IRS_norm.py"
Reads in TMT summary files from PAW processing and
computes normalizations (within each TMT and between TMT runs).

IRS normalization should be computed now. Aug. 2016, -PW
Converted for Python 3. Aug. 2017, -PW
Modified for PAW protein summary files. -PW 10/4/2017

written by Phil Wilmarth, OHSU, June 2016.
updated for use with PAW pipeline, -PW Oct. 2017.

The MIT License (MIT)

Copyright (c) 2017 Phillip A. Wilmarth and OHSU

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Direct questions to:
Technology & Research Collaborations, Oregon Health & Science University,
Ph: 503-494-8200, FAX: 503-494-4729, Email: techmgmt@ohsu.edu.
"""
# this version modified to do the yeast triple knock out data
# uses the first 9 of 11 channels and creates insted of locates a pooled channel
# correctly leaves the 131 unused channels alone (except for some slight leveling between experiments)

import os
import sys
import re
import functools
import csv

from pandas import Series, DataFrame
import pandas as pd
import numpy as np

import tkinter
from tkinter import filedialog

# any globals defined here
VERSION = 'v1.0.1'

def get_file(default_location, ext_list, title_string=None):
    """Dialog box to browse to a file.  Returns full file name.

    Usage: full_file_name = get_file(default_location, ext_list, [title]),
        where "default_location" is a starting folder location,
        ext_list is a list of (label, pattern) tuples,
        e.g. ext_list = [('Text files', '*.txt')],
        "title" is an optional message to list in the dialog box, and
        "full_file_name" is the complete name of the selected file.
    Written by Phil Wilmarth, OHSU, 2008, 2016.
    """
    # set up GUI elements
    root = tkinter.Tk()
    root.withdraw()
    try:
        root.tk.call('console', 'hide')
    except:
        pass
    
    # set default title string and ext list if not passed
    if title_string is None:   
        title_string = 'Select a single FILE'
    if not ext_list:
        ext_list =  [('All files', '*.*')]
    if not default_location:
        default_location = os.getcwd()
    
    # create dialog box for file selection
    root.update()   # helps make sure dialog box goes away after selection
    filename = filedialog.askopenfilename(parent=root, initialdir=default_location, 
                                          filetypes=ext_list, title=title_string)    
    # return full filename
    return filename 
class TMT_experiment(object):
    """Each TMT experiment has several things to keep track of."""
    def __init__(self, label):
        self.label = label          # TMT experiment name
        self.plex = 0               # number of TMT channels
        self.channels = []          # channel designations
        self.sample_key = {}        # channel designation to biological sample dictionary
        self.grand_totals = []      # column totals for each channel
        self.average_total = 0.0    # average of the column totals
        self.target_total = 0.0     # study-wide global normalization target value
        self.loading_factors = []   # sample loading correction factors
        return

class PAW_TMT_results(object):
    """Data container for PAW TMT protein summary files."""
    def __init__(self):
        self.prefix_lines = []      # any lines before table
        self.suffix_lines = []      # any lines after table
        self.pre_headers = []       # the headers in the row before the headers
        self.headers = []           # the headers from the main header row
        self.frame = None           # pandas DataFrame with main table
        self.num_cols = 0           # number of columns in main table
        self.num_rows = 0           # number of rows in main table
        self.num_prot = 0           # actual number of proteins to normalize
        self.results_file = ''      # original data file's name
        self.TMT_exps = []          # list of TMT_experiment objects for each TMT experiment
        return

    def save_prefix(self, lines):
        """Saves contents of lines before table"""
        lines = [x.rstrip() for x in lines]
        for line in lines:
            if line:
                self.prefix_lines.append(line)
            else:
                self.prefix_lines.append(' ')
        return

    def save_suffix(self, lines):
        """Saves contents of lines after table"""
        lines = [x.rstrip() for x in lines]
        for line in lines:
            if line:
                self.suffix_lines.append(line)
            else:
                self.suffix_lines.append(' ')
        return
        
    def get_sample_keys(self):
        """Makes TMT channel to sample key dictionary."""
        # get the indexes for each experiment
        for exp in self.TMT_exps:            
            exp.channels = [head for head in self.headers if ('TotInt_' in head) and (exp.label in head)]
            exp.plex = len(exp.channels)
            exp.sample_key = {head: self.pre_headers[j] for (j, head) in enumerate(self.headers) if head in exp.channels}            
        return
      
    def load_table(self, file_name):
        """Loads a table from a file into a pandas data frame."""
        self.results_file = file_name
        # read file into memory, remove EOL from lines
        self.contents = open(file_name, 'rt').readlines()
        for i, line in enumerate(self.contents):
            self.contents[i] = '\t'.join([x.strip() if x.strip() else ' ' for x in line.split('\t')])
        
        # find table start and table end indexes
        table_start, table_end = 0, 0
        # start is after header line
        for i, row in enumerate(self.contents):
            if row.startswith('ProtGroup\tCounter'): # header line
                table_start = i + 1
                
                # replace spaces with underscores in header lines
                self.pre_headers = [re.sub(r' ', r'_', x) for x in self.contents[i-1].split('\t')]
                self.headers = [re.sub(r' ', r'_', x) for x in row.split('\t')]
                
                self.num_cols = len(self.headers)
                labels = [x[len('Total_'):] for x in self.headers if x.startswith('Total_')]
                self.TMT_exps = [TMT_experiment(label) for label in labels] # create the TMT exp containers
                break
            
        # end is after the actual data
        for i, row in enumerate(self.contents):
            if i > table_start:
                try:
                    # lines after table are "short" in original TXT files. When Excel is used as an
                    # editor to add sample keys and flag contaminants, the text export has padded lines.
                    # can test for a ProtGroup number in the first cell.
                    x = float(row.split('\t')[0])
                except ValueError:
                    table_end = i
                    break
                
        # trap case where there are no lines after the data
        if table_end == 0:
            table_end = i
        
        # save table length, get sample keys, save prefix and suffix lines
        self.num_rows = table_end - table_start
        self.get_sample_keys()
        self.save_prefix(self.contents[:table_start-2])
        self.save_suffix(self.contents[table_end:])
                   
        # parse columns into dictionaries (header: list) for data frame loading
        table_dict = {header: [] for header in self.headers}
        for line in self.contents[table_start:table_end]:
            for i, cell in enumerate(line.split('\t')): 
                table_dict[self.headers[i]].append(cell)
                
        # cast intensity columns to correct np array data types
        for exp in self.TMT_exps:
            for header in exp.channels:
                table_dict[header] = np.array(table_dict[header], dtype='float64')
                
        # make a data frame from the table column dictionary and return it
        self.frame = DataFrame(table_dict, columns=self.headers)
        return

    def find_study_intersection(self):
        """Finds the intersection of protein identifications."""
        # compute average columns for each TMT experiment and for all
        print('proteins in results file:', len(self.frame))
        total_proteins = len(self.frame['Filter'][self.frame['Filter'] == ' '])
        print('total number of non-contaminant proteins:', total_proteins)
        contams = len(self.frame) - total_proteins
        print('contamaminant and decoy count:', contams)
        self.exp_frames = []
        self.frame['Average_all'] = 0.0
        for exp in self.TMT_exps:
            exp_frame = self.frame[exp.channels].copy()
            exp_frame['Average_' + exp.label] = exp_frame.mean(axis=1)
            self.frame['Average_' + exp.label] = exp_frame.mean(axis=1)
            self.frame['Average_all'] += exp_frame.mean(axis=1) # accumulate running total
            self.exp_frames.append(exp_frame)
        self.frame['Average_all'] /= len(self.exp_frames) # make running total into average

        # get experiment labels for tagging missing data
        self.exp_labels = [x.label for x in self.TMT_exps]
        self.all = 0

        def _make_missing_flag(row):
            """Helper function to make flag strings."""
            flags = []
            for i, intensity in enumerate(row):
                if intensity == 0.0:
                    flags.append(self.exp_labels[i])
            if flags:
                if len(flags) == len(self.exp_labels):
                    flag = 'missing in all'
                    self.all += 1
                else:
                    flag = 'missing in ' + '+'.join(flags)
            else:
                flag = ' '
            return flag
       
        # get the average intensities for each experiment in separate frame       
        ave_cols = [x for x in self.frame.columns if x.startswith('Average_') and not x.endswith('_all')]
        missing_df = self.frame[ave_cols]
        missing_flags = missing_df.apply(_make_missing_flag, axis=1)
        self.frame.insert(7, 'Missing', missing_flags)

        # individual experiment averages of zero are missing in that experiment
        for exp in self.TMT_exps:
            not_quant = len(self.frame.loc[self.frame['Average_'+exp.label] == 0.0])
            print('..%s proteins were not quantifiable in %s' % (not_quant, exp.label))
        print('..proteins with no intensities experiment-wide: %s' % self.all)
        return

    def move_rows_to_bottom(self):
        """Sorts contams and missing reporter ion proteins to bottom of table.
        Also order row by decreasing average intensity."""
        self.frame.loc[self.frame['Filter'] != ' ', 'Average_all'] = 0.0
#        self.frame.loc[self.frame['Missing'] != ' ', 'Average_all'] = 0.0
        self.frame.sort_values(by=['Filter', 'Missing', 'Average_all'], ascending=[True, True, False], inplace=True)

        # set the index
        self.frame.index = list(range(len(self.frame)))

        # get number of rows without the contaminants
        self.num_prot = len(self.frame[self.frame['Filter'] == ' ']) - 1     # dataframe indexing is inclusive
        return

    def compute_loading_norm_factors(self):
        """Compute reporter ion columns totals and loading correction factors."""
        for exp in self.TMT_exps:
            exp.grand_totals = self.frame.loc[:self.num_prot, exp.channels].sum(axis=0)
            exp.average_total = exp.grand_totals[:-2].mean()
        for exp in self.TMT_exps:
            exp.target_total = np.mean([x.average_total for x in self.TMT_exps])
            exp.loading_factors = exp.target_total / exp.grand_totals
            exp.loading_factors[-2:] = exp.target_total / exp.average_total
        print()
        for exp in self.TMT_exps:
            print('%s SL Norm Factors:' % exp.label)
            for ch, fac in zip(exp.channels, exp.loading_factors):
                print('  %s: %0.6f' % (ch, fac))
        return

    def add_loading_normalized_columns(self):
        """Add new columns with sample laoding corrected intensities."""
        all_int_cols = [x.channels for x in self.TMT_exps]
        self.all_int_cols = [item for sublist in all_int_cols for item in sublist]   # get list of all reporter ion columns
        self.work_frame = self.frame.loc[:self.num_prot, self.all_int_cols].copy()    # get sub-frame of just the data
        

        # loop over TMT experiments and correct each channel by single factor
        for j, exp in enumerate(self.TMT_exps):
            exp.SLNorm_cols = []
            for i, channel in enumerate(exp.channels):
                SLNorm_col = 'SLNorm_' + exp.sample_key[channel] + '_' + str(j+1)
                exp.SLNorm_cols.append(SLNorm_col)
                self.work_frame[SLNorm_col] = self.work_frame[channel] * exp.loading_factors[i]

            # compute the average of the pool(s) in each TMT experiment
#            pool_list = [x for x in exp.SLNorm_cols if 'POOL' in x.upper()]
            self.work_frame['AvePool_%d' % (j+1,)] = self.work_frame[exp.SLNorm_cols[:-2]].mean(axis=1)
        return

    def compute_IRS_norm_factors(self):
        """Now do the IRS stuff."""                
        # find the two pooled standard channels and make temporary frame for IRS computations
        ave_pool_list = [x for x in self.work_frame.columns if 'AVEPOOL' in x.upper()]
        temp_frame = self.work_frame.loc[:, ave_pool_list].copy()

        def _geo_mean(row):
            """Computes gemoetric mean, skipping any zeros."""
            net_row = row[row != 0.0]
            try:
                power = 1 / len(net_row)
            except ZeroDivisionError:
                return 0.0
            return net_row.product() ** power        

        # compute the geometric mean of the pools and the correction factor columns        
        temp_frame['GeoMeanPools'] = temp_frame.apply(_geo_mean, axis=1)
        self.IRSFactors = ['IRS_Fac_' + x for x in ave_pool_list]        
        for i, col in enumerate(self.IRSFactors):
            ave_pool = temp_frame[ave_pool_list[i]]
            temp_frame[col] = np.where(ave_pool > 0.0, temp_frame['GeoMeanPools'] / temp_frame[ave_pool_list[i]], 0.0)

#        temp_frame.to_csv(os.path.join(self.here, 'IRS_dump.txt'), quoting=csv.QUOTE_NONE, index=False)

        # add the columns back to the original data frame
        col_list = ['GeoMeanPools'] + self.IRSFactors
        for col in col_list:
            self.work_frame[col] = temp_frame[col].copy()
        return

    def add_IRS_normalized_columns(self):
        """Add the IRS normalized columns."""
        # loop over TMT experiments and correct each protein by IRS factors
        for j, exp in enumerate(self.TMT_exps):
            exp.IRSNorm_cols = []
            for i, channel in enumerate(exp.channels):
                IRSNorm_col = 'IRSNorm_' + exp.sample_key[channel] + '_' + str(j+1)
                exp.IRSNorm_cols.append(IRSNorm_col)
                self.work_frame[IRSNorm_col] = self.work_frame[exp.SLNorm_cols[i]] * self.work_frame[self.IRSFactors[j]]                
        return

######## starts here ##############
print('=================================================================')
print(' pandas_TMT_IRS_norm.py, %s, written by Phil Wilmarth 2017-8 ' % VERSION)    
print('=================================================================')

print('select a protein results file with TMT intensities and sample key:')
default_location = os.getcwd()
results_file = get_file(default_location, [('Text files', '*.txt')],
                        'Select a protein results file with sample keys')
if not results_file:
    sys.exit()

# load the results file
print('processing:', os.path.split(results_file)[1])
frame = PAW_TMT_results()               # create a container for results
frame.load_table(results_file)          # load the results file into container

frame.here = os.path.split(results_file)[0]

frame.find_study_intersection()         # flag proteins not seen in all TMT experiments
frame.move_rows_to_bottom()             # move contaminants and excluded proteins to botton of table
frame.compute_loading_norm_factors()    # compute the sample loading correction factors
frame.add_loading_normalized_columns()  # add the sample-loading-normalized columns
frame.compute_IRS_norm_factors()        # add the IRS normalized column factors
frame.add_IRS_normalized_columns()      # finally add the IRS normalized columns

# dump table and see what we have
"""Output matches by hand analysis.
Need to get column header key added to suffix
Need to get proper prefix and suffix framing the table
"""
# need to drop the raw instensity columns from the work_frame before merging
frame.work_frame.drop(frame.all_int_cols, axis = 1, inplace = True)
final = frame.frame.merge(frame.work_frame, how='left', left_index=True, right_index=True)
path, ext = os.path.splitext(results_file)
new_path = path + '_IRS_normalized' + ext
final.to_csv(new_path, sep='\t', quoting=csv.QUOTE_NONE, index=False)

print('\nnumber of non-contaminant/decoy proteins:', frame.num_prot + 1)
print('number of proteins with full intensity sets:', len(final[(final['Missing'] == ' ') & (final['Filter'] == ' ')]))

# end
    

    
