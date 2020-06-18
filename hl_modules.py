#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modules for reading .rpf files for headlands project.

@author: giuliabronzi
"""


# =============================================================================
# Takes in a dataframe, finds average temp between any duplicates,
# removes the duplicates, and returns the dataframe. 
# =============================================================================
def removeDuplicates(df):
    
    #sets temperature to the average (if not a duplicate, it doesn't affect temp)
    df['temperature'] = ((df.temperature.resample('H').max() + 
                              df.temperature.resample('H').min())/2)

    #makes datetime a column so its easy to delete duplicates
    df = df.reset_index()

    #drops datetime duplicates
    df = df.drop_duplicates(subset='datetime')
    df = df.set_index('datetime')

    return df











