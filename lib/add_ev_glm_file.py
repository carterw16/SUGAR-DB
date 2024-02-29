#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 01:03:05 2021

@author: emfoster
"""
import pandas as pd
import numpy as np

def add_ev_glm_file(oldfile, newfile, ev_file, lf = None):
    ev_testcases = []
    old_lines = []
    new_lines = []
    
    with open(ev_file,'r') as f:
        for line in f:
            if 'name' in line:
                line_parts = line.split(" ")
                old_name_as_string = line_parts[6]
                old_name_without_semi = old_name_as_string[:-1]
                new_name = old_name_without_semi + '_EV'
                new_line = ' '+ ' '+ ' '+ ' '+ ' '+ line_parts[5] + " " + new_name + ';\n'
                new_lines.append(new_line)
            elif 'parent' in line:
                line_parts = line.split(" ")
                old_name_as_string = line_parts[6]
                old_name_without_semi = old_name_as_string[:-1]
                new_name = old_name_without_semi + '_EV'
                new_line = ' '+ ' '+ ' '+ ' '+ ' '+ line_parts[5] + " " + new_name + ';\n'
                new_lines.append(new_line)
            elif 'power_12' in line:
                line_parts = line.split(" ")
                new_power = 7200+0.0j
                new_power_as_string = str(new_power)              
                new_power_as_string = new_power_as_string[1:-1]
                new_line = ' '+ ' '+ ' '+ ' '+ ' '+ line_parts[5] + " " + str(new_power_as_string) + ';\n' 
                new_lines.append(new_line)    
            elif 'power_2' in line or 'Power_1' in line:
                line_parts = line.split(" ")
                new_power = 1600+0.0j
                new_power_as_string = str(new_power)              
                new_power_as_string = new_power_as_string[1:-1]
                new_line = ' '+ ' '+ ' '+ ' '+ ' '+ line_parts[5] + " " + str(new_power_as_string) + ';\n' 
                new_lines.append(new_line) 
            elif 'to' in line:
                line_parts = line.split(" ")
                old_name_as_string = line_parts[6]
                old_name_without_semi = old_name_as_string[:-1]
                new_name = old_name_without_semi + '_EV'
                new_line = ' '+ ' '+ ' '+ ' '+ ' '+ line_parts[5] + " " + new_name + ';\n'
                new_lines.append(new_line)     
            else:
                new_lines.append(line)
                
                
    NEW_FILENAME = newfile #'/'.join(oldfile.split('/')[:-1]) + '/test.glm'
    with open(NEW_FILENAME, 'a') as f:
        for line in new_lines:
            f.write(line)
            
    return

if __name__ == "__main__":
    oldpath = '/Users/emfoster/Documents/Research/SUGAR3/testcases/gridlabd/R4-12.47-1/R4-12.47-1.glm'
    newpath = '/Users/emfoster/Documents/Research/SUGAR3/testcases/gridlabd/R4-12.47-1_EV_nodes/R4-12.47-1_EV_nodes.glm'
    ev_nodes = '/Users/emfoster/Documents/Research/SUGAR3/testcases/R4-12.47-1_EV_50.glm'
    lf = 1
    add_ev_glm_file(oldpath, newpath, ev_nodes, lf)
    

            