#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 17:27:10 2021

@author: emfoster
"""

def change_glm_file(oldfile, newfile, lf):
    new_lines = []
    with open(oldfile,'r') as f:
        for line in f:
            if 'constant_power' in line:
                line_parts = line.split(" ")
                old_power_as_string = line_parts[6]
                old_power_stripped = old_power_as_string.strip()
                old_power_without_semi = old_power_stripped[:-1]
                old_power = complex(old_power_without_semi)
                new_power = lf*old_power
                new_power_as_string = str(new_power)
                if new_power == 0:
                    new_power_as_string = '0+0j'
                else:
                    new_power_as_string = new_power_as_string[1:-1]
                
                new_line = ' '+ ' '+ ' '+ ' '+ ' '+ line_parts[5] + " " + str(new_power_as_string) + ';\n' 
                new_lines.append(new_line)
            elif 'power_1' in line or 'power_2' in line:
                line_parts = line.split(" ")
                old_power_as_string = line_parts[6]
                old_power_stripped = old_power_as_string.strip()
                old_power_without_semi = old_power_stripped[:-1]
                old_power = complex(old_power_without_semi)
                new_power = lf*old_power
                new_power_as_string = str(new_power)
                if new_power == 0:
                    new_power_as_string = '0+0j'
                else:
                    new_power_as_string = new_power_as_string[1:-1]
                
                new_line = ' '+ ' '+ ' '+ ' '+ ' '+ line_parts[5] + " " + str(new_power_as_string) + ';\n' 
                new_lines.append(new_line)                  
            else:
                new_lines.append(line)
    
    NEW_FILENAME = newfile #'/'.join(oldfile.split('/')[:-1]) + '/test.glm'
    with open(NEW_FILENAME, 'w') as f:
        for line in new_lines:
            f.write(line)
            
    return

if __name__ == "__main__":
    oldpath = '/Users/emfoster/Documents/Research/SUGAR3/testcases/gridlabd/R4-12.47-1/R4-12.47-1.glm'
    newpath = '/Users/emfoster/Documents/Research/SUGAR3/testcases/gridlabd/R4-12.47-1_LF4/R4-12.47-1_LF4.glm'
    change_glm_file(oldpath, newpath, 3.2)
    
    