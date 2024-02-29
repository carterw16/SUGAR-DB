#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 13:53:39 2021

@author: emfoster
"""

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from classes.Nodes import Nodes


class Map():
    # def map_RMSS(nodes, xfmr, reg, switch, load, oh_lines, ug_lines, tplx_line, variable_to_plot, fuses, V):
    #     # This will take awhile to run FYI 
    #     sys_graph = nx.Graph()
    #     n_of_nodes = len(nodes)   
    #     max_val = np.empty(n_of_nodes)
    #     node_size = []
        
    #     for index in range(0, n_of_nodes):
    #         max_val[index] = # Here, pick the maximum value at that node
    #         sys_graph.add_node(nodes[index].name)
    #         nx.set_node_attributes(sys_graph, {nodes[index].name: max_val[index]}, name = 'infeas') 
        
    #     n_of_oh_lines = len(oh_lines)
    #     for index in range(0, n_of_oh_lines):
    #         sys_graph.add_edge(oh_lines[index].from_node, oh_lines[index].to_node, length = oh_lines[index].length)
            
    #     n_of_ug_lines = len(ug_lines)
    #     for index in range(0, n_of_ug_lines):
    #         sys_graph.add_edge(ug_lines[index].from_node, ug_lines[index].to_node, length = ug_lines[index].length)
        
    #     n_of_tplx_lines = len(tplx_line)
    #     for index in range(0, n_of_tplx_lines):
    #         sys_graph.add_edge(tplx_line[index].from_node, tplx_line[index].to_node, length = tplx_line[index].length)
        
    #     n_of_fuses = len(fuses)
    #     for index in range(0, n_of_fuses):
    #         if fuses[index]._BLOWN == 0:
    #             sys_graph.add_edge(fuses[index].from_node, fuses[index].to_node, length = 0.0001)
        
    #     n_of_xfmr = len(xfmr)
    #     for index in range(0, n_of_xfmr):
    #         sys_graph.add_edge(xfmr[index].from_node, xfmr[index].to_node)
            
    #     n_of_reg = len(reg)
    #     for index in range(0, n_of_reg):
    #         sys_graph.add_edge(reg[index].from_node, reg[index].to_node)
        
    #     n_of_switch = len(switch)
    #     for index in range(0, n_of_switch):
    #         if switch[index]._CLOSED == 1:
    #             sys_graph.add_edge(switch[index].to_node, switch[index].from_node) 

    #     colors = plt.cm.Reds #plt.cm.coolwarm
        
    #     np.random.seed(5)
    #     max_val_node_values = nx.get_node_attributes(sys_graph,'infeas').values()
    #     low = 0
    #     high = max(infeas)
    #     center = 1/2*high 
    #     offset = mcolors.DivergingNorm(vcenter = center, vmin = low, vmax = high)
    #     sm = plt.cm.ScalarMappable(cmap=colors, norm= offset)
    #     pos = nx.nx_agraph.graphviz_layout(sys_graph)
        
    #     # If network is big, change the figure size to something larger:
    #     # plt.figure(3,figsize=(8,12)) 
        
    #     # When drawing with nx.draw:
    #     # - You might want to change the node_size parameter to something smaller
    #     # - You can modify the edge size, but I'm not using that parameter right now
    #     nx.draw(sys_graph, pos, node_size = 50, 
    #             node_color = [sm.to_rgba(i) for i in max_val_node_values], with_labels = False)
    #     ax = plt.gca() # to get the current axis
    #     ax.collections[0].set_edgecolor("black") # Sets the color of the outline of nodes, limes
    #     sm._A = []
    #     cbar = plt.colorbar(sm, fraction = 0.046, pad = 0.04, label = "ADD LABEL")
    #     cbar.set_label("ADD LABEL", fontsize=15)
    #     plt.title("ADD TITLE", fontsize = 15)
    #     plt.show()
    #     # plt.savefig("L2LF5mismatch1.png", format="PNG", dpi = 1000, bbox_inches='tight')   
        

    
    def map_infeasibility(nodes, xfmr, reg, switch, load, oh_lines, ug_lines, tplx_line, infeasibility, fuses, V):
        sys_graph = nx.Graph()
        n_of_nodes = len(nodes)   
        max_curr = np.empty(n_of_nodes)
        node_size = []
        if infeasibility.obj == 1:
            for index in range(0, n_of_nodes):
                if nodes[index].isTriplex == False and nodes[index].bustype != 3:
                    infeas_curr_Ar = V[nodes[index].nodeA_if_r_plus] -  V[nodes[index].nodeA_if_r_minus]
                    infeas_curr_Br = V[nodes[index].nodeB_if_r_plus] -  V[nodes[index].nodeB_if_r_minus]                    
                    infeas_curr_Cr = V[nodes[index].nodeC_if_r_plus] -  V[nodes[index].nodeC_if_r_minus]                     
                    infeas_curr_Ai = V[nodes[index].nodeA_if_i_plus] -  V[nodes[index].nodeA_if_i_minus]
                    infeas_curr_Bi = V[nodes[index].nodeB_if_i_plus] -  V[nodes[index].nodeB_if_i_minus]                    
                    infeas_curr_Ci = V[nodes[index].nodeC_if_i_plus] -  V[nodes[index].nodeC_if_i_minus] 
                    
                    magA = np.sqrt(infeas_curr_Ar**2 + infeas_curr_Ai**2)
                    magB = np.sqrt(infeas_curr_Br**2 + infeas_curr_Bi**2)
                    magC = np.sqrt(infeas_curr_Cr**2 + infeas_curr_Ci**2)
                    infeas_curr = [magA, magB, magC]
                    max_curr[index] = max(np.abs(infeas_curr)) 
                    if abs(max_curr[index]) >= 1e-3:
                        node_size.append(50)
                    else:
                        node_size.append(50)
                    node_name = nodes[index].name.replace('R1-12.47-3_', "")
                    sys_graph.add_node(node_name)
                    nx.set_node_attributes(sys_graph, {node_name: max_curr[index]}, name = 'infeas')
                elif nodes[index].isTriplex == True and nodes[index].bustype != 3:
                    infeas_curr_1r = V[nodes[index].node1_if_r_plus] -  V[nodes[index].node1_if_r_minus]
                    infeas_curr_2r = V[nodes[index].node2_if_r_plus] -  V[nodes[index].node2_if_r_minus]                                        
                    infeas_curr_1i = V[nodes[index].node1_if_i_plus] -  V[nodes[index].node1_if_i_minus]
                    infeas_curr_2i = V[nodes[index].node2_if_i_plus] -  V[nodes[index].node2_if_i_minus]                    
                    
                    mag1 = np.sqrt(infeas_curr_1r**2 + infeas_curr_1i**2)
                    mag2 = np.sqrt(infeas_curr_2r**2 + infeas_curr_1i**2)
                    
                    infeas_curr = [mag1, mag2]
                    max_curr[index] = max(np.abs(infeas_curr))
                    if abs(max_curr[index]) >= 1e-3:
                        node_size.append(50)
                    else:
                        node_size.append(50)
                    node_name = nodes[index].name.replace('R1-12.47-3_', "")
                    sys_graph.add_node(node_name)
                    nx.set_node_attributes(sys_graph, {node_name: max_curr[index]}, name = 'infeas')
                else:
                    node_name = nodes[index].name.replace('R1-12.47-3_', "")
                    sys_graph.add_node(node_name)
                    nx.set_node_attributes(sys_graph, {node_name: 0}, name = 'infeas')
                    node_size.append(50)
        
        if infeasibility.obj == 2:
            for index in range(0, n_of_nodes):
                if nodes[index].isTriplex == False and nodes[index].bustype != 3:
                    magA = np.sqrt(V[nodes[index].nodeA_if_r]**2 + V[nodes[index].nodeA_if_i]**2)
                    magB = np.sqrt(V[nodes[index].nodeB_if_r]**2 + V[nodes[index].nodeB_if_i]**2)
                    magC = np.sqrt(V[nodes[index].nodeC_if_r]**2 + V[nodes[index].nodeC_if_i]**2)
                    infeas_curr = [magA, magB, magC]
                    max_curr[index] = max(np.abs(infeas_curr)) 
                    node_name = nodes[index].name.replace('R1-12.47-3_', "")
                    sys_graph.add_node(node_name)
                    nx.set_node_attributes(sys_graph, {node_name: max_curr[index]}, name = 'infeas')
                elif nodes[index].isTriplex == True and nodes[index].bustype != 3:
                    mag1 = np.sqrt(V[nodes[index].node1_if_r]**2 + V[nodes[index].node1_if_i]**2)
                    mag2 = np.sqrt(V[nodes[index].node2_if_r]**2 + V[nodes[index].node2_if_i]**2)
                    infeas_curr = [mag1, mag2]
                    max_curr[index] = max(np.abs(infeas_curr)) 
                    
                    node_name = nodes[index].name.replace('R1-12.47-3_', "")
                    sys_graph.add_node(nodes[index].name)
                    nx.set_node_attributes(sys_graph, {nodes[index].name: max_curr[index]}, name = 'infeas') 
                else:
                    sys_graph.add_node(nodes[index].name)
                    nx.set_node_attributes(sys_graph, {nodes[index].name: 0}, name = 'infeas')    

        n_of_oh_lines = len(oh_lines)
        for index in range(0, n_of_oh_lines):
            sys_graph.add_edge(oh_lines[index].from_node, oh_lines[index].to_node, length = oh_lines[index].length)
            
        n_of_ug_lines = len(ug_lines)
        for index in range(0, n_of_ug_lines):
            sys_graph.add_edge(ug_lines[index].from_node, ug_lines[index].to_node, length = ug_lines[index].length)
        
        n_of_tplx_lines = len(tplx_line)
        for index in range(0, n_of_tplx_lines):
            sys_graph.add_edge(tplx_line[index].from_node, tplx_line[index].to_node, length = tplx_line[index].length)
        
        n_of_fuses = len(fuses)
        for index in range(0, n_of_fuses):
            if fuses[index]._BLOWN == 0:
                sys_graph.add_edge(fuses[index].from_node, fuses[index].to_node, length = 0.0001)
        
        n_of_xfmr = len(xfmr)
        for index in range(0, n_of_xfmr):
            sys_graph.add_edge(xfmr[index].from_node, xfmr[index].to_node)
            
        n_of_reg = len(reg)
        for index in range(0, n_of_reg):
            sys_graph.add_edge(reg[index].from_node, reg[index].to_node)
        
        n_of_switch = len(switch)
        for index in range(0, n_of_switch):
            if switch[index]._CLOSED == 1:
                sys_graph.add_edge(switch[index].to_node, switch[index].from_node) 

        # Below you can set the array of colors in the heatmap to a number of options
        colors = plt.cm.Reds #plt.cm.coolwarm
        colors.set_bad("black")

        # I set a seed here so that the network layout looks the same everytime. You might play
        # around with the seed to find one that you like.
        np.random.seed(5)

        # Here I'm creating the values that will eventually go to the heatmap. I'm pulling from 
        # the attributes I created earlier.
        infeas = nx.get_node_attributes(sys_graph,'infeas').values()
        low = 0
        high = max(infeas)
        # np.mean(max_infeas_curr)
        center = 1/2*high 
        mfeas = np.array(max_curr)
        # Here is where I'm setting my cutoff parameter
        mfeas[abs(mfeas) < 1e-4] = 0
        mfeas = np.where(mfeas > 0, mfeas, np.nan)
        
        # Here I'm setting the bounds of the heatmap
        offset = mcolors.LogNorm(vmin = np.nanmin(mfeas), vmax = np.nanmax(mfeas))
       
        # The below line is an alternative method for the offset
        # offset = mcolors.DivergingNorm(vcenter = center, vmin = low, vmax = high)

        sm = plt.cm.ScalarMappable(cmap=colors, norm= offset)

        # You can play around with network layouts to find what works best, but here
        # this was the option I liked most
        pos = nx.nx_agraph.graphviz_layout(sys_graph)
        
        # If you're saving your figure, you might want to declare the size
        # plt.figure(3,figsize=(8,12)) 
        
        # The parameters below are the graph, network layout, size of nodes, node colors, and labels
        # You can also do edge color, edge size, font, etc.
        nx.draw(sys_graph, pos, node_size = 50, 
                node_color = [sm.to_rgba(i) for i in infeas], with_labels = False)
        ax = plt.gca() # to get the current axis
        ax.collections[0].set_edgecolor("black")
        sm._A = []
        cbar = plt.colorbar(sm, fraction = 0.046, pad = 0.04, label = r"$|i_{f,k}^{\Omega}|$ in Amps") #, ticks = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2])
        cbar.set_label(r"$|i_{f,k}^{\Omega}|$ in Amps", fontsize=15)
        cbar.ax.tick_params(labelsize = 15)
        #bounds = np.power(10.0, np.arange(-3,2))
        #cbar.ax.set_yticklabels([f'$10^{{{np.log10(b): .0f}}}$' for b in bounds])
        plt.title(r"Greatest Infeasibility Current Magnitude per Node", fontsize = 15)
        plt.show()

        # If you want to save the figure, use the below line and comment out "plt.show()" above
        # plt.savefig("L1ev.png", format="PNG", dpi = 1000, bbox_inches='tight')      
        
        
    def map_normal_infeasibility(nodes, xfmr, reg, switch, load, oh_lines, ug_lines, tplx_line, infeasibility, fuses, V, max_infeas_curr):
        name_to_strip = ''
        sys_graph = nx.Graph()
        n_of_nodes = len(nodes)   

        for index in range(0, n_of_nodes):
            if nodes[index].isTriplex == False and nodes[index].bustype != 3:                
                node_name = nodes[index].name.replace(name_to_strip, "")
                sys_graph.add_node(node_name)
                nx.set_node_attributes(sys_graph, {node_name: max_infeas_curr[index]}, name = 'infeas')
            elif nodes[index].isTriplex == True and nodes[index].bustype != 3:
                node_name = nodes[index].name.replace(name_to_strip, "")
                sys_graph.add_node(node_name)
                nx.set_node_attributes(sys_graph, {node_name: max_infeas_curr[index]}, name = 'infeas')
            else:
                node_name = nodes[index].name.replace(name_to_strip, "")
                sys_graph.add_node(node_name)
                nx.set_node_attributes(sys_graph, {node_name: 0}, name = 'infeas')
        
        n_of_oh_lines = len(oh_lines)
        for index in range(0, n_of_oh_lines):
            from_node_name = oh_lines[index].from_node.replace(name_to_strip, "")
            to_node_name = oh_lines[index].to_node.replace(name_to_strip, "")
            sys_graph.add_edge(from_node_name, to_node_name, length = oh_lines[index].length)
            
        n_of_ug_lines = len(ug_lines)
        for index in range(0, n_of_ug_lines):
            from_node_name = ug_lines[index].from_node.replace(name_to_strip, "")
            to_node_name = ug_lines[index].to_node.replace(name_to_strip, "")
            sys_graph.add_edge(from_node_name, to_node_name, length = ug_lines[index].length)
        
        n_of_tplx_lines = len(tplx_line)
        for index in range(0, n_of_tplx_lines):
            from_node_name = tplx_line[index].from_node.replace(name_to_strip, "")
            to_node_name = tplx_line[index].to_node.replace(name_to_strip, "")
            sys_graph.add_edge(from_node_name, to_node_name, length = tplx_line[index].length)
        
        n_of_fuses = len(fuses)
        for index in range(0, n_of_fuses):
            if fuses[index]._BLOWN == 0:
                from_node_name = fuses[index].from_node.replace(name_to_strip, "")
                to_node_name = fuses[index].to_node.replace(name_to_strip, "")
                sys_graph.add_edge(from_node_name, to_node_name, length = 0.0001)
        
        n_of_xfmr = len(xfmr)
        for index in range(0, n_of_xfmr):
            node_name = xfmr[index].name.replace('name_to_strip', "")
            sys_graph.add_node(node_name)
            nx.set_node_attributes(sys_graph, {node_name: 0}, name = 'infeas')
            # color_map.append('orange')
            from_node_name = xfmr[index].from_node.replace(name_to_strip, "")
            to_node_name = xfmr[index].to_node.replace(name_to_strip, "")
            sys_graph.add_edge(from_node_name, node_name)
            sys_graph.add_edge(to_node_name, node_name)
                 
        n_of_reg = len(reg)
        for index in range(0, n_of_reg):
            node_name = reg[index].name.replace(name_to_strip, "")
            sys_graph.add_node(node_name)
            nx.set_node_attributes(sys_graph, {node_name: 0}, name = 'infeas')
            # color_map.append('orange')
            from_node_name = reg[index].from_node.replace(name_to_strip, "")
            to_node_name = reg[index].to_node.replace(name_to_strip, "")
            sys_graph.add_edge(from_node_name, node_name)
            sys_graph.add_edge(to_node_name, node_name) 
        
        n_of_switch = len(switch)
        for index in range(0, n_of_switch):
            if switch[index]._CLOSED == 1:
                node_name = switch[index].name.replace(name_to_strip, "")
                sys_graph.add_node(node_name)
                nx.set_node_attributes(sys_graph, {node_name: 0}, name = 'infeas')
                # color_map.append('yellow')
                from_node_name = switch[index].from_node.replace(name_to_strip, "")
                to_node_name = switch[index].to_node.replace(name_to_strip, "")
                sys_graph.add_edge(from_node_name, node_name)
                sys_graph.add_edge(to_node_name, node_name)   

        colors = plt.cm.Reds #plt.cm.coolwarm
        colors.set_bad("black")
        np.random.seed(5)
        infeas = nx.get_node_attributes(sys_graph,'infeas').values()
        low = 0
        #np.mean(max_infeas_curr)
        high = max(infeas)
        center = 1/2*high 
        mfeas = np.array(max_infeas_curr)
        mfeas = np.where(mfeas > 0, mfeas, np.nan)
        
        offset = mcolors.LogNorm(vmin = np.nanmin(mfeas), vmax = np.nanmax(mfeas))
        #offset = mcolors.DivergingNorm(vcenter = center, vmin = low, vmax = high)
        sm = plt.cm.ScalarMappable(cmap=colors, norm= offset)
        pos = nx.nx_agraph.graphviz_layout(sys_graph)
        nx.draw(sys_graph, pos, node_size = 50, node_color = [sm.to_rgba(i) for i in infeas], with_labels = False)
        ax = plt.gca() # to get the current axis
        ax.collections[0].set_edgecolor("black")
        sm._A = []
        plt.colorbar(sm)
        plt.show()
        #plt.savefig("L1_LF5_normal.png", format="PNG", dpi = 1000, bbox_inches='tight')
 
    def highlight_specific_node(nodes, xfmr, reg, switch, load, oh_lines, ug_lines, tplx_line, fuses, V, list_of_names):
        sys_graph = nx.Graph()
        color_map = []
        n_of_nodes = len(nodes)   
        
        for index in range(0, n_of_nodes):
            if nodes[index].name == list_of_names:
                sys_graph.add_node(nodes[index].name)
                color_map.append("red")
            else:
                sys_graph.add_node(nodes[index].name)
                color_map.append("grey")
        
        n_of_oh_lines = len(oh_lines)
        for index in range(0, n_of_oh_lines):
            sys_graph.add_edge(oh_lines[index].from_node, oh_lines[index].to_node, length = oh_lines[index].length)
            
        n_of_ug_lines = len(ug_lines)
        for index in range(0, n_of_ug_lines):
            sys_graph.add_edge(ug_lines[index].from_node, ug_lines[index].to_node, length = ug_lines[index].length)
        
        n_of_tplx_lines = len(tplx_line)
        for index in range(0, n_of_tplx_lines):
            sys_graph.add_edge(tplx_line[index].from_node, tplx_line[index].to_node, length = tplx_line[index].length)
        
        n_of_fuses = len(fuses)
        for index in range(0, n_of_fuses):
            if fuses[index]._BLOWN == 0:
                sys_graph.add_edge(fuses[index].from_node, fuses[index].to_node, length = 0.0001)
        
        n_of_xfmr = len(xfmr)
        for index in range(0, n_of_xfmr):
            sys_graph.add_edge(xfmr[index].from_node, xfmr[index].to_node)
            
        n_of_reg = len(reg)
        for index in range(0, n_of_reg):
            sys_graph.add_edge(reg[index].from_node, reg[index].to_node)
        
        n_of_switch = len(switch)
        for index in range(0, n_of_switch):
            if switch[index]._CLOSED == 1:
                sys_graph.add_edge(switch[index].to_node, switch[index].from_node)  
        
        np.random.seed(5)
        pos = nx.nx_agraph.graphviz_layout(sys_graph)
        nx.draw(sys_graph, pos, node_size = 50, node_color = color_map, with_labels = False) 
        plt.show()
             
    def make_map(nodes, oh_lines, ug_lines, tplx_line, xfmr, reg, fuses, switch):
        # Make_map should be run for the overall system (i.e. ieee 13 bus will have 13 bus map)
    
        sys_graph = nx.Graph()
        # For color map, we are going to assign a specific color to go along with each element
        color_map = []
        n_of_nodes = len(nodes)
        
        for index in range(0, n_of_nodes):
            if nodes[index].bustype == 3: # Slack bus
                sys_graph.add_node(nodes[index].name)
                color_map.append('blue')
            elif nodes[index].bustype == 2: # PV bus
                sys_graph.add_node(nodes[index].name)
                color_map.append('red')
            elif nodes[index].bustype == 1: #PQ bus
                sys_graph.add_node(nodes[index].name)
                color_map.append('green')
        
        n_of_oh_lines = len(oh_lines)
        for index in range(0, n_of_oh_lines):
            sys_graph.add_edge(oh_lines[index].from_node, oh_lines[index].to_node, length = oh_lines[index].length)
            
        n_of_ug_lines = len(ug_lines)
        for index in range(0, n_of_ug_lines):
            sys_graph.add_edge(ug_lines[index].from_node, ug_lines[index].to_node, length = ug_lines[index].length)
        
        n_of_tplx_lines = len(tplx_line)
        for index in range(0, n_of_tplx_lines):
            sys_graph.add_edge(tplx_line[index].from_node, tplx_line[index].to_node, length = tplx_line[index].length)
            
        n_of_fuses = len(fuses)
        for index in range(0, n_of_fuses):
            if fuses[index]._BLOWN == 0:
                sys_graph.add_edge(fuses[index].from_node, fuses[index].to_node, length = 0.0001)
        
        n_of_xfmr = len(xfmr)
        for index in range(0, n_of_xfmr):
            sys_graph.add_node(xfmr[index].name)
            color_map.append('orange')
            sys_graph.add_edge(xfmr[index].from_node, xfmr[index].name)
            sys_graph.add_edge(xfmr[index].to_node, xfmr[index].name)

            # An alternative approach here would to not add the transformer as a node, and then just 
            # connect from_node to to_node (this is best if you want to focus on just the original
            # identified nodes i.e. you want just 13 nodes for the ieee13 case) [see map_infeasibility]
                 
        n_of_reg = len(reg)
        for index in range(0, n_of_reg):
            sys_graph.add_node(reg[index].name)
            color_map.append('orange')
            sys_graph.add_edge(reg[index].from_node, reg[index].name)
            sys_graph.add_edge(reg[index].to_node, reg[index].name)

            # An alternative approach here would to not add the transformer as a node, and then just 
            # connect from_node to to_node (this is best if you want to focus on just the original
            # identified nodes i.e. you want just 13 nodes for the ieee13 case) [see map_infeasibility]
        
        n_of_switch = len(switch)
        for index in range(0, n_of_switch):
            sys_graph.add_node(switch[index].name)
            color_map.append('yellow')
            sys_graph.add_edge(switch[index].from_node, switch[index].name)
            sys_graph.add_edge(switch[index].to_node, switch[index].name)
            
            # An alternative approach here would to not add the transformer as a node, and then just 
            # connect from_node to to_node (this is best if you want to focus on just the original
            # identified nodes i.e. you want just 13 nodes for the ieee13 case) [see map_infeasibility]
            
        np.random.seed(5)
        pos = nx.nx_agraph.graphviz_layout(sys_graph)
        nx.draw(sys_graph, pos, node_size = 10, node_color = color_map, with_labels = False)
        plt.show()
        
        return sys_graph

    def make_KVA_map(nodes, oh_lines, ug_lines, tplx_line, xfmr, reg, fuses, switch):
        # Make_map should be run for the overall system (i.e. ieee 13 bus will have 13 bus map)
        base_voltage = []
        n_of_nodes = len(nodes)
        for index in range(0, n_of_nodes):
            node_base_voltage = round(nodes[index].Vbase_LN*1000,2)
            if node_base_voltage not in base_voltage:
                base_voltage.append(node_base_voltage)
        
        base_voltage.sort(reverse = True)
        color_dict = ["red", "orange", "yellow", "yellowgreen", "green", "aqua", "blue", "purple"]
        print('Base Voltages', base_voltage)
        print('Colors', color_dict)
        
        base_voltage = np.array(base_voltage)
        sys_graph = nx.Graph()
        color_map = []
        
        for index in range(0, n_of_nodes):
            sys_graph.add_node(nodes[index].name)
            color_for_node = color_dict[int(np.where(round(nodes[index].Vbase_LN*1000,2) == base_voltage)[0])]
            color_map.append(color_for_node)
            if color_for_node == "blue":
                print("blue", nodes[index].name)
            if color_for_node == "green":
                print("green", nodes[index].name)

        
        n_of_oh_lines = len(oh_lines)
        for index in range(0, n_of_oh_lines):
            sys_graph.add_edge(oh_lines[index].from_node, oh_lines[index].to_node, length = oh_lines[index].length)
            
        n_of_ug_lines = len(ug_lines)
        for index in range(0, n_of_ug_lines):
            sys_graph.add_edge(ug_lines[index].from_node, ug_lines[index].to_node, length = ug_lines[index].length)
        
        n_of_tplx_lines = len(tplx_line)
        for index in range(0, n_of_tplx_lines):
            sys_graph.add_edge(tplx_line[index].from_node, tplx_line[index].to_node, length = tplx_line[index].length)
            
        n_of_fuses = len(fuses)
        for index in range(0, n_of_fuses):
            if fuses[index]._BLOWN == 0:
                sys_graph.add_edge(fuses[index].from_node, fuses[index].to_node, length = 0.0001)
        
        n_of_xfmr = len(xfmr)
        for index in range(0, n_of_xfmr):
            sys_graph.add_node(xfmr[index].name)
            color_map.append('orange')
            sys_graph.add_edge(xfmr[index].from_node, xfmr[index].name)
            sys_graph.add_edge(xfmr[index].to_node, xfmr[index].name)
                 
        n_of_reg = len(reg)
        for index in range(0, n_of_reg):
            sys_graph.add_node(reg[index].name)
            color_map.append('orange')
            sys_graph.add_edge(reg[index].from_node, reg[index].name)
            sys_graph.add_edge(reg[index].to_node, reg[index].name)
        
        n_of_switch = len(switch)
        for index in range(0, n_of_switch):
            sys_graph.add_node(switch[index].name)
            color_map.append('yellow')
            sys_graph.add_edge(switch[index].from_node, switch[index].name)
            sys_graph.add_edge(switch[index].to_node, switch[index].name)
            
            
        np.random.seed(5)
        pos = nx.nx_agraph.graphviz_layout(sys_graph)
        nx.draw(sys_graph, pos, node_size = 10, node_color = color_map, with_labels = False) 
        plt.show()
        
        return sys_graph        
    def feas_map(Ylin_row, Ylin_col, Ynlin_row, Ynlin_col, Ylin_h_row, Ylin_h_col, nodes, V_error):
            sys_graph = nx.Graph()
            color_map = []
            for index in range(0, len(Ylin_row)):
                sys_graph.add_edge(Ylin_row[index], Ylin_col[index])
            for index in range(0, len(Ynlin_row)):
                sys_graph.add_edge(Ynlin_row[index], Ynlin_col[index])
            for index in range(0, len(Ylin_h_row)):
                sys_graph.add_edge(Ylin_h_row[index], Ylin_h_col[index])
            

            for index in range(0, len(V_error)):
                if V_error[index] <= 1e6:
                    sys_graph.add_node(index)
                    color_map.append('green')
                elif V_error[index] > 1e6 and V_error[index] <= 1e3:
                    sys_graph.add_node(index)
                    color_map.append('yellow')
                elif V_error[index] > 1e3:
                    sys_graph.add_node(index)
                    color_map.append('red')
        
            nx.draw(sys_graph, node_color=color_map, width = 3, font_weight='bold', with_labels=True)
            plt.show()
            
            sys_graph.clear()
    
    def make_dual_map(nodes, oh_lines, ug_lines, xfmr, reg, switch):
        # Only Phase A
        sys_graph = nx.Graph()
        color_map = []
        
        n_of_nodes = len(nodes)
        for index in range(0, n_of_nodes):
            if nodes[index].bustype == 3: # Slack bus
                sys_graph.add_node(nodes[index].nodeA_Vr)
                color_map.append('#4682B4')
                sys_graph.add_node(nodes[index].nodeA_Vi)
                color_map.append('#F08080')
                sys_graph.add_node(nodes[index].nodeA_Lr)
                color_map.append('#CC99CC')
                sys_graph.add_node(nodes[index].nodeA_Li)
                color_map.append('#FFA54F')
            elif nodes[index].bustype == 1: #PQ bus
                sys_graph.add_node(nodes[index].nodeA_Vr)
                color_map.append('#0D4F8B')
                sys_graph.add_node(nodes[index].nodeA_Vi)
                color_map.append('#EE0000')
                sys_graph.add_node(nodes[index].nodeA_Lr)
                color_map.append('#6B238E')
                sys_graph.add_node(nodes[index].nodeA_Li)
                color_map.append('#FF6600')
                
                if 'l' in nodes[index].ID:
                    nodeA_Vr = nodes[index].nodeA_Vr 
                    nodeA_Vi = nodes[index].nodeA_Vi 
                    nodeA_Lr = nodes[index].nodeA_Lr 
                    nodeA_Li = nodes[index].nodeA_Li                                        

                    sys_graph.add_edge(nodeA_Vi, nodeA_Vr)  
                    sys_graph.add_edge(nodeA_Lr, nodeA_Vr)
                    sys_graph.add_edge(nodeA_Lr, nodeA_Vi)
                    sys_graph.add_edge(nodeA_Li, nodeA_Vr)
                    sys_graph.add_edge(nodeA_Li, nodeA_Vi)
                    sys_graph.add_edge(nodeA_Lr, nodeA_Li)
                    
                    
                
        
        n_of_oh_lines = len(oh_lines)
        for index in range(0, n_of_oh_lines):
            nodeA_Vr_from = nodes[Nodes.nodeKey[oh_lines[index].from_node]].nodeA_Vr 
            nodeA_Vi_from = nodes[Nodes.nodeKey[oh_lines[index].from_node]].nodeA_Vi 
            nodeA_Lr_from = nodes[Nodes.nodeKey[oh_lines[index].from_node]].nodeA_Lr 
            nodeA_Li_from = nodes[Nodes.nodeKey[oh_lines[index].from_node]].nodeA_Li 
            
            nodeA_Vr_to = nodes[Nodes.nodeKey[oh_lines[index].to_node]].nodeA_Vr 
            nodeA_Vi_to = nodes[Nodes.nodeKey[oh_lines[index].to_node]].nodeA_Vi 
            nodeA_Lr_to = nodes[Nodes.nodeKey[oh_lines[index].to_node]].nodeA_Lr 
            nodeA_Li_to = nodes[Nodes.nodeKey[oh_lines[index].to_node]].nodeA_Li 
            
            sys_graph.add_edge(nodeA_Vr_from, nodeA_Vr_to)
            sys_graph.add_edge(nodeA_Vr_from, nodeA_Vi_to)
            sys_graph.add_edge(nodeA_Vi_from, nodeA_Vr_to)
            sys_graph.add_edge(nodeA_Vi_from, nodeA_Vi_to)
            sys_graph.add_edge(nodeA_Vi_to, nodeA_Vr_to)
            sys_graph.add_edge(nodeA_Vi_from, nodeA_Vr_from)
            
            sys_graph.add_edge(nodeA_Lr_from, nodeA_Lr_to)
            sys_graph.add_edge(nodeA_Lr_from, nodeA_Li_to)
            sys_graph.add_edge(nodeA_Li_from, nodeA_Lr_to)
            sys_graph.add_edge(nodeA_Li_from, nodeA_Li_to)
            sys_graph.add_edge(nodeA_Li_to, nodeA_Lr_to)
            sys_graph.add_edge(nodeA_Li_from, nodeA_Lr_from)
            
        n_of_ug_lines = len(ug_lines)
        for index in range(0, n_of_ug_lines):
            sys_graph.add_edge(ug_lines[index].from_node, ug_lines[index].to_node)
        
        n_of_xfmr = len(xfmr)
        for index in range(0, n_of_xfmr):
            nodeA_Vr_from = nodes[Nodes.nodeKey[xfmr[index].from_node]].nodeA_Vr
            nodeA_Vi_from = nodes[Nodes.nodeKey[xfmr[index].from_node]].nodeA_Vi 
            nodeA_Vr_to = nodes[Nodes.nodeKey[xfmr[index].to_node]].nodeA_Vr 
            nodeA_Vi_to = nodes[Nodes.nodeKey[xfmr[index].to_node]].nodeA_Vi 
            
            nodeA_Vr_primary = xfmr[index].nodeA_Vr_primary
            nodeA_Vi_primary = xfmr[index].nodeA_Vi_primary
            nodeA_Vr_pos_secondary = xfmr[index].nodeA_Vr_pos_secondary
            nodeA_Vi_pos_secondary = xfmr[index].nodeA_Vi_pos_secondary

            nodeA_Lr_from = nodes[Nodes.nodeKey[xfmr[index].from_node]].nodeA_Lr
            nodeA_Li_from = nodes[Nodes.nodeKey[xfmr[index].from_node]].nodeA_Li 
            nodeA_Lr_to = nodes[Nodes.nodeKey[xfmr[index].to_node]].nodeA_Lr 
            nodeA_Li_to = nodes[Nodes.nodeKey[xfmr[index].to_node]].nodeA_Li 
            
            nodeA_Lr_primary = xfmr[index].nodeA_Lr_primary
            nodeA_Li_primary = xfmr[index].nodeA_Li_primary
            nodeA_Lr_pos_secondary = xfmr[index].nodeA_Lr_pos_secondary
            nodeA_Li_pos_secondary = xfmr[index].nodeA_Li_pos_secondary          

            sys_graph.add_node(nodeA_Vr_primary)
            color_map.append('#B7C3D0')
            
            sys_graph.add_node(nodeA_Vi_primary)
            color_map.append('#6C7B8B')

            sys_graph.add_node(nodeA_Vr_pos_secondary)
            color_map.append('#B7C3D0')

            sys_graph.add_node(nodeA_Vi_pos_secondary)
            color_map.append('#6C7B8B')
            
            sys_graph.add_node(nodeA_Lr_primary)
            color_map.append('#B7C3D0')

            sys_graph.add_node(nodeA_Li_primary)
            color_map.append('#6C7B8B')

            sys_graph.add_node(nodeA_Lr_pos_secondary)
            color_map.append('#B7C3D0')

            sys_graph.add_node(nodeA_Li_pos_secondary)
            color_map.append('#6C7B8B')
            
            sys_graph.add_edge(nodeA_Vr_from, nodeA_Vr_primary)
            sys_graph.add_edge(nodeA_Vi_from, nodeA_Vi_primary)
            sys_graph.add_edge(nodeA_Vr_primary, nodeA_Vr_pos_secondary)
            sys_graph.add_edge(nodeA_Vi_primary, nodeA_Vi_pos_secondary)  
            sys_graph.add_edge(nodeA_Vr_primary, nodeA_Vi_pos_secondary)
            sys_graph.add_edge(nodeA_Vi_primary, nodeA_Vr_pos_secondary) 
            sys_graph.add_edge(nodeA_Vr_primary, nodeA_Vi_primary)
            sys_graph.add_edge(nodeA_Vr_pos_secondary, nodeA_Vi_pos_secondary) 
            sys_graph.add_edge(nodeA_Vr_pos_secondary, nodeA_Vr_to)
            sys_graph.add_edge(nodeA_Vi_pos_secondary, nodeA_Vi_to)  

            sys_graph.add_edge(nodeA_Lr_from, nodeA_Lr_primary)
            sys_graph.add_edge(nodeA_Li_from, nodeA_Li_primary)
            sys_graph.add_edge(nodeA_Lr_primary, nodeA_Lr_pos_secondary)
            sys_graph.add_edge(nodeA_Li_primary, nodeA_Li_pos_secondary)  
            sys_graph.add_edge(nodeA_Lr_primary, nodeA_Li_pos_secondary)
            sys_graph.add_edge(nodeA_Li_primary, nodeA_Lr_pos_secondary) 
            sys_graph.add_edge(nodeA_Lr_primary, nodeA_Li_primary)
            sys_graph.add_edge(nodeA_Lr_pos_secondary, nodeA_Li_pos_secondary) 
            sys_graph.add_edge(nodeA_Lr_pos_secondary, nodeA_Lr_to)
            sys_graph.add_edge(nodeA_Li_pos_secondary, nodeA_Li_to)            
            
            
            # sys_graph.add_edge(nodeA_Vr_primary, nodeA_Vi_pos_secondary)
            # sys_graph.add_edge(nodeA_Vr_primary, nodeA_Vr_pos_secondary)
                 
        # n_of_reg = len(reg)
        # for index in range(0, n_of_reg):
        #     sys_graph.add_node(reg[index].name)
        #     color_map.append('orange')
        #     sys_graph.add_edge(reg[index].from_node, reg[index].name)
        #     sys_graph.add_edge(reg[index].to_node, reg[index].name)
        
        # n_of_switch = len(switch)
        # for index in range(0, n_of_switch):
        #     sys_graph.add_node(switch[index].name)
        #     color_map.append('yellow')
        #     sys_graph.add_edge(switch[index].from_node, switch[index].name)
        #     sys_graph.add_edge(switch[index].to_node, switch[index].name)
            
        
        nx.drawing.layout.bipartite_layout(sys_graph, Nodes.voltage_index)
        nx.draw_kamada_kawai(sys_graph, node_color=color_map, width = 3, font_weight='bold', with_labels=True)
        # nx.draw_kamada_kawai(sys_graph, node_color=color_map, width = 3, font_weight='bold', with_labels=True)
        plt.show()
            
