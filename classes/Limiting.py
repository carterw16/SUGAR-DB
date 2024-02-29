#!/usr/bin/env python3z
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 08:21:53 2021

@author: emfoster
"""

import numpy as np
from classes.Nodes import Nodes

import ipdb


class Limiting:
    def __init__(self, limiting_type, features):
        '''
        xmin and xmax are the minimum and maximum allowed values for voltage V,
        either a scalar that is applicable for all voltages or a vector 
        
        dual and x index are the dual and primal values
        
        Limiting is initialized in respective pairs and then later called with apply diode limiting
        '''
        if limiting_type == 'diode':
            try: 
                self.xmin = features['xmin']
                self.mu_min_index = features['mu_min_index']
            except:
                self.xmin = None
                self.mu_min_index = None
            
            try:
                self.xmax = features['xmax']
                self.mu_max_index = features['mu_max_index']
            except:
                self.xmax = None
                self.mu_max_index = None

            self.x_index = features['x_index']
            self.cs_tol = features['cs_tol']
            
        if limiting_type == 'voltage':
            self.maxStep__ = features['maxstep']
            self.minStep__ = features['minstep']
            self.VrLimit__ = features['Vr Limit']
            self.ViLimit__ = features['Vi Limit']
            self.stamp_dual = features['dual']
            self.Vr_norm = features['Vr Norm']
            self.Vi_norm = features['Vi Norm']
            
        if limiting_type == 'variable':
            self.Vmax_lim = features['Vmax limit']
            self.Vmin_lim = features['Vmin limit']
            self.sigma = features['sigma']
            self.sigma_count = features['sigma_count']
            self.sigma_max = 1
            
    def apply_diode_limiting_upper_bound_only(self, V_new, V_prior, d = 0.95):
        # Inspired by SUGAR Transmission Code Diode Limiting
        cs_tol = 1e-6
        
        high_x_index = np.where(V_prior[self.x_index] > self.xmax)[0]
        V_prior[self.x_index][high_x_index] = self.xmax[high_x_index] - cs_tol
        
        x_k = V_prior[self.x_index]
        x_k1 = V_new[self.x_index]
        delta_x = x_k1 - x_k
        
        mu_max_k = V_prior[self.mu_max_index]
        mu_max_k1 = V_new[self.mu_max_index]
        delta_mu_max = mu_max_k1 - mu_max_k

        pos_delta_mu_max = np.where(delta_mu_max > 0)[0]
        neg_delta_mu_max = np.where(delta_mu_max < 0)[0]
        zero_delta_mu_max = np.where(delta_mu_max == 0)[0]         
        
        neg_delta_x = np.where(delta_x < 0)[0]
        pos_delta_x = np.where(delta_x > 0)[0]
        zero_delta_x = np.where(delta_x == 0)[0] 

        delta_mu_max[zero_delta_mu_max] = self.cs_tol**2 #1e-9
        delta_x[zero_delta_x] = self.cs_tol**2 #1e-9
        
        e_vec = 1e-2*cs_tol*np.ones(len(delta_mu_max)).reshape(-1,1)
        
        amax_mu = np.ones(len(delta_mu_max)).reshape(-1,1)
        amax_mu[neg_delta_mu_max] = d*np.divide((e_vec[neg_delta_mu_max] - mu_max_k[neg_delta_mu_max]), delta_mu_max[neg_delta_mu_max])
         
        if len(self.xmax) == len(x_k):
            amax_x = d*np.divide((self.xmax - x_k), delta_x)
            a_x = np.fmin(1, amax_x)
            a_x[neg_delta_x] = 1
        
        delta_x[zero_delta_x] = 0
        x_lim = (x_k + np.multiply(delta_x, a_x)).reshape(-1,1)
        
        amax_mu = np.fmin(1, amax_mu)
        
        mu_max_lim = mu_max_k + np.multiply(delta_mu_max, amax_mu)
        
        V_lim = np.copy(V_new)
        V_lim[self.x_index] = np.copy(x_lim)
        V_lim[self.mu_max_index] = np.copy(mu_max_lim)
        
        return V_lim
                
    def apply_diode_limiting(self, V_new, V_old):
        #ipdb.set_trace()
        # Inspired by SUGAR Transmission Code Diode Limiting
        cs_tol = 1e-6
        

        if (self.mu_min_index != None):

            # high_x_index = np.where(V_old[self.x_index] > self.xmax)
            low_x_index = np.where(V_new[self.x_index] < self.xmin)[0]
            # V_old[self.x_index][high_x_index] = self.xmax[high_x_index] - 1e-6
            V_new[self.x_index][low_x_index] = self.xmin[low_x_index] + 1e-6
            
            x_k = V_old[self.x_index]
            x_k1 = V_new[self.x_index]
            delta_x = x_k1 - x_k
            
            # mu_max_k = V_old[self.mu_max_index]
            mu_min_k = V_old[self.mu_min_index]
            # mu_max_k1 = V_new[self.mu_max_index]
            mu_min_k1 = V_new[self.mu_min_index]   
            # delta_mu_max = mu_max_k1 - mu_max_k
            delta_mu_min = mu_min_k1 - mu_min_k
            
            # neg_delta_mu_max = np.where(delta_mu_max < 0)
            # zero_delta_mu_max = np.where(delta_mu_max == 0)         
            # delta_mu_max[zero_delta_mu_max] = 1e-9
            
            neg_delta_mu_min = np.where(delta_mu_min < 0)      
            zero_delta_mu_min = np.where(delta_mu_min == 0)
            delta_mu_min[zero_delta_mu_min] = 1e-9
        
            neg_delta_x = np.where(delta_x < 0)
            pos_delta_x = np.where(delta_x > 0)
            zero_delta_x = np.where(delta_x == 0)    
            delta_x[zero_delta_x] = 1e-9
            
            e_vec = 1e-2*cs_tol*np.ones(len(delta_mu_min)).reshape(-1,1)
            
            # amax_mu = np.ones(len(delta_mu_max)).reshape(-1,1)
            # amax_mu[neg_delta_mu_max] = 0.95*np.divide((e_vec[neg_delta_mu_max] - mu_max_k[neg_delta_mu_max]), delta_mu_max[neg_delta_mu_max])
            
            amin_mu = np.ones(len(delta_mu_min)).reshape(-1,1)
            amin_mu[neg_delta_mu_min] = 0.95*np.divide((e_vec[neg_delta_mu_min] - mu_min_k[neg_delta_mu_min]), delta_mu_min[neg_delta_mu_min])   
            
            # if len(self.xmin) == len(x_k) and len(self.xmax) == len(x_k):
            #     amax_x = 0.95*np.divide((self.xmax - x_k), delta_x)
            #     amin_x = 0.95*np.divide((self.xmin - x_k), delta_x)
            #     a_x = np.fmin(1, np.maximum(amin_x, amax_x))
            # elif len(self.xmin) == len(x_k):
            amin_x = 0.95*np.divide((self.xmin - x_k), delta_x)
            a_x = np.fmin(1, amin_x)
            a_x[pos_delta_x] = 1
            # elif len(self.xmax) == len(x_k):
            #     amax_x = 0.95*np.divide((self.xmax - x_k), delta_x)
            #     a_x = np.fmin(1, amax_x)
            #     a_x[neg_delta_x] = 1
            
            x_lim = (x_k + np.multiply(delta_x, a_x)).reshape(-1,1)
            
            # amax_mu = np.fmin(1, amax_mu)
            amin_mu = np.fmin(1, amin_mu)
            
            # mu_max_lim = mu_max_k + np.multiply(delta_mu_max, amax_mu)
            mu_min_lim = mu_min_k + np.multiply(delta_mu_min, amin_mu)
            
            V_lim = np.copy(V_new)
            V_lim[self.x_index] = np.copy(x_lim)
            # V_lim[self.mu_max_index] = np.copy(mu_max_lim)
            V_lim[self.mu_min_index] = np.copy(mu_min_lim)


        if (self.mu_max_index != None):
            print("performing maximum diode limiting")
            # Inspired by SUGAR Transmission Code Diode Limiting
            cs_tol = 1e-6
            
            high_x_index = np.where(V_old[self.x_index] > self.xmax)
            V_old[self.x_index][high_x_index] = self.xmax[high_x_index] - 1e-6
            
            x_k = V_old[self.x_index]
            x_k1 = V_new[self.x_index]
            delta_x = x_k1 - x_k
            
            mu_max_k = V_old[self.mu_max_index]
            mu_max_k1 = V_new[self.mu_max_index] 
            delta_mu_max = mu_max_k1 - mu_max_k
            
            neg_delta_mu_max = np.where(delta_mu_max < 0)
            zero_delta_mu_max = np.where(delta_mu_max == 0)         
            delta_mu_max[zero_delta_mu_max] = 1e-9
            
        
            neg_delta_x = np.where(delta_x < 0)
            pos_delta_x = np.where(delta_x > 0)
            zero_delta_x = np.where(delta_x == 0)    
            delta_x[zero_delta_x] = 1e-9
            
            e_vec = 1e-2*cs_tol*np.ones(len(delta_mu_max)).reshape(-1,1)
            
            amax_mu = np.ones(len(delta_mu_max)).reshape(-1,1)
            amax_mu[neg_delta_mu_max] = 0.95*np.divide((e_vec[neg_delta_mu_max] - mu_max_k[neg_delta_mu_max]), delta_mu_max[neg_delta_mu_max])
            
            
            amax_x = 0.95*np.divide((self.xmax - x_k), delta_x)
            a_x = np.fmin(1, amax_x)
            a_x[neg_delta_x] = 1
            
            x_lim = (x_k + np.multiply(delta_x, a_x)).reshape(-1,1)
            
            amax_mu = np.fmin(1, amax_mu)
            
            mu_max_lim = mu_max_k + np.multiply(delta_mu_max, amax_mu)
            
            V_lim = np.copy(V_new)
            V_lim[self.x_index] = np.copy(x_lim)

            V_lim[self.mu_max_index] = np.copy(mu_max_lim)

        
        return V_lim

    def apply_voltage_limiting(self, Vsol, V):
        # Apply voltage limiting module here
        deltaVr = Vsol[Nodes.Vr_index] - V[Nodes.Vr_index]
        deltaVi = Vsol[Nodes.Vi_index] - V[Nodes.Vi_index]
        # Convert the norm vectors into array
        if self.Vr_norm is None:
            Vr_norm = np.array(Nodes.Vr_mag, dtype=float, ndmin=2).T
        if self.Vi_norm is None:
            Vi_norm = np.array(Nodes.Vi_mag, dtype=float, ndmin=2).T
    
        # Normalize the solution deltaVr and deltaVi vectors
        delta_Vr_norm = np.divide(deltaVr, Vr_norm)
        delta_Vi_norm = np.divide(deltaVi, Vi_norm)
    
        # Min-max of the deltaVr and deltaVi
        delta_Vr_norm = np.minimum(delta_Vr_norm, self.maxStep__)
        delta_Vr_norm = np.maximum(delta_Vr_norm, self.minStep__)
        delta_Vi_norm = np.minimum(delta_Vi_norm, self.maxStep__)
        delta_Vi_norm = np.maximum(delta_Vi_norm, self.minStep__)
    
        # Convert back to deltaVr
        delta_Vr = np.multiply(delta_Vr_norm, Vr_norm)
        delta_Vi = np.multiply(delta_Vi_norm, Vi_norm)
    
        # Update V
        V[Nodes.Vr_index] += delta_Vr
        V[Nodes.Vi_index] += delta_Vi
    
    
        if self.stamp_dual:
            dual_max_step = 2*self.maxStep__
            dual_min_step = 2*self.minStep__
            deltaL = Vsol[Nodes.L_index] - V[Nodes.L_index]
            np.fmin(deltaL, dual_min_step * np.ones((np.shape(deltaL)[0], 1)), out=deltaL)
            np.fmax(deltaL, dual_max_step * np.ones((np.shape(deltaL)[0], 1)), out=deltaL)    
            V[Nodes.L_index] += deltaL
            

        return V       

    def apply_variable_limiting(self, Vsol, errStore, nodes):
       
        def find_max_min(Vsol, Vnom, Nodes):
            V = []
        
            V = abs(Vsol[Nodes.Vr_index] + 1j * Vsol[Nodes.Vi_index])
            Vr_mag = np.array(Nodes.Vr_mag, dtype=float, ndmin=2).T
        
            Vsol_pu = np.divide(V, Vr_mag)
        
            # find the maximum node voltage and its location
            Vmax = np.amax(Vsol_pu)
            loc_max = np.argmax(Vsol_pu)
        
            # find the minimum node voltage and its location
            Vmin = np.amin(Vsol_pu)
            loc_min = np.argmin(Vsol_pu)
        
            return Vmax, Vmin, loc_max, loc_min
    
    
        # update nominal voltages/ setpoints
        Vnom = []
    
        for node in nodes:
            if not node.isTriplex:
                node.update_Vnom(Vsol)
            Vnom.append(node.Vnom)
    
        # find maximum and minimum voltages in the solution vector
        Vmax, Vmin, loc_max, loc_min = find_max_min(Vsol, Vnom, Nodes)
    
        # check if the step is too large
        if Vmax > self.Vmax_lim and self.sigma > 1e-2:
            self.sigma *= 0.5
        elif Vmin < self.Vmin_lim and self.sigma > 1e-2:
            self.sigma *= 0.5
    
        self.sigma_count += 1
    
        # if variable limiting applied at least once, check if error is monotonically decreasing
        # if monotonically decreasing and sigma is greater than 1 then scale back up
        if errStore[-1] < errStore[-2]:
            self.sigma_count += 1
        else:
            self.sigma_count = 0
    
        if self.sigma_count > 5:
            self.sigma /= 0.75
            min(self.sigma, self.sigma_max)
        
        