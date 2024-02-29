import pandas as pd
import numpy as np

def read_open_dss_results(file_name, V, nodes, node_key):
    @staticmethod
    def sind(x):
        return np.sin(float(x) * np.pi / 180)

    @staticmethod
    def cosd(x):
        return np.cos(float(x) * np.pi / 180)
    
    df = pd.read_csv(file_name)

    for bus_index in range(0, len(df['Bus'])):
        name_of_bus = df['Bus'][bus_index].lower()
        if 'load' in name_of_bus:
            name_of_bus = df['Bus'][bus_index].lower().replace("_","")
        node_index_from_SUGAR = node_key[name_of_bus]

        if df[' Node1'][bus_index] == 1:
            V[nodes[node_index_from_SUGAR].nodeA_Vr] = df[' Magnitude1'][bus_index]*cosd(df[' Angle1'][bus_index])
            V[nodes[node_index_from_SUGAR].nodeA_Vi] = df[' Magnitude1'][bus_index]*sind(df[' Angle1'][bus_index])

            if df[' Node2'][bus_index] == 2:
                V[nodes[node_index_from_SUGAR].nodeB_Vr] = df[' Magnitude2'][bus_index]*cosd(df[' Angle2'][bus_index])
                V[nodes[node_index_from_SUGAR].nodeB_Vi] = df[' Magnitude2'][bus_index]*sind(df[' Angle2'][bus_index])
                
                V[nodes[node_index_from_SUGAR].nodeC_Vr] = df[' Magnitude3'][bus_index]*cosd(df[' Angle3'][bus_index])
                V[nodes[node_index_from_SUGAR].nodeC_Vi] = df[' Magnitude3'][bus_index]*sind(df[' Angle3'][bus_index])
            
            else:
                V[nodes[node_index_from_SUGAR].nodeC_Vr] = df[' Magnitude2'][bus_index]*cosd(df[' Angle2'][bus_index])
                V[nodes[node_index_from_SUGAR].nodeC_Vi] = df[' Magnitude2'][bus_index]*sind(df[' Angle2'][bus_index])
                
                V[nodes[node_index_from_SUGAR].nodeB_Vr] = df[' Magnitude3'][bus_index]*cosd(df[' Angle3'][bus_index])
                V[nodes[node_index_from_SUGAR].nodeB_Vi] = df[' Magnitude3'][bus_index]*sind(df[' Angle3'][bus_index])
                
        elif df[' Node1'][bus_index] == 2:
            V[nodes[node_index_from_SUGAR].nodeB_Vr] = df[' Magnitude1'][bus_index]*cosd(df[' Angle1'][bus_index])
            V[nodes[node_index_from_SUGAR].nodeB_Vi] = df[' Magnitude1'][bus_index]*sind(df[' Angle1'][bus_index])

            V[nodes[node_index_from_SUGAR].nodeC_Vr] = df[' Magnitude2'][bus_index]*cosd(df[' Angle2'][bus_index])
            V[nodes[node_index_from_SUGAR].nodeC_Vi] = df[' Magnitude2'][bus_index]*sind(df[' Angle2'][bus_index])
            
            V[nodes[node_index_from_SUGAR].nodeA_Vr] = df[' Magnitude3'][bus_index]*cosd(df[' Angle3'][bus_index])
            V[nodes[node_index_from_SUGAR].nodeA_Vi] = df[' Magnitude3'][bus_index]*sind(df[' Angle3'][bus_index])
        
        elif df[' Node1'][bus_index] == 3:
            V[nodes[node_index_from_SUGAR].nodeC_Vr] = df[' Magnitude1'][bus_index]*cosd(df[' Angle1'][bus_index])
            V[nodes[node_index_from_SUGAR].nodeC_Vi] = df[' Magnitude1'][bus_index]*sind(df[' Angle1'][bus_index])
        
            V[nodes[node_index_from_SUGAR].nodeA_Vr] = df[' Magnitude2'][bus_index]*cosd(df[' Angle2'][bus_index])
            V[nodes[node_index_from_SUGAR].nodeA_Vi] = df[' Magnitude2'][bus_index]*sind(df[' Angle2'][bus_index])
            
            V[nodes[node_index_from_SUGAR].nodeB_Vr] = df[' Magnitude3'][bus_index]*cosd(df[' Angle3'][bus_index])
            V[nodes[node_index_from_SUGAR].nodeB_Vi] = df[' Magnitude3'][bus_index]*sind(df[' Angle3'][bus_index])
        else:
            print('Error in lib/read_open_dss_results: no initial voltage identified for:', name_of_bus)
    
    return V

