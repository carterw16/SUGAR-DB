import pandapower as pn
import ditto.readers.csv.read as djson
import os
from ditto.store import Store

net = pn.networks.create_cigre_network_mv(with_der=False)

pn.plotting.simple_plot(net, respect_switches=False, line_width=1.0, bus_size=1.0, ext_grid_size=1.0, trafo_size=1.0, 
						plot_loads=False, plot_sgens=False, load_size=1.0, sgen_size=1.0, switch_size=2.0, switch_distance=1.0, 
						plot_line_switches=False, scale_size=True, bus_color='b', line_color='grey', trafo_color='k', 
						ext_grid_color='y', switch_color='k', library='igraph', show_plot=True, ax=None)

print("network", net)
print("network methods", dir(net))

pn.to_excel(net, "cigre_testcase/cigre_testcase.xlsx")

#pn.converter.to_mpc(net, "cigre_case.mat", init = "flat")

"""
model = Store()
r = djson.Reader(input_file=os.path.join("cigre_testcase","cigre_testcase.json"))

r.parse_lines(model)


r.parse(model)

print("r", r)

print("reader methods", dir(r))
print("reader vars", vars(r))
"""