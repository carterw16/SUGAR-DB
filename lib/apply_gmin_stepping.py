from classes.Loads import Constant
from classes.TriplexLoads import TriplexLoads


def apply_gmin_stepping(node, triplexload, load, gmin_load, connected_nodes,
                        gmin_stepping):
    # Stamp High Impedance Conductance to ground for floating nodes
    gmin_load_count = 0
    for ele in node:
        _gmin_val = []
        if gmin_stepping:
            # Check the phases
            _gmin_val.append(gmin_load) if ele.phases & 0x01 == int(
                0x01) else _gmin_val.append(None)
            _gmin_val.append(gmin_load) if ele.phases & 0x02 == int(
                0x02) else _gmin_val.append(None)
            _gmin_val.append(gmin_load) if ele.phases & 0x04 == int(
                0x04) else _gmin_val.append(None)
            # Connect gmin to floating nodes
            if ele.degree == 1:
                if ele.ID not in connected_nodes and ele.bustype != 3:
                    if not ele.isTriplex:
                        for ld in load:
                            if ld.Type == "Constant":
                                load.append(
                                    Constant(ele.ID,
                                             ['Gmin load :' + str(ele.name)],
                                             None, ele.phases, None, None, None,
                                             None, None, None, None, None, None,
                                             None, None, None, None,
                                             _gmin_val[0], _gmin_val[1],
                                             _gmin_val[2]))
                            else:
                                load.append(
                                    ZIPLoad(ele.ID,
                                            ['Gmin load :' + str(ele.name)],
                                            None, ele.phases, None, None, None,
                                            None, None, None, None, None, None,
                                            None, None, None, None,
                                            _gmin_val[0], _gmin_val[1],
                                            _gmin_val[2]))
                            gmin_load_count += 1
                    else:
                        triplexload.append(
                            TriplexLoads(
                                ['Gmin Triplex Load :' + str(ele.name)], ele.ID,
                                None, None, None, None, None, gmin_load, None,
                                gmin_load, None, None, ele.nominal_voltage,
                                None, ele.phases))
                        gmin_load_count += 1
    return load, triplexload, gmin_load_count
