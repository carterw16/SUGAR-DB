from classes.Methods.base import Methods

class PowerFlow(Methods):
    def __init__(self, testcase, casedata, features, settings, path_to_output, node_key, node_index_):
        super(Methods, self).__init__(testcase, casedata, features, settings, path_to_output, node_key, node_index_)
