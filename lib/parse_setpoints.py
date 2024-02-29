import csv


def parse_setpoints(config_file):
    '''
    Parse generator setpoints into a list of dictonaries.
    :param config_file: File with generator settings.
    :return:
    '''
    reader = csv.DictReader(open(config_file, 'r'))
    config_list = []
    for line in reader:
        config_list.append(line)
    return config_list
