import os


def get_root_directory():
    if os.path.basename((os.path.abspath(os.getcwd()))) != 'SUGAR-DB':
        root_dir = os.path.abspath('.').replace(
            '{sep}{dirname}'.format(dirname=os.path.basename(
                os.path.abspath(os.getcwd())),
                                    sep=os.path.sep), '')
    else:
        root_dir = os.path.abspath('.')

    return root_dir
