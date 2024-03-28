from ditto.readers.windmil_ascii.read import Reader
from ditto.store import Store
from ditto.writers.opendss.write import Writer

m = Store()

Path = {'network_folder': r'C:\Users\alatif\Desktop\Test-Exports\ASCII\Basalt'}

windmil_reader = Reader(**Path)
windmil_reader.parse(m)
OpenDSS_writer = Writer(
    output_path=r'C:\Users\alatif\Desktop\Test-Exports\OpenDSS\Basalt')
OpenDSS_writer.write(m)
