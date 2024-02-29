from ditto.readers.windmil.read import Reader
from ditto.store import Store
from ditto.writers.opendss.write import Writer

m = Store()
windmil_reader = Reader()
windmil_reader.parse(m)
OpenDSS_writer = Writer()
OpenDSS_writer.write(m)
