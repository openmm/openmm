from simtk.openmm import XmlSerializer
from simtk.openmm.openmm import State, System, Integrator


class XMLFile(object):

    def __new__(cls, filename):
        obj = cls._load_xml_(filename)
        assert isinstance(obj, cls.xmlcls), "File does not contain %s." % cls.xmlcls
        return obj

    @staticmethod
    def _load_xml_(filename):
        with open(filename, 'rb') as f:
            xmlfile = XmlSerializer.deserialize(f.read())
        return xmlfile


class XMLStateFile(XMLFile):
    xmlcls = State


class XMLSystemFile(XMLFile):
    xmlcls = System


class XMLIntegratorFile(XMLFile):
    xmlcls = Integrator
