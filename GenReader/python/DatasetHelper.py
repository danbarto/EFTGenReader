import os
import subprocess
import json

# Helper class for managing collections of datasets and getting the corresponding root files
class DatasetHelper(object):
    def __init__(self,fn=None):
        self.__datasets = {}
        self.__json_files = []
        self.root_redirect = "root://ndcms.crc.nd.edu/"
        self.hadoop_protocol = "file://"
        if fn: self.load(fn)

    # Load the dataset info from a json file
    def load(self,fn,update=False):
        # update: If set to false, will remove any pre-existing datasets and only include ones
        #         found in the specified file
        if not os.path.exists(fn):
            print("ERROR: Unknown file/directory: {fn}".format(fn=fn))
            return
        if not update:
            self.__datasets = {}
            self.__json_files = []
        print("Loading JSON: {fn}".format(fn=fn))
        self.__json_files.append(fn)
        with open(fn) as f:
            d = json.load(f)
            for k,ds in d.items():
                self.updateDataset(k,**ds)

    # Save the dataset info to a json file
    def save(self,fn):
        print("Saving JSON: {fn}".format(fn=fn))
        with open(fn,'w') as f:
            d = self.toDict()
            json.dump(d,f,indent=2,sort_keys=True)

    # Return a list of all dataset names currently loaded
    def list(self):
        return list(self.__datasets.keys())

    # Return if the dataset name is known
    def exists(self,name):
        return name in self.__datasets

    # Print a dataset as formatted json
    def dump(self,name):
        if not self.exists(name):
            return
        ds = self.__datasets[name]
        print(json.dumps(ds.toDict(),indent=2,sort_keys=True))

    # Remove all datasets
    def clear(self):
        self.__datasets.clear()

    # Create a new dataset entry, only overwrites if force=True
    def newDataset(self,name,force=False,**kwargs):
        if force:
            self.updateDataset(name,**kwargs)
        else:
            if self.exists(name):
                print("ERROR: Skipping {name} since it already exists!".format(name=name))
                return
            self.__datasets[name] = DSContainer(name,**kwargs)

    # Modify a dataset (creates one if none exists)
    def updateDataset(self,name,**kwargs):
        if not self.exists(name):
            self.__datasets[name] = DSContainer(name)
        self.__datasets[name].setData(**kwargs)

    # Remove a specific dataset
    def removeDataset(self,name):
        if self.exists(name):
            del self.__datasets[name]

    # Renames a dataset entry
    def renameDataset(self,old,new):
        if not self.exists(old):
            return
        self.__datasets[new] = self.__datasets.pop(old)

    # Returns a list of root file locations corresponding to the specified dataset
    def getFiles(self,name):
        if not self.exists(name):
            return []
        ds = self.__datasets[name]
        if ds.getData('on_das'):
            lst = self.findOnDAS(ds)
        else:
            lst = self.findOnHadoop(ds)
        return lst

    # Passthrough to underlying dataset getter
    def getData(self,name,data_field):
        if not self.exists(name):
            return None
        return self.__datasets[name].getData(data_field)

    # Try and find dataset root files from DAS
    # NOTE: Requires access to the 'dasgoclient' script, which should be available
    #       in any CMSSW release after doing 'cmsenv'
    def findOnDAS(self,ds):
        lst = []
        das_query = 'file dataset={ds} | grep file.name, file.size, file.nevents'.format(ds=ds.getData('loc'))
        ret = subprocess.check_output(['dasgoclient','--query',das_query])
        ret = ret.split('\n')
        max_files = 10      # currently arbitrary
        for idx,l in enumerate(ret):
            if idx >= max_files: break
            fn,sz,evts = l.split()
            lst.append(self.root_redirect + fn)
        return lst

    # Try and find dataset root files on the local hadoop cluster
    def findOnHadoop(self,ds):
        lst = []
        loc = str(ds.getData('loc'))    # converts from a unicode string object to a normal str
        if not os.path.exists(loc):
            print("Unknown fpath: {loc}".format(loc=loc))
            return []
        idx = 0
        max_files = 100     # currently arbitrary
        for fn in os.listdir(loc):
            if not ".root" in fn: continue
            if idx >= max_files: break
            fpath = self.hadoop_protocol + os.path.join(loc,fn)
            lst.append(fpath)
            idx += 1
        return lst

    # Returns all datasets in a format suitable to be passed to json.dump()
    # NOTE: The nested dictionaries will be references to the ones stored in the actual
    #       DSContainer() objects
    def toDict(self):
        d = {}
        for k,ds in self.__datasets.items():
            d[k] = ds.toDict()
        return d

class DSContainer(object):
    def __init__(self,name,**kwargs):
        self.__name = name
        self.__data = {
            'datatier': '',
            'dataset': '',
            'description': '',
            'loc': None,
            'on_das': False,
            'is_eft': False,
            'central_xsec': 0.0,
            'files_per_task': 1
        }
        self.setData(**kwargs)

    # Generic setter
    def setData(self,**kwargs):
        for k,v in kwargs.items():
            if k not in self.__data:
                print("Unknown key in {name}: '{key}'".format(name=self.__name,key=k))
                continue
            self.__data[k] = v

    # Generic getter
    def getData(self,data_field):
        if data_field not in self.__data:
            return None
        return self.__data[data_field]

    # Convert self to dictionary suitable to be passed to json.dump()
    def toDict(self):
        return self.__data