# -*- coding: utf-8 -*-

import time
from Bio import Entrez, SeqIO, SeqRecord
from tqdm.auto import tqdm

class DownloadRecords:
    
    """
    Downloads a maximum of <num_ids> records from the NCBI's database <database>, using
    the search expression 'txid<taxid>[ORGN]'.
    """
    
    def __init__(self, database: str, base_dir: str, taxid: str, num_ids: int) -> None:
        """
        Initializes an instance of DownloadRecords.
        
        Parameters
        ----------
        database: str
            The NCBI's database from which the records are to be retrieved
        base_dir: str
            The directory where the records are to be saved
        taxid: str
            The txid identifier of the sought-after taxa (as in NCBI's Taxonomy database)
        num_ids: int
            The number of IDs to be inspected in the database
        """
        self.database = database
        self.base_dir = base_dir
        self.taxid = taxid
        self.num_ids = num_ids
    
    def __repr__(self) -> str:
        """
        Returns the string representation of the object.
        """
        class_ = self.__class__.__name__
        return f"{class_}({database!r}, {base_dir!r}, {taxid!r}, {num_ids})"
    
    @staticmethod
    def _try_entrez(e_util: callable, max_tries: int) -> object:
        """
        Tries to get a record using the callable <e_util> a maximum number of <max_tries>
        times.
        
        Parameters
        ----------
        e_util: callable
            An E-utility ready to be called
        max_tries: int
            The maximum number of tries
        """
        # receive warnings in case of excessive usage of the E-utilities
        Entrez.email = "pg45464@alunos.uminho.pt"
        # initialize 'e_msgs'
        e_msgs = []
        # main loop -> try get record a maximum of <max_tries> times
        for t in range(max_tries):
            # try calling 'e_util' and assign to 'handle'
            try:
                handle = e_util()
            # in case of error, save message and sleep for 12 seconds
            except Exception as exc:
                e_msgs.append((t+1, exc))
                time.sleep(60)
            # otherwise, return the handle
            else:
                return handle
        # raise exception if the number of allowed tries is exceeded
        e_msgs_str = "\n".join([f"{i}. {e_msg}" for (i, e_msg) in e_msgs])
        raise Exception("Unsuccessfully tried to get the record a total number of "
                        f"{max_tries} times.\n{e_msgs_str}")
        
    def _get_ids(self, max_tries: int) -> list:
        """
        Retrieves a list of IDs from NCBI's database <self.database> according to the
        parameter <self.taxid>.
        
        Parameters
        ----------
        max_tries: int
            The maximum number of tries
        """
        # construct and display search expression
        search = f"txid{self.taxid}[ORGN]"
        print(f"NCBI search: {search}")
        # construct callable wrapping 'Entrez.esearch'
        e_util = lambda: Entrez.esearch(db=self.database,
                                        term=search,
                                        retmax=self.num_ids)
        # try read record
        handle = DownloadRecords._try_entrez(e_util=e_util, max_tries=max_tries)
        record = Entrez.read(handle)
        handle.close()
        # return list of IDs
        return record["IdList"]
                
    def _read_record(self, id_: int, max_tries: int) -> SeqRecord:
        """
        Reads the record of ID <id_> from the NCBI database <self.database>.
        
        Parameters
        ----------
        id_: int
            The ID of the record to be fetched
        max_tries: int
            The maximum number of tries
        """
        # construct callable wrapping 'Entrez.efetch'
        e_util = lambda: Entrez.efetch(db=self.database,
                                       id=id_,
                                       rettype="gb",
                                       retmode="text")
        # try read record
        handle = DownloadRecords._try_entrez(e_util=e_util, max_tries=max_tries)
        record = SeqIO.read(handle, format="gb")
        handle.close()
        # return SeqRecord object
        return record
               
    def download(self, max_tries: int) -> int:
        """
        Downloads the records. A maximum number of <max_tries> attempts is performed
        when fetching each record. Returns the number of records downloaded.
        
        Parameters
        ----------
        max_tries: int
            The maximum number of tries
        """
        # get IDs from search
        id_list = self._get_ids(max_tries=max_tries)   
        # main loop -> iterate through IDs
        for id_ in tqdm(id_list):
            # read record
            record = self._read_record(id_=id_, max_tries=max_tries)
            # save record to .gb file
            with open(f"{self.base_dir}/record_{id_}.gb", "w") as f:
                SeqIO.write(record, f, "gb")
        # return number of records
        return len(id_list)
