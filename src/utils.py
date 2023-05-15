# -*- coding: utf-8 -*-

import time
from Bio import Entrez, SeqIO, SeqRecord
from urllib.error import HTTPError

def get_ids(search: str, num_ids: int) -> list:
    """
    Retrieves a list of IDs from NCBI's Nucleotide database according to the parameters
    <search> and <num_ids>.
    
    Parameters
    ----------
    search: str
        The query to be performed in NCBI's Nucleotide
    num_ids: int
        The maximum number of IDs to be retrieved
    """
    # print NCBI search
    print(f"NCBI search: {search}")
    # idtype: by default, ESearch returns GI numbers in its output
    # retmax: total number of UIDs from the retrieved set to be shown in the XML output
    with Entrez.esearch(db="nucleotide", term=search, retmax=num_ids) as handle:
        record = Entrez.read(handle)
    # list of IDs found for the expression <search>
    return record["IdList"]

def read_record(id_: int, max_tries: int) -> SeqRecord:
    """
    Reads the NCBI record corresponding to the ID <id_>.

    Parameters
    ----------
    id_: int
        The ID of the record to be fetched
    max_tries: int
        The maximum number of tries
    """
    # initialize local variables (try_ and e_msgs)
    try_ = 0
    e_msgs = []
    # main loop -> try fetching the record a maximum of <max_tries> times
    while try_ < max_tries:
        # try fetching record (try number <try_>)
        try:
            handle = Entrez.efetch(db="nucleotide",
                                   id=id_,
                                   rettype="gb",
                                   retmode="text")
        # update local variables and sleep for 10 seconds
        except HTTPError as err:
            try_ += 1
            e_msgs.append(err)
            time.sleep(10)
        # read and return the record
        else:
            record = SeqIO.read(handle, format="gb")
            handle.close()
            return record
    # raise Exception whenever the record isn't successfuly fetched
    e_msgs_str = "\n".join(e_msgs)
    raise Exception(f"Unsuccessfully tried to fetch the record with id={id_} a total "
                    f"of {max_tries} times.\n{e_msgs_str}")
