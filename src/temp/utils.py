# -*- coding: utf-8 -*-

import time
from Bio import Entrez, SeqIO, SeqRecord

def get_ids(search: str, num_ids: int, max_tries: int) -> list:
    """
    Retrieves a list of IDs from NCBI's Nucleotide database according to the parameters
    <search> and <num_ids>.
    
    Parameters
    ----------
    search: str
        The query to be performed in NCBI's Nucleotide database
    num_ids: int
        The maximum number of IDs to be retrieved
    max_tries: int
        The maximum number of tries
    """
    # receive warnings in case of excessive usage of the E-utilities
    Entrez.email = "pg45464@alunos.uminho.pt"
    # print NCBI search
    print(f"NCBI search: {search}")
    # main loop -> try searching for the IDs a maximum of <max_tries> times
    for _ in range(max_tries):
        try:
            # idtype: by default, ESearch returns GI numbers in its output
            # retmax: total number of UIDs from the retrieved set
            handle = Entrez.esearch(db="nucleotide", term=search, retmax=num_ids)
        except:
            time.sleep(12)
        else:
            record = Entrez.read(handle)
            handle.close()
            # list of IDs found for the expression <search>
            return record["IdList"]
    # raise Exception whenever the IDs are not successfully collected
    raise Exception("Unsuccessfully tried to collect the IDs associated with the "
                    f"expression {search!r} ({max_tries} tries).")

def read_record(id_: int, max_tries: int) -> SeqRecord:
    """
    Reads the NCBI Nucleotide's record having the ID <id_>.

    Parameters
    ----------
    id_: int
        The ID of the record to be fetched
    max_tries: int
        The maximum number of tries
    """
    # receive warnings in case of excessive usage of the E-utilities
    Entrez.email = "pg45464@alunos.uminho.pt"
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
        except Exception as exc:
            try_ += 1
            e_msgs.append(exc)
            time.sleep(12)
        # read and return the record
        else:
            record = SeqIO.read(handle, format="gb")
            handle.close()
            return record
    # raise Exception whenever the record isn't successfuly fetched
    e_msgs_str = "\n".join(e_msgs)
    raise Exception(f"Unsuccessfully tried to fetch the record with ID {id_} "
                    f"({max_tries} tries).\n{e_msgs_str}")
