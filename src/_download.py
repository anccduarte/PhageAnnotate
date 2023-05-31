# -*- coding: utf-8 -*-

if __name__ == "__main__":
    
    import utils
    from download_records import DownloadRecords
    from pathlib import Path
        
    args = utils.get_args(("-database", "nucleotide"),
                          ("-base_dir", "../records"),
                          ("-taxid",),
                          ("-num_ids",))
    
    database = args.database
    base_dir = args.base_dir
    taxid = args.taxid # 28883 OR 2731619
    num_ids = args.num_ids
    
    Path(base_dir).mkdir(exist_ok=True)
    
    if taxid is None or num_ids is None:
        raise ValueError("'taxid' and 'num_ids' have no default values. Please do:\n"
                         ">>> python _download.py -taxid <taxid> -num_ids <num_ids>")
        
    num_recs = DownloadRecords(database=database,
                               base_dir=base_dir,
                               taxid=taxid,
                               num_ids=int(num_ids)).download(max_tries=5)
    
    print(f"Number of records downloaded: {num_recs}")
