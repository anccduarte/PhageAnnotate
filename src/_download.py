# -*- coding: utf-8 -*-

if __name__ == "__main__":
    
    import argparse
    from download_records import DownloadRecords
    from pathlib import Path
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-database", default="nucleotide")
    parser.add_argument("-base_dir", default="../records")
    parser.add_argument("-taxid")
    parser.add_argument("-num_ids")
    
    args = parser.parse_args()
    database = args.database
    base_dir = args.base_dir
    taxid = args.taxid # 28883 OR 2731619
    num_ids = args.num_ids
    
    Path(base_dir).mkdir(exist_ok=True)
    
    if taxid is None:
        raise ValueError("<taxid> and <num_ids> have no default values.")
        
    num_recs = DownloadRecords(database=database,
                               base_dir=base_dir,
                               taxid=taxid,
                               num_ids=int(num_ids)).download(max_tries=5)
    print(f"Number of records downloaded: {num_recs}")
