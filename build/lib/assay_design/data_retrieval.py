# assay_design/data_retrieval.py

import logging
from Bio import Entrez, SeqIO

# Logger for the module
logger = logging.getLogger(__name__)

def fetch_sequence_by_accession(accession_id: str, email: str, db: str = "nucleotide") --> SeqIO.SeqRecord:
    """
    Fetch a single sequence by accession from NCBI.
    
    Args:
        accession_id (str): NCBI accession ID (e.g., 'NC_000913.3').
        email (str): User's email address (required by NCBI to track usage).
        db (str): NCBI database from which to fetch (default 'nucleotide').
        
    Returns:
        SeqIO.SeqRecord: A Biopython SeqRecord containing the requested sequence.
        
    Raises:
        ValueError: If no email is provided.
        RuntimeError: If retrieval fails or sequence is empty.
    """
    if not email:
        raise ValueError("You must provide a valid email address to compy with NCBI's usage policies.")
        
    # Set NCBI Entrez parameters
    Entrez.email = email
    
    try:
        logger.info(f"Fetching sequence for accession {accession_id} from NCBI {db} database.")
        handle = Entrez.efetch(db=db, id=accession_id, rettype="fasta", retmode="text")
        seq_record = SeqIO.read(handle, "fasta")
        handle.close()
        
        if not seq_record or len(seq_record.seq) == 0:
            raise RuntimeError(f"No sequence data returned for accession {accession_id}")
            
        logger.info(f"Successfully retrieved sequence: {seq_record.id}, length: {len(seq_record.seq)}")
        return seq_record
        
    except Exception as e:
        logger.error(f"Error fetching {accession_id} from NCBI: {str(e)}")
        raise
        
def fetch_sequences_for_taxid(taxid: str, email: str, db: str = "nucleotide", max_records: int = 10):
    """
    Example function: search for sequences by NCBI TaxID, then fetch up to max_records results. (Often you'd like to refine your search query further, e.g., to a certain gene.)
    
    Args:
        taxid (str): NCBI TaxID, e.g., '562' for E. coli.
        email (str): User's email address (required by NCBI).
        db (str): NCBI database, default 'nucleotide'.
        max_records (int): Maximum number of sequences to fetch.
        
    Returns:
        List[SeqIO.SeqRecord]: A list of up to max_records retrieved SeqRecords.
    """
    if not email:
        raise ValueError("Email address is required.")
        
    Entrez.email = email
    
    try:
        # Example: searching for sequences with that taxid in NCBI
        logger.info(f"Searching for up to {max_records} sequences in taxid:{taxid}")
        search_handle = Entrez.esearch(db=db, term=f"txid{taxid}[Organism]", retmax=max_records)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        ids = search_results.get("IdList", [])
        logger.info(f"Found {len(ids)} records. Fetching details...")
        
        records = []
        if len(ids) > 0:
            fetch_handle = Entrez.efetch(db=db, id=",".join(ids), rettype="fasta", retmode="text")
            records = list(SeqIO.parse(fetch_handle, "fasta"))
            fetch_handle.close()
            logger.info(f"Fetched {len(records)} sequences for TaxID {taxid}.")
            
        return records
        
    except Exception as e:
        logger.error(f"Error fetching sequences for TaxID {taxid}: {str(e)}")
        raise
        
def load_local_fasta(fasta_path: str):
    """
    Load sequences from a local FASTA file using Biopython.
    
    Args:
        fasta_path (str): Path to the local FASTA file.
        
    Returns:
        List[SeqIO.SeqRecord]: A list of sequences parsed from the FASTA file.
    """
    logger.info(f"Loading local FASTA file from {fasta_path}")
    try:
        seq_records = list(SeqIO.parse(fasta_path, "fasta"))
        logger.info(f"Loaded {len(seq_records)} sequences from {fasta_path}")
        return seq_records
    except Exception as e:
        logger.error(f"Error loading local FASTA file {fasta_path}: {str(e)}")
        raise