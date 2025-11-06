# assay_design/cache_manager.py

import os
import json
import logging
import hashlib
import datetime
from typing import List, Dict, Any, Optional
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

class SequenceCacheManager:
    """Manages caching of sequence data to avoid redundant downloads."""
    
    def __init__(self, cache_dir: str = None):
        """
        Initialize the cache manager.
        
        Args:
            cache_dir: Directory to store cached sequences (default: ~/.assay_design/cache)
        """
        if cache_dir is None:
            # Default cache directory in user's home folder
            home_dir = os.path.expanduser("~")
            cache_dir = os.path.join(home_dir, ".assay_design", "cache")
            
        self.cache_dir = cache_dir
        self.sequence_dir = os.path.join(cache_dir, "sequences")
        self.metadata_dir = os.path.join(cache_dir, "metadata")
        
        # Create cache directories if they don't exist
        os.makedirs(self.sequence_dir, exist_ok=True)
        os.makedirs(self.metadata_dir, exist_ok=True)
        
        logger.info(f"Initialized sequence cache at {cache_dir}")
    
    def _generate_cache_key(self, params: Dict[str, Any]) -> str:
        """
        Generate a unique cache key from parameters.
        
        Args:
            params: Dictionary of parameters that define the sequence query
            
        Returns:
            Unique hash key for the cache entry
        """
        # Create a sorted string representation of the parameters
        param_str = json.dumps(params, sort_keys=True)
        
        # Hash the parameters
        key = hashlib.md5(param_str.encode()).hexdigest()
        return key
    
    def get_cached_sequences(self, 
                           query_type: str, 
                           params: Dict[str, Any],
                           max_age_days: int = 30) -> Optional[List[SeqRecord]]:
        """
        Get sequences from cache if available and not expired.
        
        Args:
            query_type: Type of query ('taxid', 'gene', 'cds', etc.)
            params: Parameters used for the original query
            max_age_days: Maximum age of cached data in days
            
        Returns:
            List of cached sequences or None if not available/expired
        """
        # Add query_type to parameters for unique key generation
        cache_params = params.copy()
        cache_params['query_type'] = query_type
        
        # Generate cache key
        cache_key = self._generate_cache_key(cache_params)
        
        # Check if metadata exists
        metadata_path = os.path.join(self.metadata_dir, f"{cache_key}.json")
        if not os.path.exists(metadata_path):
            logger.debug(f"No cache metadata found for key {cache_key}")
            return None
        
        # Load metadata
        try:
            with open(metadata_path, 'r') as f:
                metadata = json.load(f)
                
            # Check cache freshness
            timestamp = metadata.get('timestamp')
            cache_date = datetime.datetime.fromisoformat(timestamp)
            age = datetime.datetime.now() - cache_date
            
            if age.days > max_age_days:
                logger.info(f"Cache expired (age: {age.days} days, max: {max_age_days})")
                return None
                
            # Check if sequence file exists
            sequence_path = os.path.join(self.sequence_dir, f"{cache_key}.fasta")
            if not os.path.exists(sequence_path):
                logger.warning(f"Cache metadata exists but sequence file missing for {cache_key}")
                return None
                
            # Load sequences
            try:
                sequences = list(SeqIO.parse(sequence_path, "fasta"))
                logger.info(f"Loaded {len(sequences)} sequences from cache (key: {cache_key})")
                return sequences
            except Exception as e:
                logger.error(f"Error loading cached sequences: {str(e)}")
                return None
                
        except Exception as e:
            logger.error(f"Error reading cache metadata: {str(e)}")
            return None
    
    def cache_sequences(self, 
                      query_type: str, 
                      params: Dict[str, Any], 
                      sequences: List[SeqRecord]) -> str:
        """
        Save sequences to cache.
        
        Args:
            query_type: Type of query ('taxid', 'gene', 'cds', etc.)
            params: Parameters used for the original query
            sequences: Sequence data to cache
            
        Returns:
            Cache key used for storage
        """
        if not sequences:
            logger.warning("Attempted to cache empty sequence list")
            return ""
            
        # Add query_type to parameters for unique key generation
        cache_params = params.copy()
        cache_params['query_type'] = query_type
        
        # Generate cache key
        cache_key = self._generate_cache_key(cache_params)
        
        # Save sequences to FASTA file
        sequence_path = os.path.join(self.sequence_dir, f"{cache_key}.fasta")
        try:
            SeqIO.write(sequences, sequence_path, "fasta")
            
            # Create metadata
            metadata = {
                'query_type': query_type,
                'params': params,
                'timestamp': datetime.datetime.now().isoformat(),
                'count': len(sequences),
                'source': 'NCBI'
            }
            
            # Save metadata
            metadata_path = os.path.join(self.metadata_dir, f"{cache_key}.json")
            with open(metadata_path, 'w') as f:
                json.dump(metadata, f, indent=2)
                
            logger.info(f"Cached {len(sequences)} sequences with key {cache_key}")
            return cache_key
            
        except Exception as e:
            logger.error(f"Error caching sequences: {str(e)}")
            return ""
    
    def clear_cache(self, older_than_days: Optional[int] = None) -> int:
        """
        Clear the cache, optionally only entries older than specified days.
        
        Args:
            older_than_days: If provided, only clear entries older than this many days
            
        Returns:
            Number of cache entries cleared
        """
        cleared_count = 0
        
        # If no age specified, clear all
        if older_than_days is None:
            # Get all files in metadata directory
            for filename in os.listdir(self.metadata_dir):
                if filename.endswith('.json'):
                    cache_key = filename.split('.')[0]
                    
                    # Remove metadata file
                    os.remove(os.path.join(self.metadata_dir, filename))
                    
                    # Remove sequence file if it exists
                    sequence_path = os.path.join(self.sequence_dir, f"{cache_key}.fasta")
                    if os.path.exists(sequence_path):
                        os.remove(sequence_path)
                        
                    cleared_count += 1
        else:
            # Clear only old entries
            cutoff_date = datetime.datetime.now() - datetime.timedelta(days=older_than_days)
            
            for filename in os.listdir(self.metadata_dir):
                if filename.endswith('.json'):
                    metadata_path = os.path.join(self.metadata_dir, filename)
                    cache_key = filename.split('.')[0]
                    
                    try:
                        with open(metadata_path, 'r') as f:
                            metadata = json.load(f)
                            
                        timestamp = metadata.get('timestamp')
                        cache_date = datetime.datetime.fromisoformat(timestamp)
                        
                        if cache_date < cutoff_date:
                            # Remove metadata file
                            os.remove(metadata_path)
                            
                            # Remove sequence file if it exists
                            sequence_path = os.path.join(self.sequence_dir, f"{cache_key}.fasta")
                            if os.path.exists(sequence_path):
                                os.remove(sequence_path)
                                
                            cleared_count += 1
                    except Exception as e:
                        logger.error(f"Error processing metadata file {filename}: {str(e)}")
                        
        logger.info(f"Cleared {cleared_count} cache entries")
        return cleared_count
    
    def get_cache_info(self) -> Dict[str, Any]:
        """
        Get information about the current cache state.
        
        Returns:
            Dictionary with cache statistics
        """
        metadata_files = [f for f in os.listdir(self.metadata_dir) if f.endswith('.json')]
        sequence_files = [f for f in os.listdir(self.sequence_dir) if f.endswith('.fasta')]
        
        # Calculate total size
        total_size = 0
        for filename in sequence_files:
            file_path = os.path.join(self.sequence_dir, filename)
            total_size += os.path.getsize(file_path)
            
        # Get age information
        ages = []
        for filename in metadata_files:
            try:
                with open(os.path.join(self.metadata_dir, filename), 'r') as f:
                    metadata = json.load(f)
                timestamp = metadata.get('timestamp')
                cache_date = datetime.datetime.fromisoformat(timestamp)
                age = (datetime.datetime.now() - cache_date).days
                ages.append(age)
            except:
                continue
                
        # Compile statistics
        stats = {
            'entry_count': len(metadata_files),
            'total_size_mb': total_size / (1024 * 1024),
            'oldest_entry_days': max(ages) if ages else 0,
            'newest_entry_days': min(ages) if ages else 0,
            'average_age_days': sum(ages) / len(ages) if ages else 0,
            'cache_dir': self.cache_dir
        }
        
        return stats