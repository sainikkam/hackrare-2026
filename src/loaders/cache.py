"""
File-based JSON cache. Each loader function is decorated with @file_cache(namespace).
Cache miss -> API call -> write to disk -> return.
Cache hit -> read from disk -> return. Prevents redundant API calls on re-runs.
"""

import json
import hashlib
import functools
import logging
from pathlib import Path

CACHE_DIR = Path(__file__).resolve().parents[2] / "data" / "cache"
logger = logging.getLogger(__name__)


def file_cache(namespace: str):
    """Decorator factory. namespace is used as a subdirectory under data/cache/."""
    def decorator(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            key_data = json.dumps({"args": args, "kwargs": kwargs}, sort_keys=True, default=str)
            key = hashlib.md5(key_data.encode()).hexdigest()[:12]
            cache_path = CACHE_DIR / namespace / f"{key}.json"

            if cache_path.exists():
                logger.debug("Cache hit: %s/%s", namespace, key)
                return json.loads(cache_path.read_text())

            logger.debug("Cache miss: %s/%s — calling API", namespace, key)
            result = fn(*args, **kwargs)

            if result is not None:
                cache_path.parent.mkdir(parents=True, exist_ok=True)
                cache_path.write_text(json.dumps(result, indent=2))

            return result
        return wrapper
    return decorator


def clear_cache(namespace: str = None):
    """Clear all cache, or just one namespace."""
    target = CACHE_DIR / namespace if namespace else CACHE_DIR
    if target.exists():
        import shutil
        shutil.rmtree(target)
        target.mkdir(parents=True, exist_ok=True)
