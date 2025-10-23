"""
Memory registry for creating virtual file paths for in-memory objects.
This allows run_from_json to work with in-memory AnnData objects.
"""
import uuid

from spac.templates import template_utils

# Global registry for in-memory objects
_MEMORY_REGISTRY = {}


def register_memory_object(obj) -> str:
    """
    Register an in-memory object and return a virtual path.
    
    Parameters
    ----------
    obj : Any
        The object to register (typically AnnData)
        
    Returns
    -------
    str
        Virtual path that can be used with load_input
    """
    # Create unique identifier
    virtual_path = f"memory://{uuid.uuid4().hex}"
    _MEMORY_REGISTRY[virtual_path] = obj
    return virtual_path


def unregister_memory_object(virtual_path: str):
    """Remove object from memory registry."""
    _MEMORY_REGISTRY.pop(virtual_path, None)


def get_memory_object(virtual_path: str):
    """Get object from memory registry."""
    return _MEMORY_REGISTRY.get(virtual_path)


# Patch load_input to handle memory:// paths
_original_load_input = template_utils.load_input


def patched_load_input(file_path):
    """Enhanced load_input that handles memory:// virtual paths."""
    if isinstance(file_path, str) and file_path.startswith("memory://"):
        obj = get_memory_object(file_path)
        if obj is not None:
            return obj
        else:
            raise FileNotFoundError(f"Memory object not found: {file_path}")
    else:
        # Use original function for real file paths
        return _original_load_input(file_path)


# Apply the patch globally
template_utils.load_input = patched_load_input