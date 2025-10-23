"""
Minimal scimap module to prevent import errors.
This module provides just enough interface to prevent the spac package from failing.
"""

__version__ = "2.1.3"

# Minimal stub class that can be called as a function or accessed as an attribute
class _SciMapStub:
    def __init__(self, *args, **kwargs):
        pass
    
    def __call__(self, *args, **kwargs):
        return self
    
    def __getattr__(self, name):
        return _SciMapStub()

# Create stub submodules and functions that spac might try to use
tl = _SciMapStub()
pl = _SciMapStub() 
pp = _SciMapStub()
hl = _SciMapStub()

# Add any specific functions that might be called
def phenotype(*args, **kwargs):
    return _SciMapStub()

def cluster(*args, **kwargs):
    return _SciMapStub()

def spatial(*args, **kwargs):
    return _SciMapStub()

# Make this module callable and attribute-accessible
def __getattr__(name):
    return _SciMapStub()