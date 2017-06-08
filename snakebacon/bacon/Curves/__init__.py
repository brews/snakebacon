"""
This is a big hack to keep track of where the Curves directory is for the bacon 
C/C++ scripts.
"""

from os import path as _path

here = _path.abspath(_path.dirname(__file__))
