import sys
import os
import versioneer
from setuptools import setup

LONG_DESCRIPTION_SRC = 'README.rst'

def read(file):
    with open(os.path.abspath(file), 'r', encoding='utf-8') as f:
        return f.read()

setup(
    long_description=read(LONG_DESCRIPTION_SRC),
    long_description_content_type="text/x-rst; charset=UTF-8",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)