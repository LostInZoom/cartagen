import sys
from setuptools import setup

if sys.version_info[:2] < (3, 8):
    error = ('cartagen4py requires Python 3.8 or later (%d.%d detected).' % sys.version_info[:2])
    sys.stderr.write(error + "\n")
    sys.exit(1)

# General informations
name = 'cartagen4py'
version = '0.1.2'
description = 'Python package to generalise geographic objects for cartographic purposes'
url = 'https://github.com/LostInZoom/cartagen4py'
author = 'Guillaume Touya'
author_email = 'guillaume.touya@ign.fr'
lic = 'GPLv3'
packages = [
    'cartagen4py',
    'cartagen4py.algorithms',
    'cartagen4py.algorithms.buildings'
]

# Requirements and dependencies
python_requires = '>=3.8'
install_requires = [
    'numpy',
    'shapely',
    'geopandas',
]

# Meta informations
keywords = [
    'Generalization',
    'Cartography',
    'cartographic generalization',
]
platforms = ['Linux', 'Mac OSX', 'Windows', 'Unix']
classifiers = [
    'Development Status :: 1 - Planning',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
    'Programming Language :: Python :: 3 :: Only',
]

if __name__ == '__main__':
    setup(
        name=name,
        version=version,    
        description=description,
        url=url,
        author=author,
        author_email=author_email,
        license=lic,
        packages=packages,
        python_requires=python_requires,
        install_requires=install_requires,
        keywords = keywords,
        platforms = platforms,
        classifiers=classifiers,
    )
