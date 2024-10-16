from setuptools import setup
import setuptools
import codecs
import os.path

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")
    

setup(
    name='AtomicPotentialTester',
    version=get_version("potential_tester/__init__.py"),    
    description='Test atomic potentials',
    url='https://github.com/andrew551/AtomicPotentialTester',
    author='Andrew Smith',
    author_email='andrew.d.p.smith@gmail.com',
    license='GPL-3.0',
    packages=setuptools.find_packages(),
    install_requires=[
                      'numpy',  
                      'tqdm',                   
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: GPL-3.0 License',  
        'Operating System :: Linux',        
        'Programming Language :: Python :: 3.10',
    ],
)