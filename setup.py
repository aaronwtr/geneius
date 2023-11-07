from setuptools import setup
import io

setup(
    name='geneius',
    version='1.2',
    packages=['src'],
    url='https://github.com/aaronwtr/geneius',
    description='A tool for disease-gene evidence search and explanation',
    long_description=io.open('README.md', 'r', encoding='utf-8').read(),
    license='MIT',
    author='Aaron Wenteler',
    email='a.wenteler@qmul.ac.uk',
    install_requires=[
        'anthropic',
        'biopython',
    ],
    entry_points={
        'console_scripts': [
            'geneius = src.geneius:main'
        ]
    }
)
