from setuptools import setup, find_packages

setup(
    name='geneius',
    version='1.4',
    packages=['src'],
    url='https://github.com/aaronwtr/geneius',
    description='A tool for disease-gene evidence search and explanation',
    long_description='More information on the GitHub repository: https://github.com/aaronwtr/geneius',
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
