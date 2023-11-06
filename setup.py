from setuptools import setup, find_packages

setup(
    name='geneius',
    version='1.0',
    packages=find_packages(),
    install_requires=[
        'anthropic',
        'biopython',
    ],
    entry_points={
        'console_scripts': [
            'geneius = geneius.__main__:main'
        ]
    }
)
