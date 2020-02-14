from setuptools import setup

setup(
    name = "pyzfs",
    version = "0.0.0",
    author = "He Ma, Marco Govoni",
    author_email = "mahe@uchicago.edu, mgovoni@anl.gov",
    description = "A python code to compute zero-field splitting tensor",
    packages=["pyzfs"],
    classifiers=[],
    install_requires=[
        'numpy', 'scipy', 'mpi4py', 'ase', 'lxml',
    ]
)

