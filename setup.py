# -*- coding: utf-8 -*-
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyprimer", # Replace with your own username
    version="0.0.2a0",
    author="Alina Frolova, Michał Kowalski, Witold Wydmański, Paweł Łabaj",
    author_email="fshodan@gmail.com, michal.bozydar.kowalski@gmail.com, wwydmanski@gmail.com, pawel.labaj@gmail.com",
    description="Python library for primer bechnmark and primer design for PCR and LAMP assays",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "certifi==2020.4.5.1",
        "h5py==2.10.0",
        "llvmlite==0.32.1",
        "numba==0.49.1",
        "numpy==1.18.4",
        "pandas==1.0.3",
        "python-dateutil==2.8.1",
        "pytz==2020.1",
        "six==1.14.0",
        "wincertstore==0.2"]
)