import os
from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'A package for all things linear algebra'

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, "README.md")) as f:
    long_description = f.read()

setup(
    name="lin_alg",
    version=VERSION,
    author="Mihir Aggarwal",
    author_email="mail@mihiraggarwal.me",
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=['numpy'],
    keywords=['linear', 'algebra', 'linear algebra', 'linalg', 'lin-alg', 'linearalgebra', 'eigen', 'eigenvalue', 'eigenvector'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ]
)