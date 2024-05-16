import os
from setuptools import setup, find_packages

VERSION = '0.0.5'
DESCRIPTION = 'A package for all things linear algebra'

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, "README.md")) as f:
    long_description = f.read()

setup(
    name="nullity",
    version=VERSION,
    author="Mihir Aggarwal",
    author_email="mail@mihiraggarwal.me",
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mihiraggarwal/nullity",
    packages=find_packages(),
    install_requires=['numpy'],
    keywords=['linear', 'algebra', 'linear algebra', 'linalg', 'lin-alg', 'linearalgebra', 'eigen', 'eigenvalue', 'eigenvector', 'nullity', 'null', 'null space', 'rank nullity'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ]
)