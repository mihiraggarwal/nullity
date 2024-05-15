from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'A package for all things linear algebra'

setup(
    name="lin_alg",
    version=VERSION,
    author="Mihir Aggarwal",
    author_email="mail@mihiraggarwal.me",
    description=DESCRIPTION,
    packages=find_packages(),
    install_requires=['numpy'],
    keywords=['linear', 'algebra', 'linear algebra', 'linalg', 'lin-alg', 'linearalgebra', 'eigen', 'eigenvalue', 'eigenvector'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "License :: GPL-3.0 License"
    ]
)