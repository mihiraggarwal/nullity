# nullity

A Python package for all things linear algebra, built in vanilla python3  
Built during MAT-1001 / Linear Algebra at Ashoka University.

## Installation

The package is available on PyPI:

```ps
pip install nullity
```

To start using, simply import the package in your Python file:

```py
from nullity import Matrix
```

## Usage

**Instantiation**  

```py
m = Matrix(m, n, *args)
```

The matrix must be instantiated by first passing in the number of rows (m), then the number of columns (n), and finally the m*n values of the matrix separated by commas.  

These values are arranged first in rows and then in columns.  
For example:

```py
>>> m = Matrix(2, 3, 1, 2, 3, 4, 5, 6)
>>> print(m)

      1.00       2.00       3.00
      4.00       5.00       6.00

```

**Methods**  

| Function | Output |
| --- | --- |
| `nrows` | Number of rows |
| `ncolumns` | Number of columns |
| `is_square` | Whether the matrix is square |
| `rref` | RREF of the matrix |
| `rank` | Rank of the matrix |
| `nullity` | Nullity (dimension of the null space) of the matrix |
| `is_invertible` | Whether the matrix is invertible |
| `inverse` | Inverse of the invertible matrix |
| `rank_factorization` | The two rank factorized matrices `r` and `c` |
| `row_basis` | Basis of the row space |
| `col_basis` | Basis of the column space |
| `null_basis` | Basis of the null space |
| `transpose` | Transpose of the matrix |
| `plu` | PLU decomposed matrices `p`, `l`, and `u` |
| `det` | Determinant of a square matrix |
| `qr` | QR decomposed matrices `q` and `r` |
| `charpol` | Characteristic polynomial of the matrix |
| `eigenvals` | Eigenvalues of the matrix (from numpy) |
| `eigenvecs` | Eigenvectors corresponding to the respective eigenvalues |
| `svd` | Singular value decomposed matrices `u`, `sigma`, and `v_transpose` |

Along with these, addition and multiplication of 2 matrices is also supprted using the inbuilt operators `+` and `*` respectively.

## Issues

Python's decimal arithmetic often results in incorrect answers. In some cases, approximated output might be sufficient, however in others, errors may be raised.
