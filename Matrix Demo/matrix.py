class Matrix:
    '''
    Constructor for Matrix class; constructs matrices as arrays of rows rather than as arrays of columns

    Params:
        rows: number of rows
        cols: number of columns
        init_matrix: decides if Matrix should be initialized as a 2D array of 0's with the given dimensions, if True then yes, if False then no
    '''
    def __init__(self, rows: int, cols: int, init_matrix: bool): # constructor
        self.rows = rows
        self.cols = cols
        self.initialized = False

        self._matrix = []

        if init_matrix:
            for i in range(self.rows):
                self._matrix.append([])

                for j in range(self.cols):
                    self._matrix[i].append(0)
            
            self.initialized = True


    '''
    Private method setting specific values in the Matrix, returns Value Error if non-numerical input is given

    Params:
        row: row index of value to set
        col: column index of value to set
    '''
    def _set_matrix_value(self, row, col):
        prompt = (str)(row) + ", " + (str)(col) + ": "

        while True:
            try:
                value = float(input(prompt))
                return value
            except ValueError:
                print("[Value Error] - requires numerical value")


    '''
    Gets user input to fill the values of an uninitialized matrix
    '''
    def populate_matrix(self):
        self._matrix = [] # reset matrix in case it already has values

        for i in range(self.rows):
            self._matrix.append([])

            for j in range(self.cols):
                self._matrix[i].append(self._set_matrix_value(i, j))

        self.initialized = True


    '''
    Sets the values of self to the corresponding values of another matrix as long as their dimensions are the same

    Params:
        other: Matrix whose values self will receive
    '''
    def set_matrix_to_other(self, other: 'Matrix'):
        if self.rows == other.rows and self.cols == other.cols:
            for i in range(self.rows):
                for j in range(self.cols):
                    self._matrix[i][j] = other._matrix[i][j]

            self.initialized = True
        else:
            raise Exception('Incompatible matrix dimensions')
    

    '''
    Creates a String representation of self

    Return: String representation
    '''
    def to_string(self):
        s = ''

        for i in range(self.rows):
            s += str(self._matrix[i]) + '\n'

        return s
    

    '''
    Creates the transpose of self

    Return: transpose
    '''
    def transpose(self):
        transpose = Matrix(self.cols, self.rows, False)

        for i in range(self._matrix[0].__len__()):
            col = []

            for j in range(self._matrix.__len__()):
                col.append(self._matrix[j][i])
            
            transpose._matrix.append(col)

        return transpose


    '''
    Swaps two rows of self (for Gassian Elimination)

    Params:
        row1: index of first row in the swap
        row2: index of second row in the swap
    '''
    def row_swap(self, row1: int, row2: int):
        tmp = self._matrix[row1]
        self._matrix[row1] = self._matrix[row2]
        self._matrix[row2] = tmp


    '''
    Scales a row of self (for Gaussian Elimination)

    Params:
        row: index of row to be scaled
        scalar: factor by which row is scaled
    '''
    def row_scale(self, row: int, scalar: float):
        for i in range(self._matrix[row].__len__()):
            self._matrix[row][i] *= scalar


    '''
    Scales the entries of a row of self then adds them to the corresponding entries of another row, altering the latter row (for Gaussian Elimination)

    Params:
        source: index of row whose entries are scaled
        onto: index of row whose entries are added to/altered
        scalar: factor by which entries in source row are scaled before adding to entries of onto row
    '''
    def row_combination(self, source: int, onto: int, scalar: float):
        for i in range(self._matrix[onto].__len__()):
            self._matrix[onto][i] += scalar * self._matrix[source][i]
    

    '''
    Rounds the entries of self and corrects negative zeros ("-0.0")

    Params:
        decimals: decimal place to which to round
    '''
    def round_matrix_n(self, decimals: int):
        for i in range(self._matrix.__len__()):
            for j in range(self._matrix[i].__len__()):
                if decimals >= 0:
                    self._matrix[i][j] = round(self._matrix[i][j], decimals)

                if abs(self._matrix[i][j]) == 0:
                    self._matrix[i][j] = abs(self._matrix[i][j])


    '''
    Private method used in Gaussian Elimination (used in the methods ref() and rref() ). It checks self to see if the entry in the "pivot" position of a column is a zero. If so then it searchs for entrys below, in the same column, if it finds an element beneath the pivot that is nonzero, it performs a row swap between their respective rows.

    Params:
        pivot: index of pivot position
        col: index of column to check

    Return: True if every entry including and below the pivot is zero, otherwise False
    '''
    def _sub_pivot_zero_check(self, pivot: int, col: int):
        all_zeros = False

        if self._matrix[pivot][col] == 0:
            all_zeros = True

            for i in range(pivot, self._matrix.__len__()):
                if self._matrix[i][col] != 0:
                    all_zeros = False
                    self.row_swap(pivot, i)
                    break
        
        return all_zeros
    

    '''
    Essentially an overloaded version of _sub_pivot_zero_check. This version of the method is used in calculating the matrix inverse, where row-operations are performed in parallel on two matrices: on the matrix whose inverse you're looking for and on the identity matrix.

    Params:
        pivot: index of pivot position
        col: index of column to check
        other: another matrix, if a row swap is performed on self, it's also performed on other

    Return: True if every entry including and below the pivot is zero, otherwise False
    '''
    def _sub_pivot_zero_check_other(self, pivot: int, col: int, other: 'Matrix'):
        all_zeros = False

        if self._matrix[pivot][col] == 0:
            all_zeros = True

            for i in range(pivot, self._matrix.__len__()):
                if self._matrix[i][col] != 0:
                    all_zeros = False
                    self.row_swap(pivot, i)
                    other.row_swap(pivot, i)
                    break
        
        return all_zeros


    '''
    For a column in self, it uses row combinations to make the entries in a row below the pivot entry into zeros.

    Params:
        pivot: index of pivot entry
        col: index of column to alter
    '''
    def _sub_pivot_eliminate(self, pivot: int, col: int):
        for i in range(pivot + 1, self._matrix.__len__()):
            if self._matrix[i][col] != 0:
                try:
                    scalar = -1 * (float(self._matrix[i][col]) / float(self._matrix[pivot][col]))
                except ZeroDivisionError: # necesary?
                    scalar = 0
                finally:
                    self.row_combination(pivot, i, scalar)


    '''
    Essentially an overloaded version of _sub_pivot_eliminate. This version of the method is used in calculating the matrix inverse, where row-operations are performed in parallel on two matrices: on the matrix whose inverse you're looking for and on the identity matrix.

    Params:
        pivot: index of pivot entry
        col: index of column to alter
        other: any row combinations performed on self are also performed on other
    '''
    def _sub_pivot_eliminate_other(self, pivot: int, col: int, other: 'Matrix'):
        for i in range(pivot + 1, self._matrix.__len__()):
            if self._matrix[i][col] != 0:
                try:
                    scalar = -1 * (float(self._matrix[i][col]) / float(self._matrix[pivot][col]))
                except ZeroDivisionError: # necesary?
                    scalar = 0
                finally:
                    self.row_combination(pivot, i, scalar)
                    other.row_combination(pivot, i, scalar)


    '''
    For a column in self, it uses row combinations to make the entries in a row above the pivot entry into zeros.

    Params:
        pivot: index of pivot entry
        col: index of column to alter
    '''
    def _sup_pivot_eliminate(self, pivot: int, col: int):
        for i in range(pivot - 1, -1, -1):
            if self._matrix[i][col] != 0:
                try:
                    scalar = -1 * (float(self._matrix[i][col]) / float(self._matrix[pivot][col]))
                except ZeroDivisionError: # necesary?
                    scalar = 0
                finally:
                    self.row_combination(pivot, i, scalar)
    
    
    '''
    Essentially an overloaded version of _sup_pivot_eliminate. This version of the method is used in calculating the matrix inverse, where row-operations are performed in parallel on two matrices: on the matrix whose inverse you're looking for and on the identity matrix.

    Params:
        pivot: index of pivot entry
        col: index of column to alter
        other: any row combinations performed on self are also performed on other
    '''
    def _sup_pivot_eliminate(self, pivot: int, col: int, other: 'Matrix'):
        for i in range(pivot - 1, -1, -1):
            if self._matrix[i][col] != 0:
                try:
                    scalar = -1 * (float(self._matrix[i][col]) / float(self._matrix[pivot][col]))
                except ZeroDivisionError: # necesary?
                    scalar = 0
                finally:
                    self.row_combination(pivot, i, scalar)
                    other.row_combination(pivot, i, scalar)
    
    
    '''
    Uses Gaussian Elimination to take self to row echelon form (REF)

    Params:
        decimals: decimal place to which to round REF of self before returning it; can input a negative number for no rounding.

    Return: Row echelon form of self
    '''
    def ref(self, decimals: int) -> 'Matrix':
        if self.initialized:
            ref = Matrix(self.rows, self.cols, True)
            ref.set_matrix_to_other(self)

            max_pivots = min(self.rows, self.cols) # max possible pivot positions
            pivot = -1 # initalize current pivot

            for i in range(max_pivots): # for every potential pivot column
                all_zeros = ref._sub_pivot_zero_check(pivot + 1, i)

                if not all_zeros:
                    pivot = pivot + 1
                    ref._sub_pivot_eliminate(pivot, i)

                    try:
                        scalar = 1 / float(ref._matrix[pivot][i])
                    except ZeroDivisionError:
                        scalar = 0
                    finally:
                        ref.row_scale(pivot, scalar)

            ref.round_matrix_n(decimals)
            return ref
        else:
            raise Exception('Cannot row-reduce uninitialized matrix')
    

    '''
    Essentially an overloaded version of ref. It also performs the row operations required to take self to row echelon form on another matrix

    Params:
        decimals: decimal place to which to round REF of self before returning it; can input a negative number for no rounding.
        other: while taking self to row echelon form, every row operation done to self is also done to other.

    Return: Tuple of row echelon form of self and other matrix
    '''
    def ref_other(self, decimals: int, other: 'Matrix'):
        if self.initialized:
            ref = Matrix(self.rows, self.cols, True)
            ref.set_matrix_to_other(self)

            if ref.rows == other.rows and ref.cols == other.cols:
                max_pivots = min(self.rows, self.cols) # max possible pivot positions
                pivot = -1 # initalize current pivot

                for i in range(max_pivots): # for every potential pivot column
                    all_zeros = ref._sub_pivot_zero_check_other(pivot + 1, i, other)

                    if not all_zeros:
                        pivot = pivot + 1
                        ref._sub_pivot_eliminate_other(pivot, i, other)

                        try:
                            scalar = 1 / float(ref._matrix[pivot][i])
                        except ZeroDivisionError:
                            scalar = 0
                        finally:
                            ref.row_scale(pivot, scalar)
                            other.row_scale(pivot, scalar)

                ref._round_matrix_n(decimals)
                other._round_matrix_n(decimals)
                return ref, other
            else:
                raise Exception('Other matrix dimensions not compatible with self')
        else:
            raise Exception('Cannot row-reduce uninitialized matrix')
    

    '''
    Takes dot product of self with another vector.

    Params:
        y: vector dotted with self. Self must be an n by 1 matrix and y must be a 1 by n matrix

    Return: dot product of self and y.
    '''
    def dot_prod(self, y):
        if not isinstance(y, Matrix):
            raise TypeError('Invalid argument type. Expected a Matrix')
        else:
            if y.cols != self.rows:
                raise Exception('Incompatible vector dimensions')
            elif y.rows > 1 or self.cols > 1:
                raise Exception('Invalid dimensions. Expected a vector')
            else:
                dot_prod = 0

                for i in range(self.rows):
                    dot_prod += float(self._matrix[i][0]) * float(y._matrix[0][i])

            return dot_prod
        

    '''
    Makes an n by 1 or 1 by n matrix out of a single row or single column of self

    Params:
        index: index of row or column to be made into matrix
        row: if True then index is treated as a row index, if False then index is treated as a column index

    Return: row or column as matrix
    '''
    def _row_col_to_vector(self, index: int, row: bool):
        if self.initialized:
            if row:
                row_vector = Matrix(1, self.cols, False)

                for i in range(self.cols):
                    row_vector._matrix.append([])
                    row_vector._matrix[0].append(self._matrix[index][i])

                return row_vector
            else:
                col_vector = Matrix(self.rows, 1, False)

                for i in range(self.rows):
                    col_vector._matrix.append([self._matrix[i][index]])

                return col_vector
        else:
            raise Exception('Cannot get row/column vectors from uninitialized matrix')
    
    
    '''
    Multiplies self by another Matrix.

    Params:
        b: Matrix to multiply self on the left, so b is applied to self
        decimals: decimals to place to which to round product

    Return: product of matrix multiplication
    '''
    def multiply(self, b, decimals: int) -> 'Matrix':
        if not isinstance(b, Matrix):
            raise TypeError('Invalid argument type. Expected a Matrix')
        else:
            if b.cols != self.rows:
                raise Exception('Incompatible matrix dimensions')
            else:
                product = Matrix(b.rows, self.cols, False)

                for i in range(product.rows):
                    product._matrix.append([])

                    for j in range(product.cols):
                        self_col_j = self._row_col_to_vector(j, False)
                        b_row_i = b._row_col_to_vector(i, True)

                        product._matrix[i].append(self_col_j.dot_prod(b_row_i))

                if decimals >= 0:
                    product.round_matrix_n(decimals)

                return product
            

    '''
    Method for generating identity matrices of any size

    Params:
        n: dimension of identity matrix

    Return: Identity matrix
    '''
    @staticmethod
    def identity(n: int):
        id = Matrix(n,  n, False)

        for i in range(n):
            id._matrix.append([])

            for j in range(n):
                if j == i:
                    id._matrix[i].append(1)
                else:
                    id._matrix[i].append(0)

        return id
    

    '''
    Uses Gaussian Elimination to take self to reduced row echelon form (RREF)

    Params:
        decimals: decimal place to which to round REF of self before returning it; can input a negative number for no rounding.

    Return: Reduced row echelon form of self
    '''
    def rref(self, decimals: int) -> 'Matrix':
        if self.initialized:
            rref = self.ref(decimals)
            max_pivots = min(rref.rows, rref.cols) # max possible pivot positions

            for i in range(max_pivots): # for each relevant column
                for j in range(rref.rows - 1, -1, -1): # for each entry bottom up
                    if rref._matrix[j][i] != 0:
                        rref._sup_pivot_eliminate(j, i)
                        break

            if decimals >= 0:
                rref.round_matrix_n(decimals)
            
            return rref
        else:
            raise Exception('Cannot row-reduce uninitialized matrix')
    

    '''
    Essentially an overloaded version of rref. It also performs the row operations required to take self to reduced row echelon form on another matrix

    Params:
        decimals: decimal place to which to round RREF of self before returning it; can input a negative number for no rounding.
        other: while taking self to reduced row echelon form, every row operation done to self is also done to other.

    Return: Tuple of reduced row echelon form of self and other matrix
    '''
    def rref_other(self, decimals: int, other: 'Matrix'):
        if self.initialized:
            rref, other = self.ref_other(decimals, other)
            max_pivots = min(rref.rows, rref.cols) # max possible pivot positions

            for i in range(max_pivots): # for each relevant column
                for j in range(rref.rows - 1, -1, -1): # for each entry bottom up
                    if rref._matrix[j][i] != 0:
                        rref._sup_pivot_eliminate_other(j, i, other) ## HERE -- Done
                        break

            if decimals >= 0:
                rref.round_matrix_n(decimals)
                other.round_matrix_n(decimals)
            
            return rref, other
        else:
            raise Exception('Cannot row-reduce uninitialized matrix')


    '''
    Takes self to upper triangular form, but doesn'y fully acheive row echelon form because main diagonal entries are not scaled so that they become 1.

    Params:
        decimals: decimal place to which to round utf

    Return: Upper triangular form
    '''
    def utf(self, decimals: int) -> 'Matrix':
        if self.initialized:
            utf = Matrix(self.rows, self.cols, True) # upper triangular form
            utf.set_matrix_to_other(self)

            max_pivots = min(self.rows, self.cols) # max possible pivot positions
            pivot = -1 # initalize current pivot

            for i in range(max_pivots): # for every potential pivot column
                all_zeros = utf._sub_pivot_zero_check(pivot + 1, i)

                if not all_zeros:
                    pivot = pivot + 1
                    utf._sub_pivot_eliminate(pivot, i)

            utf.round_matrix_n(decimals)
            return utf
        else:
            raise Exception('Cannot row-reduce uninitialized matrix')


    '''
    Takes the determinant of a square matrix

    Return: Determinant
    '''
    def det(self):
        if self.initialized:
            if self.rows == self.cols:
                utf = self.utf(-1) # we don't want rounding here
                
                det = 1

                for i in range(utf.cols):
                    det *= utf._matrix[i][i]

                return det
            else:
                raise Exception('Determinant requires square matrix')
        else:
            raise Exception('Cannot take determinant of uninitialized matrix')
        

    '''
    Determines whether or not a matrix is invertible

    Return: True if inivertible, False if not
    '''
    def is_inv(self):
        if self.initialized:
            if self.rows == self.cols:
                return self.det() != 0
            else:
                return False
        else:
            return False


    '''
    Uses Gaussian Elimination to transform identity matrix into inverse of self

    Params:
        decimals: decimal place to which to round inverse of self before returning; can input negative number for no rounding.

    Return: Matrix inverse of self
    '''
    def inv(self, decimals: int) -> 'Matrix':
        return self.rref_other(decimals, Matrix.identity(self.rows))[1]