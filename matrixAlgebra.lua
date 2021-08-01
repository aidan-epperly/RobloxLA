local matrixAlgebra = {}

local _sparseMatrix = {}

local _sparseMatrixFromTableOfTables = function (tableOfTables)
    local numberOfRows = #tableOfTables
    local numberOfColumns = #tableOfTables[1]

    local sparseForm = setmetatable({}, _sparseMatrix)

    for i = 1, numberOfRows do
        local rowConstant = numberOfColumns * (i - 1)
        for j = 1, numberOfColumns do
            if tableOfTables[i][j] ~= 0 then
                sparseForm[rowConstant + j] = tableOfTables[i][j]
            end
        end
    end

    rawset(sparseForm, "_dimensions", {numberOfRows, numberOfColumns})

    return sparseForm
end

local _sparseMatrixFromNumericTableOfTables = function (tableOfTables, tol)
    tol = tol or 0

    local numberOfRows = #tableOfTables
    local numberOfColumns = #tableOfTables[1]

    local sparseForm = setmetatable({}, _sparseMatrix)

    for i = 1, numberOfRows do
        local rowConstant = numberOfColumns * (i - 1)
        for j = 1, numberOfColumns do
            if math.abs(tableOfTables[i][j]) >= tol then
                sparseForm[rowConstant + j] = tableOfTables[i][j]
            end
        end
    end

    rawset(sparseForm, "_dimensions", {numberOfRows, numberOfColumns})

    return sparseForm
end

local _sparseCopy = function (matrix)
    local copy = setmetatable({}, _sparseMatrix)

    rawset(copy, "_dimensions", {matrix._dimensions[1], matrix._dimensions[2]})

    for k, v in pairs(matrix) do
        if type(k) == "number" then
            copy[k] = v
        end
    end

    return copy
end

local _sparseIdentity = function (n)
    local result = setmetatable({}, _sparseMatrix)

    rawset(result, "_dimensions", {n, n})

    for i = 1, n do
        result[n * (i - 1) + i] = 1
    end

    return result
end

local _sparseZero = function (n, m)
    local result = setmetatable({}, _sparseMatrix)

    rawset(result, "_dimensions", {n,m})

    return result
end

local _sparseRandom = function (n, m, a, b, tol)
    tol = tol or 0

    local result = setmetatable({}, _sparseMatrix)

    rawset(result, "_dimensions", {n, m})

    for i = 1, n do
        local rowConstant = m * (i - 1)
        for j = 1, m do
            local val = (b - a) * (math.random()) + a
            if math.abs(val) >= tol then
                result[rowConstant + j] = val
            end
        end
    end

    return result
end

local _sparseColumnSwap = function (matrix, n, m)
    local numberOfRows = matrix._dimensions[1]
    local numberOfColumns = matrix._dimensions[2]

    for i = 1, numberOfRows do
        local rowConstant = numberOfColumns * (i - 1)
        matrix[rowConstant + n],  matrix[rowConstant + m] =  matrix[rowConstant + m],  matrix[rowConstant + n]
    end

    return matrix
end

local _sparseRightPermutationMatrix = function (n, permutations)
    local result = _sparseIdentity(n)

    for i = 1, n do
        local destination = permutations[i]
        if destination ~= nil then
            result = _sparseColumnSwap(result, i, destination)
        end
    end

    return result
end

local _sparseRowSwap = function (matrix, n, m)
    local numberOfRows = matrix._dimensions[1]
    local numberOfColumns = matrix._dimensions[2]

    for i = 1, numberOfColumns do
        local rowConstant = numberOfColumns * (n - 1)
        local rowConstant2 = numberOfColumns * (m - 1)
        matrix[rowConstant + i],  matrix[rowConstant2 + i] =  matrix[rowConstant2 + i],  matrix[rowConstant + i]
    end

    return matrix
end

local _sparsePermutationMatrix = function (n, permutations)
    local result = _sparseIdentity(n)

    for i = 1, n do
        local destination = permutations[i]
        if destination ~= nil then
            result = _sparseRowSwap(result, i, destination)
        end
    end

    return result
end

local _sparseShear = function (matrix, n, m, c, tol)
    tol = tol or 0

    local numberOfColumns = matrix._dimensions[2]

    for i = 1, numberOfColumns do
        local val1 = matrix[numberOfColumns * (n - 1) + i]
        local val2 = matrix[numberOfColumns * (m - 1) + i]
        if val1 ~= nil and val2 ~= nil then
            local val3 = val2 + c * val1
            if math.abs(val3) < tol then
                matrix[numberOfColumns * (m - 1) + i] = nil
            else
                matrix[numberOfColumns * (m - 1) + i] = val3
            end
        elseif val1 ~= nil then
            local val3 = c * val1
            if math.abs(val3) < tol then
                matrix[numberOfColumns * (m - 1) + i] = nil
            else
                matrix[numberOfColumns * (m - 1) + i] = val3
            end
        end
    end

    return matrix
end

local _sparseRowAdd = function (matrix, row, n, c, tol)
    tol = tol or 0

    local numberOfColumns = matrix._dimensions[2]

    for i = 1, numberOfColumns do
        local val1 = matrix[numberOfColumns * (n - 1) + i]
        local val2 = row[i]
        if val1 ~= nil and val2 ~= nil then
            local val3 = val1 + c * val2
            if math.abs(val3) < tol then
                matrix[numberOfColumns * (n - 1) + i] = nil
            else
                matrix[numberOfColumns * (n - 1) + i] = val3
            end
        elseif val2 ~= nil then
            local val3 = c * val2
            if math.abs(val3) < tol then
                matrix[numberOfColumns * (n - 1) + i] = nil
            else
                matrix[numberOfColumns * (n - 1) + i] = val3
            end
        end
    end

    return matrix
end

local _sparseTranspose = function (matrix)
    local numberOfRows = matrix._dimensions[1]
    local numberOfColumns = matrix._dimensions[2]

    local result = setmetatable({}, _sparseMatrix)

    rawset(result, "_dimensions", {numberOfColumns, numberOfRows})

    for k, v in pairs(matrix) do
        if type(k) == "number" then
            local columnNumber = k % numberOfColumns
            if columnNumber == 0 then
                columnNumber = numberOfColumns
            end
            local rowNumber = (k - columnNumber) / numberOfColumns + 1
            result[numberOfRows * (columnNumber - 1) + rowNumber] = v
        end
    end

    return result
end

local _sparseColumnMax = function (matrix, n, m)
    local numberOfColumns = matrix._dimensions[2]
    local numberOfRows = matrix._dimensions[1]

    local max = 0
    local maxRow = 1

    if matrix._dimensions[1] > #matrix - m - 1 then
        for k, v in pairs(matrix) do
            if type(k) == "number" and k > m * numberOfColumns and (k - n) % numberOfColumns == 0 and v > max then
                max = v
                maxRow = math.floor((k - n) / numberOfColumns + 0.1)
            end
        end
    else
        for i = m, numberOfRows do
            
        end
    end
end

local _sparseLU = function (matrix)
    local numberOfRows = matrix._dimensions[1]
    local numberOfColumns = matrix._dimensions[2]

    if numberOfRows ~= numberOfColumns then
        error("Cannot compute LU of sparse rectangular matrix.")
    end

    local permuations = {}

    local l = _sparseIdentity(numberOfRows)

    for i = 1, numberOfColumns - 1 do
        local maxRow = i
        local max = matrix[numberOfColumns * (i - 1) + i] or 0

        for j = i, numberOfRows do
            local maxCandidate = matrix[numberOfColumns * (j - 1) + i] or 0
            maxCandidate = math.abs(maxCandidate)
            if maxCandidate ~= nil and maxCandidate > max then
                max = maxCandidate
                maxRow = j
            end
        end

        print(matrix)

        if max == 0 then
            error("Sparse matrix is not invertible")
        end

        if maxRow ~= i then
            _sparseRowSwap(matrix, i, maxRow)
            permuations[i] = maxRow
        end

        max = matrix[numberOfColumns * (i - 1) + i]

        for j = i + 1, numberOfRows do
            local val = matrix[numberOfColumns * (j - 1) + i]
            local valOverMax = val / max
            if val ~= nil then
                l[numberOfColumns * (j - 1) + i] = valOverMax
            end
            matrix[numberOfColumns * (j - 1) + i] = nil
            for k = i + 1, numberOfColumns do
                local val1 = matrix[numberOfColumns * (j - 1) + k]
                local val2 = matrix[numberOfColumns * (i - 1) + k]
                if val1 ~= nil and val2 ~= nil then
                    matrix[numberOfColumns * (j - 1) + k] = val1 - val2 * valOverMax
                elseif val2 ~= nil then
                    matrix[numberOfColumns * (j - 1) + k] = -val2 * valOverMax
                end
            end
        end
    end

    local permuationMatrix = _sparsePermutationMatrix(numberOfRows, permuations)

    return {l, matrix, permuationMatrix}
end

local _sparseInverse = function (matrix)
    local numberOfRows = matrix._dimensions[1]
    local numberOfColumns = matrix._dimensions[2]

    if numberOfRows ~= numberOfColumns then
        error("Cannot compute inverse of sparse rectangular matrix.")
    end

    local permuations = {}

    local result = _sparseIdentity(numberOfRows)

    for i = 1, numberOfColumns - 1 do
        local maxRow = i
        local max = matrix[numberOfColumns * (i - 1) + i] or 0

        for j = i, numberOfRows do
            local maxCandidate = math.abs(matrix[numberOfColumns * (j - 1) + i])
            if maxCandidate ~= nil and maxCandidate > max then
                max = maxCandidate
                maxRow = j
            end
        end

        if max == 0 then
            error("Sparse matrix is not invertible")
        end

        if maxRow ~= i then
            _sparseRowSwap(matrix, i, maxRow)
            permuations[i] = maxRow
        end

        max = matrix[numberOfColumns * (i - 1) + i]

        for j = i + 1, numberOfRows do
            local val = matrix[numberOfColumns * (j - 1) + i]
            local valOverMax = val / max
            matrix[numberOfColumns * (j - 1) + i] = nil
            for k = i + 1, numberOfColumns do
                local val1 = matrix[numberOfColumns * (j - 1) + k]
                local val2 = matrix[numberOfColumns * (i - 1) + k]
                if val1 ~= nil and val2 ~= nil then
                    matrix[numberOfColumns * (j - 1) + k] = val1 - val2 * valOverMax
                elseif val2 ~= nil then
                    matrix[numberOfColumns * (j - 1) + k] = -val2 * valOverMax
                end
                val1 = result[numberOfColumns * (j - 1) + k]
                val2 = result[numberOfColumns * (i - 1) + k]
                if val1 ~= nil and val2 ~= nil then
                    result[numberOfColumns * (j - 1) + k] = val1 - val2 * valOverMax
                elseif val2 ~= nil then
                    result[numberOfColumns * (j - 1) + k] = -val2 * valOverMax
                end
            end
        end
    end

    for i = numberOfRows, 1, -1 do
        local rowConstant = numberOfColumns * (i - 1)
        local val = matrix[rowConstant + i]
        for j = 1, numberOfRows - 1 do
            local val1 = matrix[numberOfColumns * (i - j - 1)]
            if val1 ~= nil then
                _sparseShear(result, i, i - j, val1 / val)
            end
        end
    end

    local permuationMatrix = _sparsePermutationMatrix(numberOfRows, permuations)

    return result
end

_sparseMatrix.__add = function (left, right)
    if left._dimensions[1] ~= right._dimensions[1] or left._dimensions[2] ~= right._dimensions[2] then
        error("Attempting to add sparse matrices of different sizes.")
    end
    for k, v in pairs(left) do
        if type(k) == "number" then
            local val = right[k]
            if type(val) == "number" then
                left[k] = v + val
            end
        end
    end
end

_sparseMatrix.__sub = function (left, right)
    if left._dimensions[1] ~= right._dimensions[1] or left._dimensions[2] ~= right._dimensions[2] then
        error("Attempting to add sparse matrices of different sizes.")
    end
    for k, v in pairs(left) do
        if type(k) == "number" then
            local val = right[k]
            if type(val) == "number" then
                left[k] = v - val
            end
        end
    end
end

_sparseMatrix.__mul = function (left, right)
    if left._dimensions[2] ~= right._dimensions[1] then
        error("Attempting to multiply incompatible sparse matrices.")
    end

    local result = setmetatable({}, _sparseMatrix)

    rawset(result, "_dimensions", {left._dimensions[1], right._dimensions[2]})

    for k, v in pairs(left) do
        if type(k) == "number" then
            local columnNumber = k % left._dimensions[2]
            if columnNumber == 0 then
                columnNumber = left._dimensions[2]
            end
            local rowNumber = (k - columnNumber) / left._dimensions[2]
            local rowConstant = right._dimensions[2] * (columnNumber - 1)
            local rowConstant2 = right._dimensions[2] * rowNumber
            for i = 1, right._dimensions[2] do
                local val1 = result[rowConstant2 + i]
                local val2 = right[rowConstant + i]
                if type(val1) == "number" and val2 ~= nil then
                    val1 = val1 + v * right[rowConstant + i]
                    result[rowConstant2 + i] = val1
                elseif val2 ~= nil then
                    result[rowConstant2 + i] = v * right[rowConstant + i]
                end
            end
        end
    end

    return result
end

_sparseMatrix.__tostring = function (matrix)
    local result = "{"

    local numberOfRows = matrix._dimensions[1]
    local numberOfColumns = matrix._dimensions[2]

    for i = 1, numberOfRows - 1 do
        result = result .. "{"
        for j = 1, numberOfColumns - 1 do
            local val = matrix[numberOfColumns * (i - 1) + j] or 0
            result = result .. tostring(val) .. ","
        end
        local val = matrix[numberOfColumns * (i - 1) + numberOfColumns] or 0
        result = result .. tostring(val) .. "},"
    end

    result = result .. "{"
    for j = 1, numberOfColumns - 1 do
        local val = matrix[numberOfColumns * (numberOfRows - 1) + j] or 0
        result = result .. tostring(val) .. ","
    end
    local val = matrix[numberOfColumns * (numberOfRows - 1) + numberOfColumns] or 0
    result = result .. tostring(val) .. "}"

    result = result .. "}"

    return result
end

local _matrix = {}

local _matrixFromTableOfTables = function (tableOfTables)
    local result = setmetatable(tableOfTables, _matrix)

    rawset(result, "_dimensions", {#tableOfTables,#tableOfTables[1]})

    return result
end

local _nthStandardBasisVector = function (n, m)
    local e = {}

    for i = 1, m do
        if i == n then
            e[i] = 1
        else
            e[i] = 0
        end
    end

    return e
end

local _zeroVector = function (n)
    local z = {}

    for i = 1, n do
        z[i] = 0
    end

    return z
end

local _matrixIdentity = function (n)
    local result = {}

    for i = 1, n do
        result[i] = _nthStandardBasisVector(i, n)
    end

    return _matrixFromTableOfTables(result)
end

local _matrixZero = function (n, m)
    local result = {}

    for i = 1, n do
        result[i] = _zeroVector(m)
    end

    return _matrixFromTableOfTables(result)
end

_matrix.__add = function (left, right)
    if left._dimensions[1] ~= right._dimensions[1] or left._dimensions[2] ~= right._dimensions[2] then
        error("Attempting to add matrices of different sizes.")
    end

    local result = {}

    for k, v in ipairs(left) do
        result[k] = {}
        for kk, vv in ipairs(v) do
            result[k] = vv + right[k][kk]
        end
    end

    return _matrixFromTableOfTables(result)
end

_matrix.__sub = function (left, right)
    if left._dimensions[1] ~= right._dimensions[1] or left._dimensions[2] ~= right._dimensions[2] then
        error("Attempting to add matrices of different sizes.")
    end

    local result = {}

    for k, v in ipairs(left) do
        result[k] = {}
        for kk, vv in ipairs(v) do
            result[k] = vv - right[k][kk]
        end
    end

    return _matrixFromTableOfTables(result)
end

matrixAlgebra.sparseMatrix = {}

matrixAlgebra.sparseMatrix.new = function (tableOfTables)
    return _sparseMatrixFromTableOfTables(tableOfTables)
end

matrixAlgebra.sparseMatrix.newNumeric = function (tableOfTables, tol)
    return _sparseMatrixFromNumericTableOfTables(tableOfTables, tol)
end

matrixAlgebra.sparseMatrix.identity = function (n)
    return _sparseIdentity(n)
end

matrixAlgebra.sparseMatrix.zero = function (n, m)
    return _sparseZero(n, m)
end

matrixAlgebra.sparseMatrix.random = function (n, m, a, b, tol)
    return _sparseRandom(n, m, a, b, tol)
end

matrixAlgebra.sparseMatrix.lu = function (input)
    local matrix = _sparseCopy(input)
    return _sparseLU(matrix)
end

matrixAlgebra.sparseMatrix.transpose = function (input)
    local matrix = _sparseCopy(input)
    return _sparseTranspose(matrix)
end

matrixAlgebra.matrix = {}

matrixAlgebra.matrix.new = function (tableOfTables)
    return _matrixFromTableOfTables(tableOfTables)
end

matrixAlgebra.matrix.identity = function (n)
    return _matrixIdentity(n)
end

matrixAlgebra.matrix.zero = function (n, m)
    return _matrixZero(n, m)
end

return matrixAlgebra
