local matrixAlgebra = {}

local _sparseMatrixOps = {}

local _sparseMatrixFromTableOfTables = function(tableOfTables)
    local numberOfRows = #tableOfTables
    local numberOfColumns = #tableOfTables[1]

    local sparseForm = setmetatable({}, _sparseMatrixOps)

    for i = 1, numberOfRows do
        for j = 1, numberOfColumns do
            if tableOfTables[i][j] ~= 0 then
                sparseForm[numberOfColumns * (i - 1) + j] = tableOfTables[i][j]
            end
        end
    end

    rawset(sparseForm, "_dimensions", {numberOfRows, numberOfColumns})

    return sparseForm
end

local _sparseMatrixFromNumericTableOfTables = function(tableOfTables, tol)
    local numberOfRows = #tableOfTables
    local numberOfColumns = #tableOfTables[1]

    local sparseForm = setmetatable({}, _sparseMatrixOps)

    for i = 1, numberOfRows do
        for j = 1, numberOfColumns do
            if math.abs(tableOfTables[i][j]) >= tol then
                sparseForm[numberOfColumns * (i - 1) + j] = tableOfTables[i][j]
            end
        end
    end

    rawset(sparseForm, "_dimensions", {numberOfRows, numberOfColumns})

    return sparseForm
end 

local _sparseIdentity = function(n)
    local result = setmetatable({}, _sparseMatrixOps)

    rawset(result, "_dimensions", {n, n})

    for i = 1, n do
        result[n * (i - 1) + i] = 1
    end

    return result
end

local _sparseRandom = function(n, m, a, b, tol)
    local result = setmetatable({}, _sparseMatrixOps)

    rawset(result, "_dimensions", {n, m})

    for i = 1, n do
        for j = 1, m do
            local val = (b - a) * (math.random()) + a
            if math.abs(val) >= tol then
                result[n * (i - 1) + j] = val
            end
        end
    end

    return result
end

_sparseMatrixOps.__add = function(left, right)
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

_sparseMatrixOps.__mul = function(left, right)
    if left._dimensions[2] ~= right._dimensions[1] then
        error("Attempting to multiply incompatible sparse matrices.")
    end

    local result = setmetatable({}, _sparseMatrixOps)

    rawset(result, "_dimensions", {left._dimensions[1], right._dimensions[2]})

    for k, v in pairs(left) do
        if type(k) == "number" then
            local columnNumber = k % left._dimensions[2]
            if columnNumber == 0 then
                columnNumber = left._dimensions[2]
            end
            local rowNumber = (k - columnNumber) / left._dimensions[2]
            for i = 1, right._dimensions[2] do
                local val = result[right._dimensions[2] * rowNumber + i]
                if type(val) == "number" and right[right._dimensions[2] * (columnNumber - 1) + i] ~= nil then
                    val = val + v * right[right._dimensions[2] * (columnNumber - 1) + i]
                    result[right._dimensions[2] * rowNumber + i] = val
                elseif right[right._dimensions[2] * (columnNumber -   1) + i] ~= nil then
                    result[right._dimensions[2] * rowNumber + i] = v * right[right._dimensions[2] * (columnNumber - 1) + i]
                end
            end
        end
    end

    return result
end

_sparseMatrixOps.__tostring = function(matrix)
    local result = "{ "
    for k, v in pairs(matrix) do
        if type(k) == "number" then
            result = result .. tostring(k) .. ":" .. tostring(v) .. " "
        end
    end
    result = result .. "}"
    return result
end

matrixAlgebra.sparseMatrix = {}

matrixAlgebra.sparseMatrix.new = function(tableOfTables)
    return _sparseMatrixFromTableOfTables(tableOfTables)
end

matrixAlgebra.sparseMatrix.identity = function(n)
    return _sparseIdentity(n)
end

matrixAlgebra.sparseMatrix.random = function(n, m, a, b, tol)
    return _sparseRandom(n, m, a, b, tol)
end

return matrixAlgebra
