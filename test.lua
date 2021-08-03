local ma = require "matrixAlgebra"

local id1 = ma.liSparseMatrix.identity(3)
local id2 = ma.liSparseMatrix.identity(3)

local rand = ma.liSparseMatrix.random(6,6,-1,1)
local inverse = ma.liSparseMatrix.inverse(rand)
local id4 = ma.liSparseMatrix.sparsify(inverse * rand, 0.000000000001)

local id3 = ma.liSparseMatrix.unflatten(id1)

local twoid = id1 + id2
local zero = id1 - id2

