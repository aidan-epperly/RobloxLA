print(2 / 2)

local ma = require("matrixAlgebra")

local id = ma.sparseMatrix.identity(3)
local rand = ma.sparseMatrix.random(2,3,-1,1)
local randT = ma.sparseMatrix.transpose(rand)
local randSquare = ma.sparseMatrix.random(3,3,-1,1)
local lu = ma.sparseMatrix.lu(randSquare)
local ones = ma.sparseMatrix.new({{1,1},{1,1}})
local luOnes = ma.sparseMatrix.lu(ones * ones)
local cornerOnes = ma.sparseMatrix.new({{1,1,0},{1,1,0},{0,0,0}})
local luCornerOnes = ma.sparseMatrix.lu(cornerOnes)

print(luCornerOnes[1])

