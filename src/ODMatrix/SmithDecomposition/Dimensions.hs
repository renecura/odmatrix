{-|
Module      : Model.Decomposition.SparseMat
Description : Representación de matrices del sistema en forma dispersa.

-}
module ODMatrix.SmithDecomposition.Dimensions (
      srows
    , scols
    , nFromN
    , nFromM    
    ) where

    
    -- | Número de filas de sistema de ecuaciones.
    srows :: Integral a => a -> a
    srows n = 2 * n - 1

    -- | Número de columnas de sistema de ecuaciones.
    scols :: Integral a => a -> a
    scols n = n * (n + 1) `div` 2


    -- | Given N, computes the value of n, such that N = 2n + 1
    nFromN :: Int -> Int
    nFromN bigN = (bigN + 1) `div` 2

    -- | Given M, computes the value of n, such that M = n(n+1)/2
    nFromM :: Int -> Int
    nFromM bigM = ((floor $ sqrt (1 + 8 * m)) - 1) `div` 2
      where m = fromIntegral bigM :: Double

    