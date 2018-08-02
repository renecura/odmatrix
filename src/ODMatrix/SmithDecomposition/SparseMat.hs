{-# LANGUAGE DeriveGeneric #-}

{-|
Module      : ODMatrix.SmithDecomposition.SparseMat
Description : Sparse representation of a matrix.

-}
module ODMatrix.SmithDecomposition.SparseMat (
      SparseMat(..)
    , denseSparseMat
    , sparseDenseMat
    , fixIndexes
    , makeAssoc
  ) where

    import GHC.Generics (Generic)
    --import Data.Aeson (ToJSON, FromJSON)
    
    import Numeric.LinearAlgebra as L hiding (rows, cols)
    import Numeric.LinearAlgebra.Data as LD (toList)

    
    -- | Sparse representation of a matrix.
    data SparseMat = SparseMat {
        rows :: Int
      , cols :: Int
      , cells :: AssocMatrix
      } deriving (Generic, Show)

    --instance ToJSON SparseMat
    --instance FromJSON SparseMat

    -- | Transform a sparse matrix in a regular one.
    denseSparseMat :: SparseMat -> Matrix Double
    denseSparseMat m = assoc (rows m, cols m) 0 (fixIndexes $ cells m)


    -- | Transform a regular matrix to an sparse one.
    sparseDenseMat :: Matrix Double-> SparseMat
    sparseDenseMat m = SparseMat r c cs
      where (r, c) = size m
            vals = toList . flatten $ m
            ids = [(i,j) | i <- [0..r-1], j <- [0..c-1]]
            cs = filter (\(_,v) -> v /= 0) $ zip ids vals


    
    -- * Support functions

    -- | Adjust the indexes of the association matrix from 1-starting to 0-starting.
    -- This is needed by the Linear Algebra library.
    fixIndexes :: (Num a) => [((a,a),b)] -> [((a,a),b)]
    fixIndexes = map (\((a,b),v) -> ((a-1,b-1),v))

    -- | Given a constant value and a list of indexes, construct a association matrix only containing that value.
    makeAssoc :: Double -> [(Int,Int)] -> AssocMatrix
    makeAssoc v = fixIndexes . map (flip (,) v)


    --------------------------------------------------------------------------------
    --------------------------------------------------------------------------------
    --------------------------------------------------------------------------------

    