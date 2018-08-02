{-|
Module      : ODMatrix.ElementaryDecomposition
Description : This module is responsible for the representation of Elementary operations and the process of elementary decomposition of the path between Origin-Destination matrices.

-}
module ODMatrix.ElementaryDecomposition (
    ElementaryMatrix(..)
  , getChildren
  , applyElementary
  , applyPath
  , applicables
  , elementaryValue
  , ODM
  , Bounds
  ) where

  import Numeric.LinearAlgebra
  --import Data.Tree (Tree(..), unfoldTree)

  
  import ODMatrix (ODM)


  import Debug.Trace (trace)

  -- | Compact representation of an Elementary Matrix.
  -- The sign indicates de sign of the i,j and h,k element.
  data ElementaryMatrix = Elementary {
    size    :: Int,                  -- ^ Number of rows/columns
    indexes :: (Int, Int, Int, Int), -- ^ (i,j,h,k)
    sign    :: Double                -- ^ 1 o -1. 0 Indicates the null operation.
  } deriving (Show)


  -- | Bounds of the capacity of the unit (Lower bound, Upper bound).
  type Bounds = (Double, Double) 

  
  type Pivot = ((Int, Int), Double)


  nullElementary :: Int -> ElementaryMatrix
  nullElementary n = Elementary n (0,0,0,0) 0

  
  -- Size, visited, s
  type ES = (ElementaryMatrix, ODM)


  -- elementaryTree :: ODM         -- source
  --                -> ODM         -- target
  --                -> Tree ElementaryMatrix -- paths tree
  -- elementaryTree a b = fst <$> unfoldTree _eTree (nullElementary (rows a), b - a)



  -- _eTree :: ES -> (ES, [ES])
  -- _eTree (e,s) = ((e,s), zip children apps)
  --   where children = getChildren s
  --         apps = applyElementaries s children



  applyElementaries :: ODM -> [ElementaryMatrix] -> [ODM]
  applyElementaries s es
    | s == (toMatrix . nullElementary) (rows s) = []
    | otherwise = map (applyElementary s) es
    

  -- | Apply a elementary operation to the given ODM.
  applyElementary :: ODM -> ElementaryMatrix -> ODM
  applyElementary a e = a + (toMatrix e)

  -- | Apply a sequence of elementary operations to the given ODM.
  applyPath :: ODM -> [ElementaryMatrix] -> ODM
  applyPath = foldl applyElementary


  -- | Given S, such that A + S = B, this function computes the first elementary permutation matrix E such that S' + E = S.
  getChildren :: ODM       -- ^ Full path permutation matrix
               -> [ElementaryMatrix]  -- ^ All possible elementary matrices
  getChildren s = foldl (\acc p -> acc ++ opposites s p pivots) [] pivots
    where -- Filtra solo los potenciales pivotes a los que puede agregarse o quitarse un pasajero.
          pivots = [ p | p@((i,j),x) <- toAssocList s, i <= j, x /= 0]  
        

  reciprocals :: Pivot -> [Pivot] -> [Pivot]
  reciprocals ((i,j),x) pivots = 
    [ p | p@((h,k),y) <- pivots, h < i, j < k, x/y > 0]


  opposites :: ODM -> Pivot -> [Pivot] -> [ElementaryMatrix]
  opposites s p@((i,j),x) pivots =
    [ Elementary (rows s) (i,j,h,k) (signum x) 
      | ((h,k),_) <- reciprocals p pivots, 
        atIndex s (h,j) / x < 0,  -- Opposite signs
        atIndex s (i,k) / x < 0 ] 



  -- | Filter the elementary matrices that cannot be applied to the odm
  applicables :: Bounds             -- ^ Bounds of the values of the target
              -> ODM      -- ^ Target Origin-Destination Matrix
              -> [ElementaryMatrix] -- ^ List of potential elementary operations
              -> [ElementaryMatrix] -- ^ List of applicable elementary operations
  applicables b s = filter (applicable b s)

  -- | Indicates if an elementary operation is applicable for a given origin destination matrix with the given bounds.
  applicable :: Bounds            -- ^ Bounds of the values of the target
             -> ODM     -- ^ Target Origin-Destination Matrix
             -> ElementaryMatrix  -- ^ Elementary operation to apply
             -> Bool              
  applicable (mn,mx) s (Elementary _ (i,j,h,k) x)
      | x ==   1  = p (i,j) < mx && p (h,k) < mx && p (i,k) > mn && p (h,j) > mn
      | x == (-1) = p (i,j) > mn && p (h,k) > mn && p (i,k) < mx && p (h,j) < mx
      | otherwise = error "ElementaryMatrix can only contain 1 o -1"
    where p = atIndex s



  -- | Compute the value of a elementary matrix in base of a distance functions between cells.
  elementaryValue :: (Int -> Int -> Double) -- ^ Measure function
                  -> ElementaryMatrix       -- ^ E                  
                  -> Double                 -- ^ Value of E
  elementaryValue m (Elementary _ (i,_,h,_) s) = s * 2 * m i h




  
  -- * Convertions

  -- | Transform the internal representation of a Elementary Matrix in a Matrix Double
  toMatrix :: ElementaryMatrix -> ODM
  toMatrix (Elementary n (i,j,h,k) s) = 
    assoc (n,n) 0 [((i,j),s), ((h,k),s), ((i,k),-s), ((h,j),-s)]
                  

  -- | Transform a matrix to an assoc list
  toAssocList :: ODM -> [((Int, Int), Double)]
  toAssocList s = zip [ (i,j) | i<-[0..(n-1)], j<-[0..(m-1)]] (toList . flatten $ s)
    where n = rows s
          m = cols s


  