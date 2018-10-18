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
  , opposite
  , isNullElementary
  , ODM
  , Bounds
  ) where

  import Numeric.LinearAlgebra
  --import Data.Tree (Tree(..), unfoldTree)
  
  import ODMatrix (ODM)


  data ElementaryMatrix = 
    NullElementary Int |
    Elementary Int (Int, Int) (Int, Int)
    deriving (Show)

  -- | Size of the elementary matrix
  emSize :: ElementaryMatrix -> Int
  emSize (NullElementary s) = s
  emSize (Elementary s _ _) = s


  -- | Bounds of the capacity of the unit (Lower bound, Upper bound).
  type Bounds = (Double, Double) 

  
  -- Size, visited, s
  type ES = (ElementaryMatrix, ODM)


  applyElementaries :: ODM -> [ElementaryMatrix] -> [ODM]
  applyElementaries s es = map (applyElementary s) es
    

  -- | Apply a elementary operation to the given ODM.
  applyElementary :: ODM -> ElementaryMatrix -> ODM
  applyElementary a e = a + (toMatrix e)

  -- | Apply a sequence of elementary operations to the given ODM.
  applyPath :: ODM -> [ElementaryMatrix] -> ODM
  applyPath = foldl applyElementary


  -- | Given S, such that A + S = B, this function computes the first elementary permutation matrix E such that S' + E = S.
  getChildren :: ODM                  -- ^ Full path permutation matrix
               -> [ElementaryMatrix]  -- ^ All possible elementary matrices
  getChildren odm = es ++ map opposite es
    where al = toAssocList odm
          s = rows odm
          es = [ Elementary s (r1,c1) (r2,c2) | 
                  ((r1,c1),_) <- al, 
                  ((r2,c2),_) <- al, 
                  r1 > r2, c1 < c2 ]


  opposite :: ElementaryMatrix -> ElementaryMatrix
  opposite (Elementary s (r1,c1) (r2,c2)) = Elementary s (r2,c1) (r1,c2)
  opposite (NullElementary s) = NullElementary s


  -- getChildren s = foldl (\acc p -> acc ++ opposites s p pivots) [] pivots
  --   where -- Filtra solo los potenciales pivotes a los que puede agregarse o quitarse un pasajero.
  --         pivots = [ p | p@((i,j),x) <- toAssocList s, i <= j, x /= 0]  
        

  -- reciprocals :: Pivot -> [Pivot] -> [Pivot]
  -- reciprocals ((i,j),x) pivots = 
  --   [ p | p@((h,k),y) <- pivots, h < i, j < k, x/y > 0]


  -- opposites :: ODM -> Pivot -> [Pivot] -> [ElementaryMatrix]
  -- opposites s p@(i,j) pivots =
  --   [ Elementary (rows s) (i,j) (h,k)
  --     | ((h,k),_) <- reciprocals p pivots, 
  --       atIndex s (h,j) / x < 0,  -- Opposite signs
  --       atIndex s (i,k) / x < 0 ] 



  -- | Filter the elementary matrices that cannot be applied to the odm
  applicables :: Bounds             -- ^ Bounds of the values of the target
              -> ODM      -- ^ Target Origin-Destination Matrix
              -> [ElementaryMatrix] -- ^ List of potential elementary operations
              -> [ElementaryMatrix] -- ^ List of applicable elementary operations
  applicables b s = filter (applicable b s)

  -- | Indicates if an elementary operation is applicable for a given origin destination matrix with the given bounds.
  applicable :: Bounds            -- ^ Bounds of the values of the target
             -> ODM               -- ^ Target Origin-Destination Matrix
             -> ElementaryMatrix  -- ^ Elementary operation to apply
             -> Bool              
  applicable _ _ (NullElementary _) = True
  applicable (mn,mx) s (Elementary _ (i,j) (h,k)) = 
    p (i,k) > mn && p (h,k) < mx &&
    p (i,j) < mx && p (h,j) > mn
    where p = atIndex s



  -- | Compute the value of a elementary matrix in base of a distance functions between cells.
  elementaryValue :: (Int -> Int -> Double) -- ^ Measure function
                  -> ElementaryMatrix       -- ^ E                  
                  -> Double                 -- ^ Value of E
  elementaryValue _ (NullElementary _) = 0
  elementaryValue m (Elementary _ (r1,_ ) (r2,_)) = s * 2 * m r1 r2
    where s = fromIntegral . signum $ (r1 - r2)




  
  -- * Convertions

  -- | Transform the internal representation of a Elementary Matrix in a Matrix Double
  toMatrix :: ElementaryMatrix -> ODM
  toMatrix (NullElementary n) = assoc (n,n) 0 []
  toMatrix (Elementary n idx1@(r1,c1) idx2@(r2,c2)) = 
    assoc (n,n) 0 [(idx1,1), (idx2,1), ((r1,c2),-1), ((r2,c1),-1)]
                  

  -- | Transform a matrix to an assoc list.
  -- Recovers only the valid, different to 0, elements of the ODM.
  toAssocList :: ODM -> [((Int, Int), Double)]
  toAssocList s = [ c | c@((i,j),v) <- zip idxs elems, i <= j, v /= 0]
    where n = rows s
          idxs = [ (i,j) | i <- [0..(n-1)], j <- [0..(n-1)] ]
          elems = toList . flatten $ s


  -- | Return True if the elementary matrix is a null elementary matrix.
  isNullElementary :: ElementaryMatrix -> Bool
  isNullElementary (NullElementary s) = True
  isNullElementary _ = False