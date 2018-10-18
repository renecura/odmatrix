{-|
Module      : ODMatrix.ElementaryDecomposition.ShortestPath
Description : This module is responsible for computing the shortest path between two ODMs.

Contains the necesary functions to measure and compare ODMs.

-}
module ODMatrix.ElementaryDecomposition.ShortestPath (
      Measure
    , sePath
    , ev
    , ev3
  ) where

  import Control.Monad.State.Lazy 
  --import Control.Monad.Loops (iterateUntilM)
  
  import Numeric.LinearAlgebra

  import qualified Data.HashMap.Lazy as M
  import qualified Data.Hashable as H
  import Data.Maybe

  

  import ODMatrix.ElementaryDecomposition (ElementaryMatrix(..), ODM, 
                                       getChildren, applyElementary, 
                                       applyPath,
                                       applicables, elementaryValue, 
                                       opposite, isNullElementary)


  -- import Debug.Trace (trace)
  

  -- | Measure or Metric prototype
  type Measure = Int -> Int -> Double


  -- | Elementary Value.
  -- Given a metric m, a in route metric mr and a scale factor alpha returns a elementary value metric.
  ev :: Double    -- ^ Scale factor alpha
     -> Measure   -- ^ Metric m
     -> Measure   -- ^ In route metric mr
     -> Measure   -- ^ Elementary value metric
  ev a m mr i j
    | a < 0 = error "The ajustment parameter must be greater than 0."
    | otherwise = m i j + a * mr i j



  -- | Shortest Path
  -- Computes the sequence of elementary matrices between ODMs
  sePath :: Measure             -- ^ Elementary value metric
         -> ODM                 -- ^ Source ODM
         -> ODM                 -- ^ Target ODM
         -> [ElementaryMatrix]  -- ^ Sequence of elementary matrices in the shortest path
  sePath m a b = evalState _process (_initialSEPState m a b)



  -- | Map a Elementary Matrix to it related cells indexes
  emc :: ElementaryMatrix -> (Int, Int)
  emc (Elementary _ (cell1, _) (cell2, _)) = (cell1,cell2)



  -- | Compares two ODMs
  odmc :: (Matrix Double -> Double) -- ^ Function g to compute the remainder
       -> Measure                   -- ^ Elementary value metric
       -> Double                    -- ^ Scale factor
       -> ODM                       -- ^ ODM A
       -> ODM                       -- ^ ODM B
       -> Double                    -- ^ Difference between A and B
  odmc g ev beta a b = beta * g s + sum [ev i j | (i,j) <- sp]
    where p = sePath ev a b
          sp = map emc . filter (not . isNullElementary) $ p
          s = b - applyPath a p -- Remainder




  newtype HashMatrix = HashMatrix { hm :: ODM }

  instance Show HashMatrix where
    show (HashMatrix m) = show m

  instance Eq HashMatrix where
    (==) (HashMatrix a) (HashMatrix b) = a == b
  
  instance H.Hashable HashMatrix where
    hashWithSalt salt (HashMatrix m) = floor $ m `atIndex` (salt `mod` rows m, salt `mod` cols m)


  
  







  data SEP = SEP {
      dist      :: M.HashMap HashMatrix Double
    , prev      :: M.HashMap HashMatrix ElementaryMatrix
    , visited   :: [ODM]
    , target    :: ODM
    , measure   :: Measure
    } 

  instance Show SEP where
    show = show . visited 

  
  


  
  


  -- | Elementary Value 3 (See the paper).
  -- This is a special case of elementary value computed with only one metric.
  -- Only works for regular grids.
  ev3 :: Double -> Measure -> Measure
  ev3 a m i j
    | a <= 0 || a >= 1 = error "The ajustment parameter a must be 0 < a < 1."
    | otherwise = a * m i j + (1-a) * sum [ m k (k+1) | k <- [min i j .. max i j] ]



  -- TODO: Implement metrics.

  _initialSEPState :: Measure -> ODM -> ODM -> SEP
  _initialSEPState m src trg = SEP (M.singleton (HashMatrix src) 0) M.empty [] trg m
 

  -- | TODO: take out the harcoded 60.
  -- | MAYBE IF IN THE FINAL STEP ASK FOR NEIGHBORS EMPTY.
  _process :: State SEP [ElementaryMatrix]  
  _process = do
    s <- get
    let (curr,d) = _getNext (dist s)
        neighbors = applicables (0,60) curr . getChildren $ (target s) - curr
        -- Filter out the visited neighbors
        unvisited = filter (\e -> applyElementary curr e `notElem` visited s) neighbors
    -- Calculate min distance to all neighbors (unvisited)
    modify $ minAndPrev curr d unvisited 
    -- \t -> t {dist = foldr (mindist (measure t) curr d) (dist t) unvisited}
    -- Add current to the previous of unvisited 
    
    -- Mark the current as visited.
    modify $ \t -> t {dist = M.delete (HashMatrix curr) (dist t), 
                      visited = curr:visited t}
    sn <- get
    if curr == (target sn) || null neighbors  then return . fst $ _wrapUp curr sn else _process
    -- if curr == (target sn) then return . fst $ _wrapUp curr sn else _process
  

  -- Calculate the min distance and the corresponding previous.
  minAndPrev :: ODM -> Double -> [ElementaryMatrix] -> SEP -> SEP
  minAndPrev curr d unvisited t = foldr f t nds
    where nds = map (neighborsDists t curr d) unvisited
          f (md,k,e) s = s { dist = M.insert k md (dist s),
                             prev = M.insert k e (prev s) }
                  

  
  neighborsDists :: SEP                           -- ^ Final state
                 -> ODM                           -- ^ Current
                 -> Double                        -- ^ Dist from origin
                 -> ElementaryMatrix              -- ^ Neighbor
                 -> (Double,HashMatrix,ElementaryMatrix)     -- ^ Min distance and Neighbor
  neighborsDists s m d e = (min newd oldd, k, e)
    where newd = d + elementaryValue (measure s) e
          oldd = M.lookupDefault newd k (dist s)
          k = HashMatrix $ applyElementary m e
  

  
   

  -- | Get next matrix to process and the shortest distance from src
  _getNext :: M.HashMap HashMatrix Double -- ^ List of distances by ODM
           -> (ODM,Double)                -- ^ ODM with min distance and the distance
  _getNext ds 
    | null ds = error "Not destination found"
    | otherwise = (\(m,d) -> (hm m, d)) $ M.foldrWithKey f initial ds
    where initial = head . M.toList $ ds            -- (k,v)
          f k v mn@(_, minv) | v < minv = (k,v)
                             | otherwise = mn

  
  _wrapUp :: ODM -> SEP -> ([ElementaryMatrix],SEP)
  _wrapUp m s = (es, s)
    where p = M.lookup (HashMatrix m) (prev s)
          --op e = e {sign = - sign e}
          odm = (\x -> applyElementary m (opposite x)) <$> p
          es = maybe [] (\x -> fromJust p : fst (_wrapUp x s)) odm
  


  -- _testMatrix :: ODM
  -- _testMatrix = matrix 10 [0,0,0,1,0,0,0,0,0,0,0,2,1,0,1,3,0,0,0,0,0,0,2,0,0,4,1,0,1,0,0,0,0,1,0,0,0,0,2,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

  -- _testSimMatrix :: ODM
  -- _testSimMatrix = matrix 10 [0,0,1,0,0,0,0,0,0,0,0,2,0,1,1,3,0,0,0,0,0,0,2,0,0,4,1,0,1,0,0,0,0,1,0,0,0,0,2,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

  -- _testDifMatrix :: ODM
  -- _testDifMatrix = matrix 10 [0,0,0,1,0,0,0,0,0,0,0,2,1,0,1,3,0,0,0,0,0,0,2,0,0,4,2,0,0,0,0,0,0,1,0,0,0,0,2,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]