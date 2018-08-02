{-|
Module      : SmithDecomposition
Description : Descomposición de matrices para la resoluciones de sistemas de ecuaciones diofánticas. Específico a la forma del sistema de ecuaciones asociado a la matriz de viajes.

Este módulo no provee funciones generales para encontrar la Forma Normal de Smith, sino que genera los resultados basado en el patron regular del sistema.
-}
module ODMatrix.SmithDecomposition (
    computeODM
  --
  , systemMatrix
  , snfR
  --
  , vectorToODM
  , odmToVector
  --
  , tVectorSize
  -- Publish from Dimensions
  , srows
  , scols
  , nFromN
  , nFromM
  --
  , generalSolution
  ) where

  import Prelude hiding ((<>))
  
  import Numeric.LinearAlgebra as L 

  import ODMatrix (ODM, FeatVec, toFlatList, fromFlatList)
  import Control.Monad (join)
  
  import ODMatrix.SmithDecomposition.SparseMat (SparseMat(SparseMat), denseSparseMat, makeAssoc)
  import ODMatrix.SmithDecomposition.Dimensions


  
  -- * Solution computation


  -- | Given a input vector and a t vector, computes a unique solution of the system.
  computeODM :: FeatVec      -- ^ Input Vector
             -> FeatVec      -- ^ t Vector
             -> ODM          -- ^ Solution Matrix
  computeODM b t = vectorToODM $ r1b + r2 <> asColumn t
    where (r1b, r2) = generalSolution b


  -- | Given and input vector, computes the general solution of the system.
  -- This is a partial solution. The aplicaction of a t vector is needed to find a particular one.
  generalSolution :: FeatVec      -- ^ Input Vector
                  -> (ODM, ODM)   -- ^ Solution Matrix
  generalSolution b = (r1 <> asColumn b, r2)
    where (r1, r2) = splitR $ snfR n
          n = nFromN $ size b


  -- * Smith's Decomposition
  
  -- | Given the size of the route, computes the system matrix.
  systemMatrix :: Int -> Matrix Double
  systemMatrix n = denseSparseMat $ SparseMat r c ab
    where r = srows n
          c = scols n
          a = zip (join $ map (\x -> [x..n]) [1..n]) [1..c]
          b = zip (join $ map (\x -> replicate (n-x) (n+x)) [1..n-1]) [n+1..c]
          ab = zip (a++b) $ repeat 1.0

  -- | Given the size of the route, computes the R matrix.
  snfR :: Int -> Matrix Double
  snfR n = stepA n <> stepB n <> stepC n

  
  stepA :: Int -> Matrix Double
  stepA n = accum (ident m) (+) as
    where as = makeAssoc (-1.0) $ zip sr [n+1..] 
          m = scols n
          sr = join [[k..n] | k <- [2..n]]

  
  stepB :: Int -> Matrix Double
  stepB n = accum (ident m) (+) $ as
    where m = scols n
          as = makeAssoc (-1.0) $ sbs
          sbs = join . map (\(p,k) -> zip [k,k..] [k+1..k+p]) . zip [n-2,n-3..1] $ ser n


  stepC :: Int -> Matrix Double
  stepC n = accum (ident m) (+) $ as ++ bs
    where m = scols n
          as = makeAssoc (-1.0) . join . map (\(a,b) -> [(a,a),(b,b)]) $ swp
          bs = makeAssoc ( 1.0) . join . map (\(a,b) -> [(a,b),(b,a)]) $ swp
          swp = tail . zip [n+1..] $ ser n

  
  ser :: Int -> [Int]
  ser n = foldl (\(a:ac) x -> (a-x):a:ac) [scols n] [2..n-1]


  -- * Computation of the t-vector

  -- | The size of the t vector given the size of the route.
  tVectorSize :: Int -> Int
  tVectorSize n = scols n - srows n



  -- * Support Functions


  -- | Break the right matrix R, for the computation of the solution
  splitR :: Matrix Double -> (Matrix Double, Matrix Double)
  splitR r = (takeColumns n r, dropColumns n r)
    where n = srows . nFromM $ cols r


  
  -- | Transform a vector in column matrix form to a ODM
  vectorToODM :: Matrix Double -> ODM
  vectorToODM v = fromFlatList . toList . flatten $ v  

  -- | Transform a ODM to vector in column matrix form
  odmToVector :: ODM -> Matrix Double
  odmToVector = matrix 1 . toFlatList




  ------------------------------------------------------------------------------------
  ------------------------------------------------------------------------------------
  ------------------------------------------------------------------------------------

  


  -- | Indices fila y columna de la matriz cuadrada triangular superior con la diagonal en cero, a partir del índice aplanado.
  --
  -- fromFlatIndex :: Int       -- ^ Tamaño de la matriz
  --               -> Int       -- ^ Índice plano
  --               -> (Int,Int) -- ^ (fila, columna)
  -- fromFlatIndex n k = as !! (k-1)
  --   where as = [1..n] >>= (\x -> [(x,j) | j <- [x+1..n]])


  -- | Indice aplanado de la matriz cuadrada triangular superior con la diagonal en cero.
  --
  -- toFlatIndex :: Int        -- ^ Tamaño de la matriz
  --             -> (Int, Int) -- ^ (fila, columna)
  --             -> Int        -- ^ Índice plano
  -- toFlatIndex' n (r,c)
  --   | r >= c = 0
  --   | otherwise = fr r - n + c
  --   where fr 0 = 0
  --         fr r = fr (r-1) + n - r

  -- toFlatIndex n (r,c)
  --   | r >= c = 0
  --   | otherwise = fr r - n + c
  --   where fr r = (r * n) - (r * (r + 1) `div` 2)

  -- * Construcción de matriz del sistema

  -- upperCond :: Int -> Int -> Int -> R
  -- upperCond n k c
  --   | (i <= k) && (k < j) = 1
  --   | otherwise           = 0
  --   where (i,j) = fromFlatIndex n c

  -- lowerCond :: Int -> Int -> Int -> R
  -- lowerCond n k c
  --   | (i == k+1) && (k+1 < j) = 1
  --   | otherwise               = 0
  --   where (i,j) = fromFlatIndex n c



  -- partialSystemMatrix :: (Int -> Int -> Int -> R) -> Int -> Matrix R
  -- partialSystemMatrix c n = (rs><cs) h
  --   where rs = n - 1
  --         cs = scols n
  --         l = [1..(scols n)]
  --         h = [ c n i j | i <- [1..rs], j <- l ]



  -- systemMatrix :: Int -> Matrix Double
  -- systemMatrix n = um === lm
  --   where um = partialSystemMatrix upperCond n
  --         lm = takeRows (n-2) $ partialSystemMatrix lowerCond n
  --
  --

  




  -- * System calculation

  

  -- | Fuck yeah!
  -- calculatorMachine :: Vector R             -- ^ Input vector
  --                   -> (Matrix R, Matrix R) -- ^ Two part calculator
  -- calculatorMachine b = (sub1, sub2)
  --   where n = nFromInputVector b
  --         (d,l,r) = snf n
  --         b' = asColumn b
  --         (r1,r2) = splitr n r
  --         sub1 = r1 <> l <> b'
  --         sub2 = r2

  
  -- calcODMatrix :: InputVector  -- ^ Input Vector
  --              -> TVector      -- ^ t Vector
  --              -> ODM          -- ^ Solution Vector
  -- calcODMatrix b t = s1 + s2 <> asColumn t
  --   where (s1, s2) = calculatorMachine b


  -- | Validate the solution of the system Ax = b and that all values of the solution are positive integers.
  -- solutionValidation :: Vector R  -- ^ Input Vector
  --                    -> Vector R  -- ^ Solution Vector
  --                    -> Bool      -- ^ True if validate.
  -- solutionValidation b x = p1 && p2
  --   where p1 = all (>= 0) $ toList x
  --         n = nFromInputVector b
  --         a = systemMatrix n
  --         p2 = a <> asColumn x == asColumn b






  -- tVectors :: Int         -- ^ N
  --          -> Int         -- ^ K
  --          -> [Vector R]  -- ^ All possible t vectors. Not necessarilly for a valid solution.
  -- tVectors n k = tVectors' n k v0
  --   where ss = scols n - srows n
  --         v0 = vector $ replicate ss 0


  -- tVectors' :: Int         -- ^ N
  --           -> Int         -- ^ K
  --           -> Vector R    -- ^ Initial vector for search.
  --           -> [Vector R]  -- ^ All possible t vectors. Not necessarilly for a valid solution.
  -- tVectors' n k v0 = vector <$> list
  --   where list = iterate (nextc k) $ toList v0
  --         nextc k ls =
  --           let (t:ts) = dropWhile (== fromIntegral k) ls
  --               s = length ls - length (t:ts)
  --               h = replicate s 0
  --           in h ++ [t+1] ++ ts

  -- tVector :: Vector R  -- ^ Input Vector
  --         -> Int       -- ^ K
  --         -> Vector R  -- ^ t Vector
  -- tVector b k = tVector' b k v0
  --   where n = srows' b
  --         v0 = vector $ replicate (tSize n) 0

  -- tVector' :: Vector R  -- ^ Input Vector
  --          -> Int       -- ^ K
  --          -> Vector R  -- ^ Initial vector for search
  --          -> Vector R  -- ^ t Vector
  -- tVector' b k v0 = head $ dropWhile cond $ tVectors' n k v0
  --   where n = nFromInputVector b
  --         cond t = validate k . solutionVectorToODMatrix . flatten $ calcODMatrix b t


  

  -- solutionVectorToODMatrix :: Vector R -> ODM
  -- solutionVectorToODMatrix = denseSparseMat . solutionVectorToSparseMat

  -- * System decomposition

  -- | Smith Normal Form
  -- Dado el tamaño N de una matriz de viajes A, retorna una tripla (D,L,R) donde D está en la Forma Normal de Smith y LAR = D.
  -- snf :: Int
  --     -> (Matrix R, Matrix R, Matrix R)
  -- snf n = (snfD n, snfL n, snfR n)
  -- --
  -- --
  -- snfD :: Int
  --      -> Matrix R
  -- snfD n = diagRect 0 (fromList $ replicate (srows n) 1) (srows n) (scols n)
  -- --
  -- --
  -- snfL :: Int
  --      -> Matrix R
  -- snfL = toDense . cells . sparseSnfL

  -- -- | Componente derecho de la descomposición de Smith.
  -- snfR :: Int
  --      -> Matrix R
  -- snfR = toDense . cells . sparseSnfR
