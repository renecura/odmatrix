{-|
Module      : ODMatrix
Description : Representación y procesamiento de matrices origen-destino.

Base de funciones para la operación sobre matrices origen-destino.


ODM Features
The ODM features are indicators to verify if two ODMs are elements of the same solution set.
-}
module ODMatrix (
    ODM
  , FeatVec
  --
  , boardingVector
  , alightingVector
  , onBoardVector
  --
  , inputVector
  --
  , getFeatures
  , compareFeatures
  --
  , toFlatList
  , fromFlatList
  ) where

  import Numeric.LinearAlgebra (Matrix, Vector, 
                                rows, cols, flatten, 
                                toColumns, toRows, 
                                subMatrix, vector,
                                vjoin, subVector)
  import Numeric.LinearAlgebra.Data (toList, takeDiag, fromLists, toLists, assoc)

  import Control.Monad (join)
  
  import ODMatrix.Util.TriangularNumbers (isTrg, trglvl)

  
  -- | Representación de una matriz de viajes. Debe ser cuadrada.
  type ODM = Matrix Double

  -- | The feature represents the input information.
  type FeatVec = Vector Double


  -- * Cálculo de vectores de información derivados. Feature Vectors.

  -- | Dada una matriz de viajes, calcula el vector de ascensos.
  boardingVector :: ODM      -- ^ Matriz de viajes de entrada.
                 -> FeatVec  -- ^ Vector de ascensos
  boardingVector = vector . map (sum . toList) . toRows

  -- | Dada una matriz de viajes, calcula el vector de descensos.
  alightingVector :: ODM      -- ^ Matriz de viajes de entrada.
                  -> FeatVec  -- ^ Vector de ascensos
  alightingVector = vector . map (sum . toList) . toColumns

  -- | Dada una matriz de viajes, calcula el vector de pasajeros a bordo.
  onBoardVector :: ODM      -- ^ Matriz de viajes de entrada.
                -> FeatVec  -- ^ Vector de pasajeros a bordo.
  onBoardVector m = vector ((sum . toList) . flatten <$> blocks)
    where s = cols m
          blocks = map (\x -> subMatrix (0,x) (x+1,s-x) m) [0..(s-1)]

  
  -- | Get all the feature vectors. (boarding, alighting, onBoard)
  getFeatures :: ODM -> (FeatVec, FeatVec, FeatVec)
  getFeatures m = (boardingVector m, alightingVector m, onBoardVector m)


  -- | Given two ODMs, verify if they are members of the same solution set.
  compareFeatures :: ODM -> ODM -> Bool
  compareFeatures a b = getFeatures a == getFeatures b



  -- | Given a ODM, build the vector of independent terms
  inputVector :: ODM -> Vector Double
  inputVector odm = vjoin [a,b]
    where a = alightingVector odm
          b = subVector 1 (cols odm - 1) $ boardingVector odm


  -- * Convertions

  -- | Aplana una matriz de origen destino, por filas y solo el triangulo superior.
  toFlatList :: ODM
             -> [Double]
  toFlatList = join . zipWith drop [0..] . toLists


  -- | Convierte una lista en una matriz origen destino.
  fromFlatList :: [Double]
               -> ODM
  fromFlatList ls
    | (not . isTrg) (length ls) =
      error "fromFlatList: Invalid list. The length must be a triangular number."
    | otherwise =
      assoc (n,n) 0 . zip [(i,j)|i<-[0..n-1],j<-[i..n-1]] $ ls
        where n = trglvl $ length ls
