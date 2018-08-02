{-|
Module      : Util.TriangularNumbers
Description : Cálculos asociados a los números triangulares.
-}
module ODMatrix.Util.TriangularNumbers where

  -- | The nth triangular number.
  trgn :: Integral a => a -> a
  trgn n = (n * (n + 1)) `div` 2

  -- | Triangular root of x.
  trgrt :: (Integral a, Floating b) => a -> b
  trgrt x = (sqrt (fromIntegral x * 8 + 1) - 1) / 2

  -- | Triangular level. Is the integer part of trgrt.
  trglvl :: Integral a => a -> a
  trglvl = floor . trgrt

  -- | Verifica si es un numero triangular.
  isTrg :: Integral a => a -> Bool
  isTrg x = trgrt x == (fromIntegral . floor) (trgrt x)
