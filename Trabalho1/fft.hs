import Data.List
import Data.Complex

odds [] = []
odds [x] = []
odds (e1:e2:xs) = e2 : odds xs

evens [] = []
evens [x] = []
evens (e1:e2:xs) = e1 : evens xs

f xs n = exp $ -(0:+2*pi)*n / genericLength xs

ffti [x,y] n = x + y * (exp $ -(0:+pi)*n)
ffti xs n = ffti (evens xs) n + f xs n * ffti (odds xs) n

fft xs [] = []
fft xs (y:ys) = (ffti xs y):(fft xs ys)

fft [0, 1, 4, 9] [0,1,2,3]