-- mightymandel -- GPU-based Mandelbrot Set explorer
-- Copyright (C) 2012,2013,2014,2015 Claude Heiland-Allen
-- License GPL3+ http://www.gnu.org/licenses/gpl.html

{-# LANGUAGE FlexibleInstances #-}
module Main (main) where

import Control.Monad (forM_, unless)
import Control.Monad.Trans.RWS (put, get, tell, execRWS, RWS)
import Data.Either (lefts, rights, isRight)
import qualified Data.Foldable as F
import Data.Function (on)
import Data.List (group, groupBy, sort, sortBy, intercalate, partition)
import Data.Map (Map)
import qualified Data.Map as M
import Data.Monoid (Monoid, (<>))
import Data.Set (Set)
import qualified Data.Set as S

-- new in ghc-7.10
sortOn :: Ord o => (a -> o) -> [a] -> [a]
sortOn f = sortBy (compare `on` f)

groupOn :: Eq e => (a -> e) -> [a] -> [[a]]
groupOn f = groupBy ((==) `on` f)

class Pretty e where
  pretty :: e -> String

class Simplify e where
  simplify :: e -> e

class Collect e where
  collect :: e -> e

normalize :: (Simplify e, Collect e) => e -> e
normalize = simplify . collect . simplify

data E = Sum [E] | Product [E] | N Integer | A Integer | Z | C | DZ | DC
  deriving (Read, Show, Eq, Ord)

instance Pretty E where
  pretty DC = "\\D{c}}"
  pretty DZ = "\\D{z}}"
  pretty C = "c"
  pretty Z = "z"
  pretty (A a) = "A_{" ++ show a ++ "}"
  pretty (N n) = "{" ++ show n ++ "}"
  pretty (Product es) = unwords (map pretty es)
  pretty (Sum es) = intercalate (" + ") (map pretty es)

instance Num E where
  fromInteger = N
  e + f = Sum [e, f]
  e * f = Product [e, f]
  negate e = Product [N (-1), e]

instance Simplify E where
  simplify (Sum []) = 0
  simplify (Sum [e]) = simplify e
  simplify (Sum es) = case sortBy sumCmp (map simplify es) of
    N 0 : gs -> simplify $ Sum gs
    N i : N j : gs -> simplify $ Sum (N (i + j) : gs)
    Sum fs : gs -> simplify $ Sum (fs ++ gs)
    gs -> Sum gs
    where
      sumCmp (Sum a) (Sum b) = compare a b
      sumCmp (Sum _) _ = LT
      sumCmp _ (Sum _) = GT
      sumCmp (N a) (N b) = compare a b
      sumCmp (N _) _ = LT
      sumCmp _ (N _) = GT
      sumCmp a b = compare a b
  simplify (Product []) = 1
  simplify (Product [e]) = simplify e
  simplify (Product es) = case sortBy prodCmp (map simplify es) of
    N 0 : _ -> 0
    N 1 : gs -> simplify $ Product gs
    N i : N j : gs -> simplify $ Product (N (i * j) : gs)
    Sum fs : gs -> simplify $ Sum [ Product (f : gs) | f <- fs ]
    Product fs : gs -> simplify $ Product (fs ++ gs)
    gs -> Product gs
    where
      prodCmp (Product a) (Product b) = compare a b
      prodCmp (Product _) _ = LT
      prodCmp _ (Product _) = GT
      prodCmp (Sum a) (Sum b) = compare a b
      prodCmp (Sum _) _ = LT
      prodCmp _ (Sum _) = GT
      prodCmp (N a) (N b) = compare a b
      prodCmp (N _) _ = LT
      prodCmp _ (N _) = GT
      prodCmp a b = compare a b
  simplify e = e

instance Collect E where
  collect (Sum es) = Sum (concatMap accum . groupOn snd . sortOn snd . map single $ es)
    where
      single (Product (N n : fs)) = (n, fs)
      single (Product fs) = (1, fs)
      single e = (1, [e])
      accum ps = case sum $ map fst ps of
        0 -> []
        1 -> [Product (      snd (head ps))]
        n -> [Product (N n : snd (head ps))]
  collect e = e

perturb :: E -> E
perturb Z = Z + DZ
perturb C = C + DC
perturb (Sum es) = Sum (map perturb es)
perturb (Product es) = Product (map perturb es)
perturb e = e

perturbed :: E -> E
perturbed e = perturb e - e

series :: Integer -> E
series n = sum [ A i * DC ^ i | i <- [1 .. n] ]

sub :: E -> E -> E
sub e (Sum es) = Sum (map (sub e) es)
sub e (Product es) = Product (map (sub e) es)
sub e DZ = e
sub _ f = f

collate :: E -> [(Integer, E)]
collate (Sum es) = map accum . groupOn fst . sortOn fst . map single $ es
  where
    single (Product fs) = (toInteger (length gs), Product (reverse hs))
      where
        (gs, hs) = span (DC ==) . reverse $ fs
    single e = single (Product [e])
    accum ps = (fst (head ps), sum (map snd ps))
collate e = collate (Sum [e])


data CE = CSum [CE] | CProduct [CE] | CN Integer | ARe Integer | AIm Integer | ZRe | ZIm | CRe | CIm | DZRe | DZIm | DCRe | DCIm | I
  deriving (Read, Show, Eq, Ord)

instance Num CE where
  fromInteger = CN
  e + f = CSum [e, f]
  e * f = CProduct [e, f]
  negate e = CProduct [CN (-1), e]

instance Simplify CE where
  simplify (CSum []) = 0
  simplify (CSum [e]) = simplify e
  simplify (CSum es) = case sortBy sumCmp (map simplify es) of
    CN 0 : gs -> simplify $ CSum gs
    CN i : CN j : gs -> simplify $ CSum (CN (i + j) : gs)
    CSum fs : gs -> simplify $ CSum (fs ++ gs)
    gs -> CSum gs
    where
      sumCmp (CSum a) (CSum b) = compare a b
      sumCmp (CSum _) _ = LT
      sumCmp _ (CSum _) = GT
      sumCmp (CN a) (CN b) = compare a b
      sumCmp (CN _) _ = LT
      sumCmp _ (CN _) = GT
      sumCmp a b = compare a b
  simplify (CProduct []) = 1
  simplify (CProduct [e]) = simplify e
  simplify (CProduct es) = case sortBy prodCmp (map simplify es) of
    I : I : gs -> simplify $ CProduct (CN (-1) : gs)
    I : gs -> case simplify (CProduct gs) of
      CProduct hs -> CProduct (I : hs)
      h -> CProduct [I, h]
    CN 0 : _ -> 0
    CN 1 : gs -> simplify $ CProduct gs
    CN i : CN j : gs -> simplify $ CProduct (CN (i * j) : gs)
    CSum fs : gs -> simplify $ CSum [ CProduct (f : gs) | f <- fs ]
    CProduct fs : gs -> simplify $ CProduct (fs ++ gs)
    gs -> CProduct gs
    where
      prodCmp (CProduct a) (CProduct b) = compare a b
      prodCmp (CProduct _) _ = LT
      prodCmp _ (CProduct _) = GT
      prodCmp (CSum a) (CSum b) = compare a b
      prodCmp (CSum _) _ = LT
      prodCmp _ (CSum _) = GT
      prodCmp I I = EQ
      prodCmp I _ = LT
      prodCmp _ I = GT
      prodCmp (CN a) (CN b) = compare a b
      prodCmp (CN _) _ = LT
      prodCmp _ (CN _) = GT
      prodCmp a b = compare a b
  simplify e = e

instance Collect CE where
  collect (CSum es) = CSum (concatMap accum . groupOn snd . sortOn snd . map single $ es)
    where
      single (CProduct (CN n : fs)) = (n, fs)
      single (CProduct fs) = (1, fs)
      single e = (1, [e])
      accum ps = case sum $ map fst ps of
        0 -> []
        1 -> [CProduct (       snd (head ps))]
        n -> [CProduct (CN n : snd (head ps))]
  collect e = e

complex :: E -> CE
complex Z  = ZRe  + I * ZIm
complex C  = CRe  + I * CIm
complex DZ = DZRe + I * DZIm
complex DC = DCRe + I * DCIm
complex (A n) = ARe n + I * AIm n
complex (N n) = CN n
complex (Sum es) = CSum (map complex es)
complex (Product es) = CProduct (map complex es)

newtype Pair t = Pair{ unPair :: (t, t) }
  deriving (Read, Show, Eq, Ord)

instance Functor Pair where
  fmap f (Pair (a, b)) = Pair (f a, f b)

decomplex :: CE -> Pair CE
decomplex (CSum es)= Pair (CSum res, CSum ims)
  where (res, ims) = unzip (map (unPair . decomplex) es)
decomplex (CProduct es)
  | even n    = Pair (CProduct (i : fs), 0)
  | otherwise = Pair (0, CProduct (i : fs))
  where
    fs = filter (I /=) es
    n = length (filter (I ==) es)
    i = case n `mod` 4 of
      0 ->  1
      1 ->  1
      2 -> -1
      3 -> -1
decomplex I = Pair (0, 1)
decomplex e = Pair (e, 0)

cnormalize :: E -> Pair CE
cnormalize = fmap normalize . decomplex . normalize . complex . normalize


data OSum = OSum [OProduct]
  deriving (Read, Show, Eq, Ord)
data OProduct = OProduct Integer [OPower]
  deriving (Read, Show, Eq)
instance Ord OProduct where
  compare (OProduct x xs) (OProduct y ys) = compare (abs x) (abs y) <> compare (signum y) (signum x) <> compare xs ys
data OPower = OPower OVar Integer
  deriving (Read, Show, Eq, Ord)
data OVar = OVar Integer
  deriving (Read, Show, Eq, Ord)

osum :: CE -> OSum
osum (CSum cs) = OSum (sort . map oproduct $ cs)
osum c = OSum [oproduct c]

oproduct :: CE -> OProduct
oproduct (CProduct (CN n : cs)) = OProduct n (sort . map opower . group $ cs)
oproduct (CProduct cs) = OProduct 1 (sort . map opower . group $ cs)
oproduct (CN n) = OProduct n []
oproduct c = OProduct 1 [opower [c]]

opower :: [CE] -> OPower
opower cs@(c:_) = OPower (ovar c) (toInteger $ length cs)

ovar :: CE -> OVar
ovar CRe = OVar 0
ovar CIm = OVar 1
ovar ZRe = OVar 2
ovar ZIm = OVar 3
ovar (ARe n) = OVar $ 2 * n + 2
ovar (AIm n) = OVar $ 2 * n + 3
ovar e = error (show e)

data PSum = PSum [PMul1]
  deriving (Read, Show, Eq, Ord)
data PMul1 = PMul1 Integer [PMul2]
  deriving (Read, Show, Eq, Ord)
data PMul2 = PMul2 Bool Integer [PPower]
  deriving (Read, Show, Eq, Ord)
type PPower = OPower
type PVar = OVar

pmul1s :: PMul1 -> Integer
pmul1s (PMul1 i _) = i

pmul1p :: PMul1 -> [PMul2]
pmul1p (PMul1 _ p) = p

psum :: OSum -> PSum
psum (OSum ps) = PSum . map (\pms@(PMul1 i _ : _) ->  PMul1 i (concatMap pmul1p pms)) . groupOn pmul1s . sort . map (pmul1 0) $ ps

pmul1 :: Integer -> OProduct -> PMul1
pmul1 n (OProduct i os)
  | odd i = PMul1 (abs i) [PMul2 (i < 0) n os]
  | otherwise = pmul1 (n + 1) (OProduct (i `div` 2) os)

-- outermost first
data P
  = P'Add1 P P
  | P'Mul1 P Integer
  | P'Add2 P P
  | P'Neg P
  | P'Mul2 P Integer
  | P'Mul P P
  | P'Sqr P
  | P'Var Integer
  | P'Constant Integer
  deriving (Read, Show, Eq, Ord)

class Compile p where
  compile :: p -> P

instance Compile OVar where
  compile (OVar v) = P'Var v

instance Compile OPower where
  compile (OPower v 1) = compile v
  compile (OPower v n)
    | even n = P'Sqr (compile (OPower v (n `div` 2)))
    | otherwise = P'Mul (compile (OPower v (n `div` 2))) (compile (OPower v ((n `div` 2) + 1)))

instance Compile PMul2 where
  compile (PMul2 s i []) = P'Constant ((if s then negate else id) (2^i))
  compile (PMul2 s i xs) = (if s then P'Neg else id) . (if i == 0 then id else (`P'Mul2` i)) . compile $ xs

instance Compile [OPower] where
  compile [x] = compile x
  compile [x,y] = P'Mul (compile x) (compile y)
  compile xs = case splitAt (length xs `div` 2) xs of
    (ys, zs) -> P'Mul (compile ys) (compile zs)

instance Compile PMul1 where
  compile p@(PMul1 _ []) = error (show p)
  compile (PMul1 1 xs) = compile xs
  compile (PMul1 n xs) = P'Mul1 (compile xs) n

instance Compile [PMul2] where
  compile [x] = compile x
  compile [x,y] = P'Add2 (compile x) (compile y)
  compile xs = case splitAt (length xs `div` 2) xs of
    (ys, zs) -> P'Add2 (compile ys) (compile zs)

instance Compile PSum where
  compile p@(PSum []) = error (show p)
  compile (PSum xs) = compile xs

instance Compile [PMul1] where
  compile [x] = compile x
  compile [x,y] = P'Add1 (compile x) (compile y)
  compile xs = case splitAt (length xs `div` 2) xs of
    (ys, zs) -> P'Add1 (compile ys) (compile zs)

parallel :: [P] -> RWS () [Phase] (Map P Int) ()
parallel expressions = case groupOn pType . S.toList . S.fromList $ expressions of
  [] -> return ()
  (ops:opss) -> do
    results <- mapM temporary ops
    let subs = map subexpressions ops
    arguments <- mapM (mapM temporary) subs
    unless (null (concat arguments)) $
      tell [Phase (phaseOp . pType . head $ ops) (zip results arguments)]
    parallel (concat (subs ++ opss))

data Phase = Phase Op [(Either Integer Int, [Either Integer Int])]
  deriving (Read, Show, Eq, Ord)

data T = T'Add1 | T'Mul1 | T'Add2 | T'Neg | T'Mul2 | T'Mul | T'Sqr | T'Var | T'Constant
  deriving (Read, Show, Eq, Ord, Enum, Bounded)

pType :: P -> T
pType P'Add1{} = T'Add1
pType P'Mul1{} = T'Mul1
pType P'Add2{} = T'Add2
pType P'Neg{} = T'Neg
pType P'Mul2{} = T'Mul2
pType P'Mul{} = T'Mul
pType P'Sqr{} = T'Sqr
pType P'Var{} = T'Var
pType P'Constant{} = T'Constant

subexpressions :: P -> [P]
subexpressions (P'Add1 a b) = [a, b]
subexpressions (P'Mul1 a b) = [a, P'Constant b]
subexpressions (P'Add2 a b) = [a, b]
subexpressions (P'Neg a) = [a]
subexpressions (P'Mul2 a b) = [a, P'Constant b]
subexpressions (P'Mul a b) = [a, b]
subexpressions (P'Sqr a) = [a]
subexpressions _ = []

temporary :: Monoid w => P -> RWS r w (Map P Int) (Either Integer Int)
temporary (P'Constant n) = return (Left n)
temporary p = do
  m <- get
  case M.lookup p m of
    Nothing -> do
      let n = M.size m
      put (M.insert p n m)
      return (Right n)
    Just n -> return (Right n)

data Op = OpAdd | OpAddI | OpMulI | OpNeg | OpMul2 | OpMul | OpSqr | OpVar | OpSet
  deriving (Read, Show, Eq, Ord, Enum, Bounded)

phaseOp :: T -> Op
phaseOp T'Add1 = OpAdd
phaseOp T'Mul1 = OpMulI
phaseOp T'Add2 = OpAdd
phaseOp T'Neg  = OpNeg
phaseOp T'Mul2 = OpMul2
phaseOp T'Mul  = OpMul
phaseOp T'Sqr  = OpSqr

referenceCounts :: [Phase] -> Map Int Int
referenceCounts ps = M.fromListWith (+) [ (a, 1) | Phase _ ras <- ps, (_, args) <- ras, Right a <- args ]

duplicates :: Map Int Int -> Set Int
duplicates = M.keysSet . M.filter (1 <)

inplace :: [Phase] -> [Phase]
inplace phases = renames phases . fst $ execRWS (tell () >> unifies phases) () S.empty

compact :: [Phase] -> [Phase]
compact phases = renamem phases (M.fromList (zip (S.toList (M.keysSet (referenceCounts phases))) [0..]))

renames :: [Phase] -> Set (Set Int) -> [Phase]
renames phases renamings = map f phases
  where
    rens = S.toList renamings
    f (Phase op ras) = Phase op (map g ras)
    g (res, args) = (h res, map h args)
    h (Right r) = case filter (S.member r) rens of
      [] -> Right r
      [s] -> Right (F.minimum s) -- F. not needed with ghc-7.10
      ss -> error (show ss)
    h l = l

renamem :: [Phase] -> Map Int Int -> [Phase]
renamem phases renamingm = map f phases
  where
    f (Phase op ras) = Phase op (map g ras)
    g (res, args) = (h res, map h args)
    h (Right r) = Right (renamingm M.! r)
    h l = l

unify :: Monoid w => Set Int -> RWS r w (Set (Set Int)) ()
unify xs = do
  s <- get
  let (disjoint, intersecting) = S.partition (S.null . S.intersection xs) s
  put $ S.insert (S.unions (xs : S.toList intersecting)) disjoint

unifies :: Monoid w => [Phase] -> RWS r w (Set (Set Int)) ()
unifies phases = forM_ phases $ \(Phase op ras) -> case op of
    OpNeg  -> mapM_ go ras
    OpMul2 -> mapM_ go ras
    _ -> return ()
  where
    dups = duplicates . referenceCounts $ phases
    go (res, args)
      | any (`S.member` dups) ios = return ()
      | otherwise = unify (S.fromList ios)
      where
        ios = rights (res : args)

splitAdds :: [Phase] -> [Phase]
splitAdds = concatMap f
  where
    f (Phase OpAdd ras) = case partition (all isRight . snd) ras of
      ([], ai) -> [Phase OpAddI ai]
      (a, []) -> [Phase OpAdd a]
      (a, ai) -> [Phase OpAdd a, Phase OpAddI ai]
    f p = [p]


codegen :: Integer -> String -> [(Int, Phase)] -> [(FilePath, String)]
codegen order fname ps =

  [ (stem ++ ".c",
  "#include <complex.h>\n\
  \#include <limits.h>\n\
  \#include <math.h>\n\
  \#include <stdbool.h>\n\
  \#include <stdint.h>\n\
  \#include <stdio.h>\n\
  \#include <stdlib.h>\n\
  \#include <mpfr.h>\n\
  \\n\
  \#include " ++ show (stem ++ ".h") ++ "\n\
  \\n\
  \static inline int max(int a, int b) { return a > b ? a : b; }\n\
  \\n\
  \static const struct {\n" ++ unlines (map struct ps) ++ "  " ++ int ++ " v[" ++ show order ++ "][2];\n} " ++ stem ++ "_series_spec =\n\
  \{\n" ++ intercalate "\n,\n" (map values ps ++ [valids (last ps)]) ++ "};\n\
  \\n\
  \struct " ++ stem ++ "_series {\n\
  \  mpfr_t v[" ++ show count ++ "];\n\
  \  int n;\n\
  \};\n\
  \\n\
  \struct " ++ stem ++ "_series *" ++ stem ++ "_series_new(const mpfr_t cx, const mpfr_t cy) {\n\
  \  struct " ++ stem ++ "_series *s = malloc(sizeof(*s));\n\
  \  if (! s) { return 0; }\n\
  \  mpfr_prec_t p = max(mpfr_get_prec(cx), mpfr_get_prec(cy));\n\
  \  for (int i = 0; i < " ++ show count ++ "; ++i) {\n\
  \    mpfr_init2(s->v[i], p);\n\
  \    mpfr_set_si(s->v[i], 0, MPFR_RNDN);\n\
  \  };\n\
  \  mpfr_set(s->v[0], cx, MPFR_RNDN);\n\
  \  mpfr_set(s->v[1], cy, MPFR_RNDN);\n\
  \  mpfr_set(s->v[2], cx, MPFR_RNDN);\n\
  \  mpfr_set(s->v[3], cy, MPFR_RNDN);\n\
  \  mpfr_set_si(s->v[4], 1, MPFR_RNDN);\n\
  \  s->n = 1;\n\
  \  return s;\n\
  \}\n\
  \\n\
  \void " ++ stem ++ "_series_delete(struct " ++ stem ++ "_series *s) {\n\
  \  for (int i = 0; i < " ++ show count ++ "; ++i) {\n\
  \    mpfr_clear(s->v[i]);\n\
  \  }\n\
  \  free(s);\n\
  \}\n\
  \\n\
  \int " ++ stem ++ "_series_get_n(const struct " ++ stem ++ "_series *s) {\n\
  \  return s->n;\n\
  \}\n\
  \\n\
  \bool " ++ stem ++ "_series_step(struct " ++ stem ++ "_series *s, int exponent, int threshold) {\n\n" ++
  unlines (map genphase (init ps)) ++
  "  bool valid = true;\n\
  \  int e0;\n\
  \  int e1;\n\
  \  for (int i = 0; i < " ++ show order ++ "; ++i) {\n\
  \    e1 = INT_MIN;\n\
  \    for (int j = 0; j < 2; ++j) {\n\
  \      if (! mpfr_zero_p(s->v[" ++ stem ++ "_series_spec.v[i][j]])) {\n\
  \        e1 = max(e1, mpfr_get_exp(s->v[" ++ stem ++ "_series_spec.v[i][j]]));\n\
  \      }\n\
  \    }\n\
  \    if (i > 0) {\n\
  \      valid = e0 - exponent >= e1 + threshold;\n\
  \      e0 = max(e0 - exponent, e1);\n\
  \    } else {\n\
  \      e0 = e1;\n\
  \    }\n\
  \  }\n\
  \  if (! valid) { return false; }\n\
  \\n\
  \  s->n += 1;\n\n" ++
  unlines [genphase (last ps)] ++
  "  return true;\n\
  \}\n\
  \\n\
  \struct " ++ fname ++ "_reference *" ++ stem ++ "_reference_new(const struct " ++ stem ++ "_series *s) {\n\
  \  return " ++ fname ++ "_reference_new(s->v[0], s->v[1], s->v[2], s->v[3], s->n);\n\
  \}\n\
  \\n\
  \#define FTYPE float\n\
  \#define FNAME(x) x ## f\n\
  \#define FMPFRGET mpfr_get_flt\n\
  \#include \"" ++ stem ++ "_native.c\"\n\
  \#undef FTYPE\n\
  \#undef FNAME\n\
  \#undef FMPFRGET\n\
  \\n\
  \#define FTYPE double\n\
  \#define FNAME(x) x\n\
  \#define FMPFRGET mpfr_get_d\n\
  \#include \"" ++ stem ++ "_native.c\"\n\
  \#undef FTYPE\n\
  \#undef FNAME\n\
  \#undef FMPFRGET\n\
  \\n\
  \#define FTYPE long double\n\
  \#define FNAME(x) x ## l\n\
  \#define FMPFRGET mpfr_get_ld\n\
  \#include \"" ++ stem ++ "_native.c\"\n\
  \#undef FTYPE\n\
  \#undef FNAME\n\
  \#undef FMPFRGET\n\
  \\n")

  , (stem ++ "_native.c",
  "struct FNAME(" ++ stem ++ "_approx) {\n\
  \  complex FTYPE v[" ++ show order ++ "];\n\
  \  int exponent;\n\
  \};\n\
  \\n\
  \struct FNAME(" ++ stem ++ "_approx) *FNAME("  ++ stem ++ "_approx_new)(const struct " ++ stem ++ "_series *s, int exponent) {\n\
  \  struct FNAME(" ++ stem ++ "_approx) *a = malloc(sizeof(*a));\n\
  \  mpfr_t t;\n\
  \  mpfr_init2(t, mpfr_get_prec(s->v[0]));\n\
  \  for (int i = 0; i < " ++ show order ++ "; ++i) {\n\
  \    mpfr_set(t, s->v[2 * i + 4], MPFR_RNDN);\n\
  \    mpfr_mul_2si(t, t, (i + 1) * exponent, MPFR_RNDN);\n\
  \    FTYPE re = FMPFRGET(t, MPFR_RNDN);\n\
  \    mpfr_set(t, s->v[2 * i + 5], MPFR_RNDN);\n\
  \    mpfr_mul_2si(t, t, (i + 1) * exponent, MPFR_RNDN);\n\
  \    FTYPE im = FMPFRGET(t, MPFR_RNDN);\n\
  \    a->v[i] = re + I * im;\n\
  \  }\n\
  \  mpfr_clear(t);\n\
  \  a->exponent = -exponent;\n\
  \  return a;\n\
  \}\n\
  \\n\
  \complex FTYPE FNAME(" ++ stem ++ "_approx_do)(const struct FNAME(" ++ stem ++ "_approx) *a, complex FTYPE dc) {\n\
  \  complex FTYPE z = FNAME(ldexp)(FNAME(creal)(dc), a->exponent) + I * FNAME(ldexp)(FNAME(cimag)(dc), a->exponent);\n\
  \  complex FTYPE zi = z;\n\
  \  complex FTYPE s = 0;\n\
  \  for (int i = 0; i < " ++ show order ++ "; ++i) {\n\
  \    s += a->v[i] * zi;\n\
  \    zi *= z;\n\
  \  }\n\
  \  return s;\n\
  \}\n\
  \\n\
  \int FNAME(" ++ stem ++ "_approx_get_exponent)(const struct FNAME(" ++ stem ++ "_approx) *a) {\n\
  \  return a->exponent;\n\
  \}\n\
  \const complex FTYPE *FNAME(" ++ stem ++ "_approx_get_coefficients)(const struct FNAME(" ++ stem ++ "_approx) *a) {\n\
  \  return &a->v[0];\n\
  \}\n\
  \")

  , (stem ++ ".h",
  "#ifndef " ++ stem ++ "_h\n\
  \#define " ++ stem ++ "_h 1\n\
  \\n\
  \#include <complex.h>\n\
  \#include <mpfr.h>\n\
  \\n\
  \#include " ++ show (fname ++ "_ref.h") ++ "\n\
  \\n\
  \struct " ++ stem ++ "_series;\n\
  \struct " ++ stem ++ "_series *" ++ stem ++ "_series_new(const mpfr_t cx, const mpfr_t cy);\n\
  \void " ++ stem ++ "_series_delete(struct " ++ stem ++ "_series *s);\n\
  \int " ++ stem ++ "_series_get_n(const struct " ++ stem ++ "_series *s);\n\
  \bool " ++ stem ++ "_series_step(struct " ++ stem ++ "_series *s, int exponent, int threshold);\n\
  \struct " ++ fname ++ "_reference *" ++ stem ++ "_reference_new(const struct " ++ stem ++ "_series *s);\n\
  \struct " ++ stem ++ "_approxf;\n\
  \struct " ++ stem ++ "_approx;\n\
  \struct " ++ stem ++ "_approxl;\n\
  \struct " ++ stem ++ "_approxf *"  ++ stem ++ "_approx_newf(const struct " ++ stem ++ "_series *s, int exponent);\n\
  \struct " ++ stem ++ "_approx *"  ++ stem ++ "_approx_new(const struct " ++ stem ++ "_series *s, int exponent);\n\
  \struct " ++ stem ++ "_approxl *"  ++ stem ++ "_approx_newl(const struct " ++ stem ++ "_series *s, int exponent);\n\
  \void " ++ stem ++ "_approx_deletef(struct " ++ stem ++ "_approxf *s);\n\
  \void " ++ stem ++ "_approx_delete(struct " ++ stem ++ "_approx *s);\n\
  \void " ++ stem ++ "_approx_deletel(struct " ++ stem ++ "_approxl *s);\n\
  \int " ++ stem ++ "_approx_get_exponentf(const struct " ++ stem ++ "_approxf *a);\n\
  \int " ++ stem ++ "_approx_get_exponent(const struct " ++ stem ++ "_approx *a);\n\
  \int " ++ stem ++ "_approx_get_exponentl(const struct " ++ stem ++ "_approxl *a);\n\
  \const complex float *" ++ stem ++ "_approx_get_coefficientsf(const struct " ++ stem ++ "_approxf *a);\n\
  \const complex double *" ++ stem ++ "_approx_get_coefficients(const struct " ++ stem ++ "_approx *a);\n\
  \const complex long double *" ++ stem ++ "_approx_get_coefficientsl(const struct " ++ stem ++ "_approxl *a);\n\
  \complex float " ++ stem ++ "_approx_dof(const struct " ++ stem ++ "_approxf *a, complex float dc);\n\
  \complex double " ++ stem ++ "_approx_do(const struct " ++ stem ++ "_approx *a, complex double dc);\n\
  \complex long double " ++ stem ++ "_approx_dol(const struct " ++ stem ++ "_approxl *a, complex long double dc);\n\
  \\n\
  \#endif\n")]

  where
    int
      | count <= 2^( 8 :: Int) = "uint8_t"
      | count <= 2^(16 :: Int) = "uint16_t"
      | count <= 2^(32 :: Int) = "uint32_t"
      | otherwise     = "uint64_t"
    stem = fname ++ "_" ++ show order
    struct (i, Phase _ ras@((_, args):_)) = "  " ++ int ++ " p" ++ show i ++ "[" ++ show (length ras) ++ "][" ++ show (length args + 1) ++ "];"
    values (_, Phase _ ras) = "{\n" ++ intercalate ",\n" (map value ras) ++ "\n}"
    value (Right res, args) = "{ " ++ show res ++ " , " ++ intercalate " , " (map show (rights args) ++ map show (lefts args)) ++ " }"
    count = 1 + maximum [ i | (_, Phase _ ras) <- ps, (res, args) <-ras, Right i <- res : args ]
    valids (_, Phase OpSet ras)
      = "{\n" ++ intercalate ",\n"
          [ "{ " ++ show re ++ " , " ++ show im ++ " }"
          | o <- [0 .. fromIntegral order - 1]
          , (Right sre, [ Right re ]) <- ras
          , sre == 2 * o + 4
          , (Right sim, [ Right im ]) <- ras
          , sim == 2 * o + 5
          ] ++ "\n}\n"
    genphase (p, Phase op ras) =
      "  #pragma omp parallel for\n\
      \  for (int i = 0; i < " ++ show (length ras) ++ "; ++i) {\n" ++
      ( case op of
          OpAdd  -> o3  "mpfr_add"
          OpAddI -> o3i "mpfr_add_si"
          OpMulI -> o3i "mpfr_mul_si"
          OpNeg  -> o2  "mpfr_neg"
          OpMul2 -> o3i "mpfr_mul_2si"
          OpMul  -> o3  "mpfr_mul"
          OpSqr  -> o2  "mpfr_sqr"
          OpSet  -> o2  "mpfr_set"
          _ -> error (show op) ) ++
      "  }\n"
      where
        o3  s = "  " ++ s ++ "\n\
                \    ( s->v[" ++ stem ++ "_series_spec.p" ++ show p ++ "[i][0]]\n\
                \    , s->v[" ++ stem ++ "_series_spec.p" ++ show p ++ "[i][1]]\n\
                \    , s->v[" ++ stem ++ "_series_spec.p" ++ show p ++ "[i][2]]\n\
                \    , MPFR_RNDN\n\
                \    );\n"
        o3i s = "  " ++ s ++ "\n\
                \    ( s->v[" ++ stem ++ "_series_spec.p" ++ show p ++ "[i][0]]\n\
                \    , s->v[" ++ stem ++ "_series_spec.p" ++ show p ++ "[i][1]]\n\
                \    ,      " ++ stem ++ "_series_spec.p" ++ show p ++ "[i][2]\n\
                \    , MPFR_RNDN\n\
                \    );\n"
        o2  s = "  " ++ s ++ "\n\
                \    ( s->v[" ++ stem ++ "_series_spec.p" ++ show p ++ "[i][0]]\n\
                \    , s->v[" ++ stem ++ "_series_spec.p" ++ show p ++ "[i][1]]\n\
                \    , MPFR_RNDN\n\
                \    );\n"

main' :: String -> [Integer] -> E -> [(FilePath, String)]
main' stem orders f = main''' stem f ++ concat [main'' stem order f | order <- orders ] ++ codegenWrap stem orders ++ codegenMake stem orders

main''' :: String -> E -> [(FilePath, String)]
main''' stem f = codegenRef stem ref
  where
    ref = zip [0..] . splitAdds . compact . inplace . reverse $ finalize : phases
      where
        finalize = Phase OpSet $ zip (map (Right . (variables M.!)) vs) (map ((:[]) . Right . (variables M.!)) es)
        initial = fst $ execRWS initialize () M.empty
        initialize = tell () >> mapM_ temporary ([var CRe, var CIm] ++ vs ++ es)
        (variables, phases) = execRWS (parallel es) () initial
        (vs, es) = unzip [ (var v, expr e) | (v, e) <- ves ]
        var = compile . ovar
        expr = compile . psum . osum
        ves =
          [ (ZRe, fre)
          , (ZIm, fim)
          ]
        Pair (fre, fim) = cnormalize $ f

codegenMake :: String -> [Integer] -> [(FilePath, String)]
codegenMake stem orders =
  [ ( "Makefile",
  "OBJECTS = \\\n" ++
  unlines [ "\t" ++ stem ++ "_" ++ show order ++ ".o \\" | order <- orders ] ++
  "\t" ++ stem ++ "_ref.o \\\n\t" ++ stem ++ ".o\n\
  \\n\
  \all: $(OBJECTS)\n\
  \.SUFFIXES:\n\
  \.PHONY: all\n\
  \%.o: %.c\n\
  \\tgcc -std=c99 -pedantic -Wall -Wextra -fopenmp -O3 -march=native -c $<\n\
  \\n")
  ]

codegenWrap :: String -> [Integer] -> [(FilePath, String)]
codegenWrap stem orders =
  [ ( stem ++ ".h",
  "#ifndef " ++ stem ++ "_h\n\
  \#define " ++ stem ++ "_h 1\n\
  \\n\
  \#include <complex.h>\n\
  \#include <stdbool.h>\n\
  \\n\
  \#include " ++ show (stem ++ "_ref.h") ++ "\n\
  \struct " ++ stem ++ "_series;\n\
  \struct " ++ stem ++ "_series *" ++ stem ++ "_series_new(const int order, const mpfr_t cx, const mpfr_t cy);\n\
  \void " ++ stem ++ "_series_delete(struct " ++ stem ++ "_series *s);\n\
  \bool " ++ stem ++ "_series_step(struct " ++ stem ++ "_series *s, const int exponent, const int threshold);\n\
  \int " ++ stem ++ "_series_get_n(const struct " ++ stem ++ "_series *s);\n\
  \struct " ++ stem ++ "_reference *" ++ stem ++ "_series_reference_new(const struct " ++ stem ++ "_series *s);\n\
  \\n\
  \#define FTYPE float\n\
  \#define FNAME(x) x ## f\n\
  \#include \"" ++ stem ++ "_approx_native.h\"\n\
  \#undef FTYPE\n\
  \#undef FNAME\n\
  \\n\
  \#define FTYPE double\n\
  \#define FNAME(x) x\n\
  \#include \"" ++ stem ++ "_approx_native.h\"\n\
  \#undef FTYPE\n\
  \#undef FNAME\n\
  \\n\
  \#define FTYPE long double\n\
  \#define FNAME(x) x ## l\n\
  \#include \"" ++ stem ++ "_approx_native.h\"\n\
  \#undef FTYPE\n\
  \#undef FNAME\n\
  \\n\
  \#endif\n\
  \\n")

  , ( stem ++ "_approx_native.h",
  "struct FNAME(" ++ stem ++ "_approx);\n\
  \struct FNAME(" ++ stem ++ "_approx) *FNAME(" ++ stem ++ "_series_approx_new)(const struct " ++ stem ++ "_series *s, const int exponent);\n\
  \void FNAME(" ++ stem ++ "_approx_delete)(struct FNAME(" ++ stem ++ "_approx) *a);\n\
  \int FNAME(" ++ stem ++ "_approx_get_order)(const struct FNAME(" ++ stem ++ "_approx) *a);\n\
  \int FNAME(" ++ stem ++ "_approx_get_exponent)(const struct FNAME(" ++ stem ++ "_approx) *a);\n\
  \const complex FTYPE *FNAME(" ++ stem ++ "_approx_get_coefficients)(const struct FNAME(" ++ stem ++ "_approx) *a);\n\
  \complex FTYPE FNAME(" ++ stem ++ "_approx_do)(const struct FNAME(" ++ stem ++ "_approx) *a, const complex FTYPE dc);\n\
  \\n")

  , ( stem ++ ".c",
  "#include <stdlib.h>\n\
  \\n\
  \#include \"" ++ stem ++ ".h\"\n\
  \\n" ++
  unlines [ "#include " ++ show (stem ++ "_" ++ show order ++ ".h") | order <- orders ] ++
  "\n\
  \struct " ++ stem ++ "_series {\n\
  \  int order;\n\
  \  void *series;\n\
  \};\n\
  \\n\
  \struct " ++ stem ++ "_series *" ++ stem ++ "_series_new(const int order, const mpfr_t cx, const mpfr_t cy) {\n\
  \   struct " ++ stem ++ "_series *s = calloc(1, sizeof(*s));\n\
  \  if (! s) { return 0; }\n\
  \  s->order = order;\n\
  \  switch (order) {\n" ++
  unlines [ "    case " ++ show order ++ ": s->series = " ++ stem ++ "_" ++ show order ++ "_series_new(cx, cy); break;" | order <- orders ] ++
  "  }\n\
  \  if (! s->series) {\n\
  \    free(s);\n\
  \    return 0;\n\
  \  }\n\
  \  return s;\n\
  \}\n\
  \\n\
  \void " ++ stem ++ "_series_delete(struct " ++ stem ++ "_series *s) {\n\
  \  if (! s) {\n\
  \    return;\n\
  \  }\n\
  \  switch (s->order) {\n" ++
  unlines [ "    case " ++ show order ++ ": " ++ stem ++ "_" ++ show order ++ "_series_delete(s->series); break;" | order <- orders ] ++
  "  }\n\
  \  free(s);\n\
  \}\n\
  \\n\
  \bool " ++ stem ++ "_series_step(struct " ++ stem ++ "_series *s, const int exponent, const int threshold) {\n\
  \  if (! s) {\n\
  \    return false;\n\
  \  }\n\
  \  switch (s->order) {\n" ++
  unlines [ "    case " ++ show order ++ ": return " ++ stem ++ "_" ++ show order ++ "_series_step(s->series, exponent, threshold);" | order <- orders ] ++
  "  }\n\
  \  return false;\n\
  \}\n\
  \\n\
  \int " ++ stem ++ "_series_get_n(const struct " ++ stem ++ "_series *s) {\n\
  \  if (! s) {\n\
  \    return 0;\n\
  \  }\n\
  \  switch (s->order) {\n" ++
  unlines [ "    case " ++ show order ++ ": return " ++ stem ++ "_" ++ show order ++ "_series_get_n(s->series);" | order <- orders ] ++
  "  }\n\
  \  return 0;\n\
  \}\n\
  \\n\
  \struct " ++ stem ++ "_reference *" ++ stem ++ "_series_reference_new(const struct " ++ stem ++ "_series *s) {\n\
  \  if (! s) {\n\
  \    return 0;\n\
  \  }\n\
  \  switch (s->order) {\n" ++
  unlines [ "    case " ++ show order ++ ": return " ++ stem ++ "_" ++ show order ++ "_reference_new(s->series);" | order <- orders ] ++
  " }\n\
  \  return 0;\n\
  \}\n\
  \\n\
  \#define FTYPE float\n\
  \#define FNAME(x) x ## f\n\
  \#include \"" ++ stem ++ "_approx_native.c\"\n\
  \#undef FTYPE\n\
  \#undef FNAME\n\
  \\n\
  \#define FTYPE double\n\
  \#define FNAME(x) x\n\
  \#include \"" ++ stem ++ "_approx_native.c\"\n\
  \#undef FTYPE\n\
  \#undef FNAME\n\
  \\n\
  \#define FTYPE long double\n\
  \#define FNAME(x) x ## l\n\
  \#include \"" ++ stem ++ "_approx_native.c\"\n\
  \#undef FTYPE\n\
  \#undef FNAME\n\
  \\n")

  , ( stem ++ "_approx_native.c",
  "struct FNAME(" ++ stem ++ "_approx) {\n\
  \  int order;\n\
  \  void *approx;\n\
  \};\n\
  \\n\
  \struct FNAME(" ++ stem ++ "_approx) *FNAME(" ++ stem ++ "_series_approx_new)(const struct " ++ stem ++ "_series *s, const int exponent) {\n\
  \  if (! s) {\n\
  \    return 0;\n\
  \  }\n\
  \  struct FNAME(" ++ stem ++ "_approx) *a = calloc(1, sizeof(*a));\n\
  \  a->order = s->order;\n\
  \  switch (s->order) {\n" ++
  unlines [ "    case " ++ show order ++ ": a->approx = FNAME(" ++ stem ++ "_" ++ show order ++ "_approx_new)(s->series, exponent); break;" | order <- orders ] ++
  "  }\n\
  \  if (! a->approx) {\n\
  \    free(a);\n\
  \    return 0;\n\
  \  }\n\
  \  return a;\n\
  \}\n\
  \\n\
  \void FNAME(" ++ stem ++ "_approx_delete)(struct FNAME(" ++ stem ++ "_approx) *a) {\n\
  \  if (! a) {\n\
  \    return;\n\
  \  }\n\
  \  free(a->approx);  // FIXME check this\n\
  \  free(a);\n\
  \}\n\
  \\n\
  \int FNAME(" ++ stem ++ "_approx_get_order)(const struct FNAME(" ++ stem ++ "_approx) *a) {\n\
  \  if (! a) {\n\
  \    return 0;\n\
  \  }\n\
  \  return a->order;\n\
  \}\n\
  \\n\
  \int FNAME(" ++ stem ++ "_approx_get_exponent)(const struct FNAME(" ++ stem ++ "_approx) *a) {\n\
  \  if (! a) {\n\
  \    return 0;\n\
  \  }\n\
  \  switch (a->order) {\n" ++
  unlines [ "    case " ++ show order ++ ": return FNAME(" ++ stem ++ "_" ++ show order ++ "_approx_get_exponent)(a->approx);" | order <- orders ] ++
  "  }\n\
  \  return 0;\n\
  \}\n\
  \\n\
  \const complex FTYPE *FNAME(" ++ stem ++ "_approx_get_coefficients)(const struct FNAME(" ++ stem ++ "_approx) *a) {\n\
  \  if (! a) {\n\
  \    return 0;\n\
  \  }\n\
  \  switch (a->order) {" ++
  unlines [ "    case " ++ show order ++ ": return FNAME(" ++ stem ++ "_" ++ show order ++ "_approx_get_coefficients)(a->approx);" | order <- orders ] ++
  "  }\n\
  \  return 0;\n\
  \}\n\
  \\n\
  \complex FTYPE FNAME(" ++ stem ++ "_approx_do)(const struct FNAME(" ++ stem ++ "_approx) *a, const complex FTYPE dc) {\n\
  \  if (! a) {\n\
  \    return 0;\n\
  \  }\n\
  \  switch (a->order) {\n" ++
  unlines [ "    case " ++ show order ++ ": return FNAME(" ++ stem ++ "_" ++ show order ++ "_approx_do)(a->approx, dc);" | order <- orders ] ++
  "  }\n\
  \  return 0;\n\
  \}\n\
  \\n")
  ]

codegenRef :: String -> [(Int, Phase)] -> [(FilePath, String)]
codegenRef stem ps =
  [(stem ++ "_ref.c",
  "#include <complex.h>\n\
  \#include <limits.h>\n\
  \#include <math.h>\n\
  \#include <stdbool.h>\n\
  \#include <stdint.h>\n\
  \#include <stdlib.h>\n\
  \#include <mpfr.h>\n\
  \\n\
  \#include " ++ show (stem ++ "_ref.h") ++ "\n\
  \\n\
  \static inline int max(int a, int b) { return a > b ? a : b; }\n\
  \\n\
  \static const struct {\n" ++ unlines (map struct ps) ++ "} " ++ stem ++ "_reference_spec =\n\
  \{\n" ++ intercalate "\n,\n" (map values ps) ++ "};\n\
  \\n\
  \struct " ++ stem ++ "_reference {\n\
  \  mpfr_t v[" ++ show count ++ "];\n\
  \  int n;\n\
  \};\n\
  \\n\
  \struct " ++ stem ++ "_reference *" ++ stem ++ "_reference_new(const mpfr_t cx, const mpfr_t cy, const mpfr_t zx, const mpfr_t zy, int n) {\n\
  \  struct " ++ stem ++ "_reference *r = malloc(sizeof(*r));\n\
  \  if (! r) { return 0; }\n\
  \  mpfr_prec_t p = max(max(mpfr_get_prec(cx), mpfr_get_prec(cy)), max(mpfr_get_prec(zx), mpfr_get_prec(zy)));\n\
  \  for (int i = 0; i < " ++ show count ++ "; ++i) {\n\
  \    mpfr_init2(r->v[i], p);\n\
  \    mpfr_set_si(r->v[i], 0, MPFR_RNDN);\n\
  \  };\n\
  \  mpfr_set(r->v[0], cx, MPFR_RNDN);\n\
  \  mpfr_set(r->v[1], cy, MPFR_RNDN);\n\
  \  mpfr_set(r->v[2], zx, MPFR_RNDN);\n\
  \  mpfr_set(r->v[3], zy, MPFR_RNDN);\n\
  \  r->n = n;\n\
  \  return r;\n\
  \}\n\
  \\n\
  \void " ++ stem ++ "_reference_delete(struct " ++ stem ++ "_reference *r) {\n\
  \  for (int i = 0; i < " ++ show count ++ "; ++i) {\n\
  \    mpfr_clear(r->v[i]);\n\
  \  }\n\
  \  free(r);\n\
  \}\n\
  \\n\
  \void " ++ stem ++ "_reference_step(struct " ++ stem ++ "_reference *r) {\n\n" ++
  unlines (map genphase ps) ++
  "  r->n += 1;\n\
  \}\n\
  \\n\
  \int " ++ stem ++ "_reference_get_n(const struct " ++ stem ++ "_reference *r) {\n\
  \  return r->n;\n\
  \}\n\
  \\n\
  \complex float " ++ stem ++ "_reference_get_zf(const struct " ++ stem ++ "_reference *r) {\n\
  \  return mpfr_get_flt(r->v[2], MPFR_RNDN) + I * mpfr_get_flt(r->v[3], MPFR_RNDN);\n\
  \}\n\
  \\n\
  \complex double " ++ stem ++ "_reference_get_z(const struct " ++ stem ++ "_reference *r) {\n\
  \  return mpfr_get_d(r->v[2], MPFR_RNDN) + I * mpfr_get_d(r->v[3], MPFR_RNDN);\n\
  \}\n\
  \\n\
  \complex long double " ++ stem ++ "_reference_get_zl(const struct " ++ stem ++ "_reference *r) {\n\
  \  return mpfr_get_ld(r->v[2], MPFR_RNDN) + I * mpfr_get_ld(r->v[3], MPFR_RNDN);\n\
  \}\n\
  \\n\
  \void " ++ stem ++ "_reference_get_zr(const struct " ++ stem ++ "_reference *ref, mpfr_t zx, mpfr_t zy) {\n\
  \  mpfr_set_prec(zx, mpfr_get_prec(ref->v[2]));\n\
  \  mpfr_set_prec(zy, mpfr_get_prec(ref->v[3]));\n\
  \  mpfr_set(zx, ref->v[2], MPFR_RNDN);\n\
  \  mpfr_set(zy, ref->v[3], MPFR_RNDN);\n\
  \}\n\
  \\n")

  , (stem ++ "_ref.h",
  "#ifndef " ++ stem ++ "_ref_h\n\
  \#define " ++ stem ++ "_ref_h 1\n\
  \\n\
  \#include <complex.h>\n\
  \#include <mpfr.h>\n\
  \\n\
  \struct " ++ stem ++ "_reference;\n\
  \struct " ++ stem ++ "_reference *" ++ stem ++ "_reference_new(const mpfr_t cx, const mpfr_t cy, const mpfr_t zx, const mpfr_t zy, int n);\n\
  \void " ++ stem ++ "_reference_delete(struct " ++ stem ++ "_reference *r);\n\
  \void " ++ stem ++ "_reference_step(struct " ++ stem ++ "_reference *r);\n\
  \int " ++ stem ++ "_reference_get_n(const struct " ++ stem ++ "_reference *r);\n\
  \complex float " ++ stem ++ "_reference_get_zf(const struct " ++ stem ++ "_reference *r);\n\
  \complex double " ++ stem ++ "_reference_get_z(const struct " ++ stem ++ "_reference *r);\n\
  \complex long double " ++ stem ++ "_reference_get_zl(const struct " ++ stem ++ "_reference *r);\n\
  \void " ++ stem ++ "_reference_get_zr(const struct " ++ stem ++ "_reference *ref, mpfr_t zx, mpfr_t zy);\n\
  \\n\
  \#endif\n")]

  where
    int
      | count <= 2^( 8 :: Int) = "uint8_t"
      | count <= 2^(16 :: Int) = "uint16_t"
      | count <= 2^(32 :: Int) = "uint32_t"
      | otherwise     = "uint64_t"
    struct (i, Phase _ ras@((_, args):_)) = "  " ++ int ++ " p" ++ show i ++ "[" ++ show (length ras) ++ "][" ++ show (length args + 1) ++ "];"
    values (_, Phase _ ras) = "{\n" ++ intercalate ",\n" (map value ras) ++ "\n}"
    value (Right res, args) = "{ " ++ show res ++ " , " ++ intercalate " , " (map show (rights args) ++ map show (lefts args)) ++ " }"
    count = 1 + maximum [ i | (_, Phase _ ras) <- ps, (res, args) <-ras, Right i <- res : args ]

    genphase (p, Phase op ras) =
      "  #pragma omp parallel for\n\
      \  for (int i = 0; i < " ++ show (length ras) ++ "; ++i) {\n" ++
      ( case op of
          OpAdd  -> o3  "mpfr_add"
          OpAddI -> o3i "mpfr_add_si"
          OpMulI -> o3i "mpfr_mul_si"
          OpNeg  -> o2  "mpfr_neg"
          OpMul2 -> o3i "mpfr_mul_2si"
          OpMul  -> o3  "mpfr_mul"
          OpSqr  -> o2  "mpfr_sqr"
          OpSet  -> o2  "mpfr_set"
          _ -> error (show op) ) ++
      "  }\n"
      where
        o3  s = "  " ++ s ++ "\n\
                \    ( r->v[" ++ stem ++ "_reference_spec.p" ++ show p ++ "[i][0]]\n\
                \    , r->v[" ++ stem ++ "_reference_spec.p" ++ show p ++ "[i][1]]\n\
                \    , r->v[" ++ stem ++ "_reference_spec.p" ++ show p ++ "[i][2]]\n\
                \    , MPFR_RNDN\n\
                \    );\n"
        o3i s = "  " ++ s ++ "\n\
                \    ( r->v[" ++ stem ++ "_reference_spec.p" ++ show p ++ "[i][0]]\n\
                \    , r->v[" ++ stem ++ "_reference_spec.p" ++ show p ++ "[i][1]]\n\
                \    ,      " ++ stem ++ "_reference_spec.p" ++ show p ++ "[i][2]\n\
                \    , MPFR_RNDN\n\
                \    );\n"
        o2  s = "  " ++ s ++ "\n\
                \    ( r->v[" ++ stem ++ "_reference_spec.p" ++ show p ++ "[i][0]]\n\
                \    , r->v[" ++ stem ++ "_reference_spec.p" ++ show p ++ "[i][1]]\n\
                \    , MPFR_RNDN\n\
                \    );\n"

main'' :: String -> Integer -> E -> [(FilePath, String)]
main'' stem n f = codegen n stem . zip [0..] . splitAdds . compact . inplace . reverse $ finalize : phases
  where
    finalize = Phase OpSet $ zip (map (Right . (variables M.!)) vs) (map ((:[]) . Right . (variables M.!)) es)
    initial = fst $ execRWS initialize () M.empty
    initialize = tell () >> mapM_ temporary ([var CRe, var CIm] ++ vs ++ es)
    (variables, phases) = execRWS (parallel es) () initial
    (vs, es) = unzip [ (var v, expr e) | (v, e) <- ves ]
    var = compile . ovar
    expr = compile . psum . osum
    ves =
      [ (ZRe, fre)
      , (ZIm, fim)
      ] ++ concat
      [ [ (ARe i, are)
        , (AIm i, aim)
        ]
      | (i, Pair (are, aim)) <- ies
      ]
    Pair (fre, fim) = cnormalize $ f
    ies
      = take (fromInteger n)
      . map (fmap cnormalize)
      . collate
      . normalize
      . sub (series n)
      . normalize
      . perturbed
      $ f

main :: IO ()
main = mapM_ (uncurry writeFile) . main' "z2c" [4,6,8,12,16,24,32,48,64] $ Z^(2::Int) + C
