-- SPDX-License-Identifier: BSD-3-Clause
-- | myb-generate: writes the lower layer of the Haskell binding from the public headers of myb.
--
-- The headers are run through the preprocessor of Clang, which keeps its line markers, so each
-- declaration is known by the header it comes from. They are read into types and functions,
-- and each header becomes a module of hsc2hs, Myb.C.<Header> (Myb.C.Ext.<Name> for ext/):
--
--   * a structure is a record with a Storable instance whose size, alignment and field offsets
--     hsc2hs takes from the header itself when the package is built: nothing is laid out here;
--   * an enumeration is a newtype over the integer type hsc2hs finds for it, its constants
--     pattern synonyms by their C names with the values of the header; so are the constants of
--     the anonymous enumerations and of #define;
--   * a function is a foreign import by its C name, and p'<name> its address; a type of
--     function pointer is a FunPtr, with mk'<name> ("wrapper") and call'<name> ("dynamic");
--   * what a foreign import cannot pass, a structure by value, goes through a C function
--     written into generated/cbits that takes a pointer instead and only copies the structure.
--
-- The text is the same on every system: it holds no path, no size, no offset. Anything the
-- reader does not know stops it with a message; nothing is skipped in silence.
--
-- usage: myb-generate [--check] [--clang <clang>] <repository root>
--   --check  write nothing: exit with 1 if the files of the package differ from what the
--            headers give now
--   --clang  the preprocessor (by default $MYB_CLANG, else clang)
module Main (main) where

import Control.Monad (foldM, forM, forM_, unless, when)
import Data.Char (isAlpha, isAlphaNum, isDigit, isSpace, toUpper)
import Data.List (intercalate, isInfixOf, isPrefixOf, isSuffixOf, nub, sort)
import qualified Data.Map.Strict as Map
import Data.Maybe (fromMaybe, isNothing, mapMaybe)
import qualified Data.Set as Set
import System.Directory (createDirectoryIfMissing, doesDirectoryExist, doesFileExist, listDirectory, removeDirectory, removeFile)
import System.Environment (getArgs, lookupEnv)
import System.Exit (ExitCode (..), die, exitWith)
import System.FilePath (dropExtension, takeDirectory, (</>))
import System.IO (IOMode (..), hGetContents, hPutStr, hPutStrLn, stderr, withBinaryFile)
import System.Process (readProcessWithExitCode)

-- ---- C declarations ------------------------------------------------------------------------

data CType
  = TBase String -- a scalar by its C spelling: "double", "unsigned int", "size_t", ...
  | TName String -- a type of the library: myb_...
  | TConst CType
  | TPtr CType
  | TArray CType String -- the element and the length as the header writes it
  | TFunc CType [Param]
  deriving (Eq, Show)

data Param = Param {paramName :: Maybe String, paramType :: CType} deriving (Eq, Show)

data Field = Field {fieldName :: String, fieldType :: CType} deriving (Eq, Show)

data Decl
  = DStruct String (Maybe String) (Maybe [Field]) -- typedef struct [tag] [{...}] name
  | DTagBody String [Field] -- struct tag {...}: the body of a structure declared before
  | DEnum (Maybe String) [String] -- typedef enum ... name, or enum {...} (constants)
  | DTypedef String CType
  | DFunction String CType [Param]

strip :: CType -> CType
strip (TConst t) = strip t
strip t = t

-- ---- Tokens --------------------------------------------------------------------------------

data Token = Ident String | Number String | Punct String deriving (Eq, Show)

tokenText :: Token -> String
tokenText (Ident s) = s
tokenText (Number s) = s
tokenText (Punct s) = s

tokenize :: String -> Either String [Token]
tokenize [] = Right []
tokenize s@(c : cs)
  | isSpace c = tokenize cs
  | isAlpha c || c == '_' = let (w, rest) = span isWord s in (Ident w :) <$> tokenize rest
  | isDigit c = let (w, rest) = number s in (Number w :) <$> tokenize rest
  | Just p <- punctuation s = (Punct p :) <$> tokenize (drop (length p) s)
  | otherwise = Left ("a character the reader does not know: " ++ show c)
  where
    isWord x = isAlphaNum x || x == '_'

number :: String -> (String, String)
number (x : y : rest)
  | x `elem` "eEpP" && y `elem` "+-" = let (a, b) = number rest in (x : y : a, b)
number (x : rest)
  | isAlphaNum x || x == '.' || x == '_' = let (a, b) = number rest in (x : a, b)
number rest = ("", rest)

punctuation :: String -> Maybe String
punctuation s =
  case filter (`isPrefixOf` s) (["...", "<<", ">>", "->", "&&", "||", "==", "!=", "<=", ">="] ++ map pure "(){}[];,*=+-/%<>&|^~!?:.") of
    (p : _) -> Just p
    [] -> Nothing

-- ---- A parser of declarations --------------------------------------------------------------

newtype P a = P {runP :: [Token] -> Either String (a, [Token])}

instance Functor P where
  fmap f (P g) = P $ \ts -> fmap (\(a, rest) -> (f a, rest)) (g ts)

instance Applicative P where
  pure a = P $ \ts -> Right (a, ts)
  pf <*> pa = pf >>= \f -> fmap f pa

instance Monad P where
  P g >>= k = P $ \ts -> g ts >>= \(a, rest) -> runP (k a) rest

failP :: String -> P a
failP message = P $ \ts -> Left (message ++ ", before: " ++ unwords (map tokenText (take 12 ts)))

lookahead :: Int -> P [Token]
lookahead n = P $ \ts -> Right (take n ts, ts)

next :: P Token
next = P $ \ts -> case ts of
  [] -> Left "the declaration ends too early"
  (t : rest) -> Right (t, rest)

skip :: P ()
skip = () <$ next

optSym :: String -> P Bool
optSym s = do
  ts <- lookahead 1
  case ts of
    [Punct p] | p == s -> True <$ skip
    _ -> pure False

sym :: String -> P ()
sym s = do
  found <- optSym s
  unless found (failP ("expected " ++ s))

keywords :: [String]
keywords =
  [ "typedef", "struct", "union", "enum", "const", "volatile", "restrict", "extern", "static"
  , "inline", "register", "auto", "void", "char", "short", "int", "long", "float", "double"
  , "signed", "unsigned", "_Bool"
  ]

identifier :: P String
identifier = do
  ts <- lookahead 1
  case ts of
    [Ident w] | w `notElem` keywords -> w <$ skip
    _ -> failP "expected a name"

end :: P ()
end = do
  ts <- lookahead 1
  unless (null ts) (failP "expected the end of the declaration")

sepBy1 :: P a -> String -> P [a]
sepBy1 p s = do
  x <- p
  more <- optSym s
  if more then (x :) <$> sepBy1 p s else pure [x]

baseWords :: [String]
baseWords = ["void", "char", "short", "int", "long", "float", "double", "signed", "unsigned", "_Bool"]

-- | The types of <stddef.h> and <stdint.h> the headers use.
standardTypes :: [String]
standardTypes =
  [ "size_t", "ptrdiff_t", "intptr_t", "uintptr_t", "int8_t", "int16_t", "int32_t", "int64_t"
  , "uint8_t", "uint16_t", "uint32_t", "uint64_t"
  ]

-- | The specifiers of a declaration: its base type.
specifiers :: P CType
specifiers = go False [] Nothing
  where
    go isConst ws named = do
      ts <- lookahead 2
      case ts of
        (Ident w : rest)
          | w == "const" -> skip >> go True ws named
          | w `elem` ["restrict", "extern"] -> skip >> go isConst ws named
          | w `elem` baseWords -> skip >> go isConst (ws ++ [w]) named
          | w `elem` ["struct", "enum"] -> case rest of
              [Ident tag] | null ws, isNothing named -> skip >> skip >> go isConst ws (Just tag)
              _ -> failP "a structure or an enumeration defined inside a declaration"
          | w `elem` ["union", "volatile", "static", "inline", "register", "auto"] ->
              failP ("not read: " ++ w)
          | null ws, isNothing named -> skip >> go isConst ws (Just w)
        _ -> finish isConst ws named
    finish isConst ws named = do
      base <- case (ws, named) of
        ([], Nothing) -> failP "expected a type"
        ([], Just w)
          | w `elem` standardTypes -> pure (TBase w)
          | "myb_" `isPrefixOf` w -> pure (TName w)
          | otherwise -> failP ("a type that is not of the library: " ++ w)
        (_, Nothing) -> TBase <$> canonical ws
        _ -> failP "two types in one declaration"
      pure (if isConst then TConst base else base)

-- | The C spelling of a base type given by its words.
canonical :: [String] -> P String
canonical ws
  | ws == ["void"] = pure "void"
  | ws == ["_Bool"] = pure "_Bool"
  | ws == ["float"] = pure "float"
  | ws == ["double"] = pure "double"
  | count "char" == 1 && only ["char", "signed", "unsigned"] =
      pure (if unsigned then "unsigned char" else if signed then "signed char" else "char")
  | count "short" == 1 && only ["short", "int", "signed", "unsigned"] =
      pure (if unsigned then "unsigned short" else "short")
  | count "long" == 2 && only ["long", "int", "signed", "unsigned"] =
      pure (if unsigned then "unsigned long long" else "long long")
  | count "long" == 1 && only ["long", "int", "signed", "unsigned"] =
      pure (if unsigned then "unsigned long" else "long")
  | not (null ws) && only ["int", "signed", "unsigned"] = pure (if unsigned then "unsigned int" else "int")
  | otherwise = failP ("a type not read: " ++ unwords ws)
  where
    count w = length (filter (== w) ws)
    only allowed = all (`elem` allowed) ws
    unsigned = count "unsigned" > 0
    signed = count "signed" > 0

-- | After a '*': whether the pointer is const.
qualifiers :: P Bool
qualifiers = go False
  where
    go c = do
      ts <- lookahead 1
      case ts of
        [Ident "const"] -> skip >> go True
        [Ident w] | w `elem` ["restrict", "__restrict", "__restrict__"] -> skip >> go c
        _ -> pure c

pointers :: P [Bool]
pointers = do
  star <- optSym "*"
  if star then (:) <$> qualifiers <*> pointers else pure []

-- | A declarator: its name, if any, and what it makes of the base type.
declarator :: P (Maybe String, CType -> CType)
declarator = do
  stars <- pointers
  (name, core) <- direct
  suffixes <- suffixList
  let withStars t = foldl (\acc c -> (if c then TConst else id) (TPtr acc)) t stars
  pure (name, core . (\t -> foldr ($) t suffixes) . withStars)
  where
    direct = do
      ts <- lookahead 2
      case ts of
        (Ident w : _) | w `notElem` keywords -> skip >> pure (Just w, id)
        [Punct "(", Punct "*"] -> nested
        [Punct "(", Punct "("] -> nested
        _ -> pure (Nothing, id)
    nested = do
      sym "("
      d <- declarator
      sym ")"
      pure d

suffixList :: P [CType -> CType]
suffixList = do
  ts <- lookahead 1
  case ts of
    [Punct "["] -> do
      skip
      size <- untilClose "]"
      rest <- suffixList
      pure ((\t -> TArray t size) : rest)
    [Punct "("] -> do
      skip
      ps <- parameters
      sym ")"
      rest <- suffixList
      pure ((\t -> TFunc t ps) : rest)
    _ -> pure []

-- | The tokens up to the closing bracket, as one text.
untilClose :: String -> P String
untilClose close = go (0 :: Int) []
  where
    go depth acc = do
      t <- next
      case t of
        Punct p
          | p == close && depth == 0 -> pure (unwords (map tokenText (reverse acc)))
          | p `elem` ["(", "["] -> go (depth + 1) (t : acc)
          | p `elem` [")", "]"] -> go (depth - 1) (t : acc)
        _ -> go depth (t : acc)

parameters :: P [Param]
parameters = do
  ts <- lookahead 2
  case ts of
    [Ident "void", Punct ")"] -> [] <$ skip
    (Punct ")" : _) -> pure []
    _ -> sepBy1 parameter ","
  where
    parameter = do
      base <- specifiers
      (name, f) <- declarator
      pure (Param name (adjust (f base)))
    -- an array parameter is a pointer, a function parameter a pointer to it
    adjust t = case t of
      TArray e _ -> TPtr e
      TFunc _ _ -> TPtr t
      _ -> t

structBody :: P [Field]
structBody = sym "{" >> loop []
  where
    loop acc = do
      close <- optSym "}"
      if close
        then pure (concat (reverse acc))
        else do
          base <- specifiers
          ds <- sepBy1 declarator ","
          sym ";"
          fields <- mapM (field base) ds
          loop (fields : acc)
    field base (Just n, f) = pure (Field n (f base))
    field _ (Nothing, _) = failP "a field without a name"

enumBody :: P [String]
enumBody = sym "{" >> loop []
  where
    loop acc = do
      close <- optSym "}"
      if close
        then pure (reverse acc)
        else do
          n <- identifier
          valued <- optSym "="
          when valued (value (0 :: Int))
          _ <- optSym ","
          loop (n : acc)
    -- the value is the header's: hsc2hs reads it there
    value depth = do
      ts <- lookahead 1
      case ts of
        [] -> failP "an enumeration without its end"
        [Punct p]
          | depth == 0 && p `elem` [",", "}"] -> pure ()
          | p `elem` ["(", "["] -> skip >> value (depth + 1)
          | p `elem` [")", "]"] -> skip >> value (depth - 1)
        _ -> skip >> value depth

tagThenBody :: P (Maybe String, Bool)
tagThenBody = do
  ts <- lookahead 2
  tag <- case ts of
    [Ident w, Punct "{"] | w `notElem` keywords -> Just <$> identifier
    [Ident w, Ident _] | w `notElem` keywords -> Just <$> identifier
    _ -> pure Nothing
  rest <- lookahead 1
  pure (tag, rest == [Punct "{"])

declaration :: P Decl
declaration = do
  ts <- lookahead 3
  case ts of
    (Ident "typedef" : Ident "struct" : _) -> do
      skip >> skip
      (tag, hasBody) <- tagThenBody
      body <- if hasBody then Just <$> structBody else pure Nothing
      name <- identifier
      end
      pure (DStruct name tag body)
    (Ident "typedef" : Ident "enum" : _) -> do
      skip >> skip
      (_, hasBody) <- tagThenBody
      unless hasBody (failP "an enumeration declared without its constants")
      names <- enumBody
      name <- identifier
      end
      pure (DEnum (Just name) names)
    (Ident "typedef" : _) -> do
      skip
      base <- specifiers
      (name, f) <- declarator
      end
      maybe (failP "a typedef without a name") (\n -> pure (DTypedef n (f base))) name
    [Ident "struct", Ident tag, Punct "{"] -> do
      skip >> skip
      fields <- structBody
      end
      pure (DTagBody tag fields)
    (Ident "enum" : Punct "{" : _) -> do
      skip
      names <- enumBody
      end
      pure (DEnum Nothing names)
    _ -> do
      base <- specifiers
      (name, f) <- declarator
      end
      case (name, f base) of
        (Just n, TFunc r ps) -> pure (DFunction n r ps)
        _ -> failP "a declaration that is not a function, a type or constants"

-- | The declarations of a text: split at the semicolons outside braces, then read.
declarations :: String -> Either String [Decl]
declarations text = do
  tokens <- tokenize text
  pieces <- split (0 :: Int) [] [] tokens
  forM pieces $ \piece -> case runP declaration piece of
    Right (d, []) -> Right d
    Right (_, rest) -> Left ("left over: " ++ unwords (map tokenText rest))
    Left e -> Left e
  where
    split _ [] done [] = Right (reverse done)
    split _ _ _ [] = Left "the declarations end inside a declaration"
    split depth cur done (t : rest) = case t of
      Punct ";" | depth == 0 -> split 0 [] (reverse cur : done) rest
      Punct "{" -> split (depth + 1) (t : cur) done rest
      Punct "}" -> split (depth - 1) (t : cur) done rest
      _ -> split depth (t : cur) done rest

-- ---- The headers ---------------------------------------------------------------------------

-- | The key of a header of the library from a path of a line marker: "time" for
-- include/myb/time.h, "ext/sgp4" for ext/sgp4/myb_sgp4.h; nothing for any other file.
headerKey :: String -> Maybe String
headerKey path =
  case reverse (filter (not . null) (splitOn (map slash path))) of
    (file : "myb" : "include" : _) | ".h" `isSuffixOf` file -> Just (dropExtension file)
    (file : n : "ext" : _) | file == "myb_" ++ n ++ ".h" -> Just ("ext/" ++ n)
    _ -> Nothing
  where
    slash c = if c == '\\' then '/' else c
    splitOn s = case break (== '/') s of
      (a, []) -> [a]
      (a, _ : rest) -> a : splitOn rest

-- | The text of each header of the library, in the order of their first declarations: a header
-- comes after those it includes.
byHeader :: String -> [(String, String)]
byHeader output = [(k, unlines (reverse (texts Map.! k))) | k <- reverse order]
  where
    (order, texts, _) = foldl step ([], Map.empty, Nothing) (lines output)
    step (o, m, current) line
      | Just path <- marker line = (o, m, headerKey path)
      | all isSpace line = (o, m, current)
      | Just k <- current =
          if Map.member k m then (o, Map.adjust (line :) k m, current) else (k : o, Map.insert k [line] m, current)
      | otherwise = (o, m, current)
    marker line = case line of
      ('#' : ' ' : rest@(d : _)) | isDigit d -> Just (takeWhile (/= '"') (drop 1 (dropWhile (/= '"') rest)))
      _ -> Nothing

preprocess :: FilePath -> FilePath -> [FilePath] -> String -> IO String
preprocess clang root extra source = do
  let arguments =
        ["-E", "-x", "c", "-std=c17", "-ffreestanding", "-nostdlibinc", "-I", root </> "include"]
          ++ concatMap (\d -> ["-I", d]) extra
          ++ ["-DMYB_STATIC", "-DMYB_NODISCARD=", "-"]
  (code, out, err) <- readProcessWithExitCode clang arguments source
  case code of
    ExitSuccess -> pure (filter (/= '\r') out)
    ExitFailure _ -> die (clang ++ " could not read the headers:\n" ++ err)

-- | The object-like #define of a header whose value is an integer constant: its name. A macro
-- of attributes (MYB_API, ...) is not one; a floating value stops the generator, as hsc2hs
-- reads integers only.
defines :: String -> Either String [String]
defines text = fmap (nub . concat) (mapM constant (joined (lines (filter (/= '\r') text))))
  where
    joined (l : rest)
      | "\\" `isSuffixOf` l, (m : rest') <- rest = joined ((init l ++ " " ++ m) : rest')
      | otherwise = l : joined rest
    joined [] = []
    constant line = case words (map (\c -> if c == '#' then ' ' else c) (dropWhile isSpace line)) of
      ("define" : name : body)
        | "#" `isPrefixOf` dropWhile isSpace line
        , "MYB_" `isPrefixOf` name
        , '(' `notElem` name
        , not (null body)
        , not (any (`isInfixOf` unwords body) ["__attribute__", "__declspec", "[["]) ->
            case tokenize (unwords body) of
              Right ts | not (any floating ts) -> Right [name]
              _ -> Left ("a #define that is not an integer constant: " ++ name)
      _ -> Right []
    floating t = case t of
      Number w -> not ("0x" `isPrefixOf` w || "0X" `isPrefixOf` w) && any (`elem` ".eE") w
      _ -> False

-- ---- The model -----------------------------------------------------------------------------

data Kind = KStruct [Field] | KOpaque | KEnum [String] | KFunPtr CType [Param]

data Module = Module
  { mKey :: String
  , mTypes :: [String] -- the types the module defines, in order
  , mConstants :: [String] -- of anonymous enumerations and of #define
  , mFunctions :: [(String, CType, [Param])]
  }

data Model = Model
  { types :: Map.Map String (String, Kind) -- name: the key of its module, what it is
  , modules :: [Module]
  }

data Build = Build
  { bTypes :: Map.Map String (String, Kind)
  , bOrder :: [String]
  , bTags :: Map.Map String String
  , bConstants :: Map.Map String [String]
  , bFunctions :: Map.Map String [(String, CType, [Param])]
  }

model :: [(String, [Decl], [String])] -> Either String Model
model headers = do
  built <- foldM header (Build Map.empty [] Map.empty Map.empty Map.empty) headers
  let known = bTypes built
      order = reverse (nub (bOrder built)) -- each name where it was last placed
      mods =
        [ Module
            key
            [n | n <- order, fmap fst (Map.lookup n known) == Just key]
            (Map.findWithDefault [] key (bConstants built) ++ macros)
            (reverse (Map.findWithDefault [] key (bFunctions built)))
        | (key, _, macros) <- headers
        ]
      unknown = Set.toList (Set.fromList (concatMap (references known) mods) Set.\\ Map.keysSet known)
  unless (null unknown) (Left ("types used but not declared: " ++ unwords unknown))
  pure (Model known [m | m <- mods, not (null (mTypes m) && null (mConstants m) && null (mFunctions m))])
  where
    header b (key, decls, _) = foldM (declare key) b decls
    place key name kind b = b {bTypes = Map.insert name (key, kind) (bTypes b), bOrder = name : bOrder b}
    declare key b d = case d of
      DStruct name tag body -> do
        let tagged = maybe b (\t -> b {bTags = Map.insert t name (bTags b)}) tag
        pure $ case body of
          Just fields -> place key name (KStruct fields) tagged
          Nothing
            | Map.member name (bTypes b) -> tagged
            | otherwise -> place key name KOpaque tagged
      DTagBody tag fields -> case Map.lookup tag (bTags b) of
        Just name -> pure (place key name (KStruct fields) b)
        Nothing -> Left ("struct " ++ tag ++ " has no typedef")
      DEnum (Just name) names -> pure (place key name (KEnum names) b)
      DEnum Nothing names -> pure b {bConstants = Map.insertWith (flip (++)) key names (bConstants b)}
      DTypedef name t -> case strip t of
        TPtr (TFunc r ps) -> pure (place key name (KFunPtr r ps) b)
        _ -> Left ("a typedef the generator does not read: " ++ name)
      DFunction name r ps -> pure b {bFunctions = Map.insertWith (++) key [(name, r, ps)] (bFunctions b)}

-- | The names of the library's types a module mentions.
references :: Map.Map String (String, Kind) -> Module -> [String]
references known m =
  concatMap ofType (concatMap kindTypes (mapMaybe (\n -> snd <$> Map.lookup n known) (mTypes m)))
    ++ concatMap (\(_, r, ps) -> ofType r ++ concatMap (ofType . paramType) ps) (mFunctions m)
  where
    kindTypes k = case k of
      KStruct fields -> map fieldType fields
      KFunPtr r ps -> r : map paramType ps
      _ -> []

ofType :: CType -> [String]
ofType t = case t of
  TBase _ -> []
  TName n -> [n]
  TConst x -> ofType x
  TPtr x -> ofType x
  TArray x _ -> ofType x
  TFunc r ps -> ofType r ++ concatMap (ofType . paramType) ps

-- | The modules a module imports; an order of the modules in which each comes after those it
-- imports, or the cycle that forbids one.
dependencies :: Model -> Module -> [String]
dependencies m md =
  sort (nub [owner | n <- references (types m) md, Just (owner, _) <- [Map.lookup n (types m)], owner /= mKey md])

checkOrder :: Model -> Either String ()
checkOrder m = go Set.empty (map mKey (modules m))
  where
    deps = Map.fromList [(mKey md, dependencies m md) | md <- modules m]
    go _ [] = Right ()
    go done (k : rest) =
      case filter (`Set.notMember` done) (Map.findWithDefault [] k deps) of
        [] -> go (Set.insert k done) rest
        later -> Left (k ++ " needs types of " ++ unwords later ++ ", which come after it")

-- ---- Haskell -------------------------------------------------------------------------------

moduleName :: String -> String
moduleName key = case key of
  'e' : 'x' : 't' : '/' : n -> "Myb.C.Ext." ++ capital n
  _ -> "Myb.C." ++ capital key

capital :: String -> String
capital (c : cs) = toUpper c : cs
capital [] = []

hsName :: String -> String
hsName = capital

hsBase :: String -> String
hsBase b = fromMaybe (error ("no Haskell type for " ++ b)) (lookup b baseTypes)

baseTypes :: [(String, String)]
baseTypes =
  [ ("void", "()"), ("char", "CChar"), ("signed char", "CSChar"), ("unsigned char", "CUChar")
  , ("short", "CShort"), ("unsigned short", "CUShort"), ("int", "CInt"), ("unsigned int", "CUInt")
  , ("long", "CLong"), ("unsigned long", "CULong"), ("long long", "CLLong")
  , ("unsigned long long", "CULLong"), ("float", "CFloat"), ("double", "CDouble")
  , ("_Bool", "CBool"), ("size_t", "CSize"), ("ptrdiff_t", "CPtrdiff"), ("intptr_t", "CIntPtr")
  , ("uintptr_t", "CUIntPtr"), ("int8_t", "Int8"), ("int16_t", "Int16"), ("int32_t", "Int32")
  , ("int64_t", "Int64"), ("uint8_t", "Word8"), ("uint16_t", "Word16"), ("uint32_t", "Word32")
  , ("uint64_t", "Word64")
  ]

atom :: String -> String
atom s = if ' ' `elem` s then "(" ++ s ++ ")" else s

hsType :: CType -> String
hsType t = case t of
  TConst x -> hsType x
  TBase b -> hsBase b
  TName n -> hsName n
  TPtr x -> case strip x of
    TBase "void" -> "Ptr ()"
    TBase "char" -> "CString"
    TFunc r ps -> "FunPtr (" ++ hsFun r ps ++ ")"
    _ -> "Ptr " ++ atom (hsType x)
  TArray x _ -> "[" ++ hsType x ++ "]"
  TFunc r ps -> "(" ++ hsFun r ps ++ ")"

hsFun :: CType -> [Param] -> String
hsFun r ps = intercalate " -> " (map (hsType . paramType) ps ++ ["IO " ++ atom (hsType r)])

-- | The C text of a declaration of the type with the name (an empty name: abstract).
renderC :: CType -> String -> String
renderC t d = case t of
  TBase b -> b ++ spaced d
  TName n -> n ++ spaced d
  TConst (TBase b) -> "const " ++ b ++ spaced d
  TConst (TName n) -> "const " ++ n ++ spaced d
  TConst inner -> renderC inner ("const" ++ spaced d)
  TPtr inner -> renderC inner (around inner ("*" ++ d))
  TArray inner n -> renderC inner (d ++ "[" ++ n ++ "]")
  TFunc r ps -> renderC r (d ++ "(" ++ renderParams ps ++ ")")
  where
    spaced "" = ""
    spaced x = ' ' : x
    around inner x = case strip inner of
      TArray _ _ -> "(" ++ x ++ ")"
      TFunc _ _ -> "(" ++ x ++ ")"
      _ -> x

renderParams :: [Param] -> String
renderParams [] = "void"
renderParams ps = intercalate ", " [renderC (paramType p) (fromMaybe "" (paramName p)) | p <- ps]

-- | Whether a value of the type is a structure, which a foreign import cannot pass.
byValue :: Model -> CType -> Bool
byValue m t = case strip t of
  TName n -> case snd <$> Map.lookup n (types m) of
    Just (KStruct _) -> True
    Just KOpaque -> True
    _ -> False
  _ -> False

argName :: Int -> String
argName i = 'a' : show i

-- | A C function of the binding that passes the structures of a call through pointers.
--   target: what the shim calls (a function, or a function pointer named by its first
--   parameter); prefix: the parameters before those of the call.
data Shim = Shim {shimName :: String, shimPrefix :: [Param], shimCall :: String, shimReturn :: CType, shimParams :: [Param]}

-- | The parameters of a shim and its result: a structure parameter by a const pointer, a
-- structure result by a pointer after them.
shimSignature :: Model -> CType -> [Param] -> (CType, [Param])
shimSignature m r ps =
  ( if byValue m r then TBase "void" else r
  , [ Param (Just (argName i)) (if byValue m (paramType p) then TPtr (TConst (strip (paramType p))) else paramType p)
    | (i, p) <- zip [1 ..] ps
    ]
      ++ [Param (Just "out") (TPtr (strip r)) | byValue m r]
  )

cShim :: Model -> Shim -> [String]
cShim m s =
  [ renderC (TFunc r' (shimPrefix s ++ ps')) (shimName s) ++ " {"
  , body
  , "}"
  , ""
  ]
  where
    (r', ps') = shimSignature m (shimReturn s) (shimParams s)
    call =
      shimCall s ++ "(" ++ intercalate ", " [(if byValue m (paramType p) then "*" else "") ++ argName i | (i, p) <- zip [1 ..] (shimParams s)] ++ ")"
    body
      | byValue m (shimReturn s) = "  *out = " ++ call ++ ";"
      | strip (shimReturn s) == TBase "void" = "  " ++ call ++ ";"
      | otherwise = "  return " ++ call ++ ";"

-- | The Haskell function over a shim: it passes structures through `with` and takes a
-- structure result from `alloca`.
shimWrapper :: Model -> String -> String -> [String] -> CType -> [Param] -> String
shimWrapper m name target leading r ps =
  name ++ concatMap (' ' :) (leading ++ map argName [1 .. length ps]) ++ " = " ++ withs ++ result
  where
    indexed = zip [1 :: Int ..] ps
    withs = concat ["with " ++ argName i ++ " $ \\p" ++ show i ++ " -> " | (i, p) <- indexed, byValue m (paramType p)]
    arg (i, p) = if byValue m (paramType p) then 'p' : show i else argName i
    call extra = unwords ([target] ++ leading ++ map arg indexed ++ extra)
    result
      | byValue m r = "alloca $ \\out -> " ++ call ["out"] ++ " >> peek out"
      | otherwise = call []

data Out = Out {hsLines :: [String], cLines :: [String], notes :: [String]}

instance Semigroup Out where
  Out a b c <> Out a' b' c' = Out (a ++ a') (b ++ b') (c ++ c')

instance Monoid Out where
  mempty = Out [] [] []

emitFunction :: Model -> (String, CType, [Param]) -> Out
emitFunction m (name, r, ps)
  | needsShim =
      Out
        ( doc
            ++ [ "foreign import ccall \"hs_" ++ name ++ "\" hs'" ++ name ++ " :: " ++ hsFun r' ps'
               , name ++ " :: " ++ hsFun r ps
               , shimWrapper m name ("hs'" ++ name) [] r ps
               , address
               ]
        )
        (["/* " ++ renderC (TFunc r ps) name ++ " */"] ++ cShim m (Shim ("hs_" ++ name) [] name r ps))
        []
  | otherwise = Out (doc ++ ["foreign import ccall \"" ++ name ++ "\" " ++ name ++ " :: " ++ hsFun r ps, address]) [] []
  where
    needsShim = byValue m r || any (byValue m . paramType) ps
    (r', ps') = shimSignature m r ps
    doc = ["", "-- | " ++ renderC (TFunc r ps) name]
    address = "foreign import ccall \"&" ++ name ++ "\" p'" ++ name ++ " :: FunPtr (" ++ hsFun r ps ++ ")"

-- | mk'<name> and call'<name> for a function pointer type, whether it is a typedef or written
-- inside a structure or a parameter list.
emitFunPtr :: Model -> String -> String -> CType -> [Param] -> Out
emitFunPtr m name hsTy r ps = emitFunPtrNamed m name (map (\c -> if c == '\'' then '_' else c) name) hsTy r ps

-- | name: the Haskell part of the names; cname: the same as a C identifier.
emitFunPtrNamed :: Model -> String -> String -> String -> CType -> [Param] -> Out
emitFunPtrNamed m name cname hsTy r ps
  | not structs =
      Out
        [ ""
        , "-- | A Haskell function as a " ++ name ++ ". Free it with freeHaskellFunPtr once C"
        , "-- holds it no more."
        , "foreign import ccall \"wrapper\" mk'" ++ name ++ " :: (" ++ sig ++ ") -> IO " ++ atom hsTy
        , ""
        , "-- | Calls a " ++ name ++ "."
        , "foreign import ccall \"dynamic\" call'" ++ name ++ " :: " ++ atom hsTy ++ " -> " ++ sig
        ]
        []
        []
  | otherwise = maybe noContext trampoline context <> dynamic
  where
    sig = hsFun r ps
    structs = byValue m r || any (byValue m . paramType) ps
    (r', ps') = shimSignature m r ps
    pointerSig = hsFun r' ps'
    context = lookup True [(strip (paramType p) == TPtr (TBase "void"), i) | (i, p) <- zip [1 :: Int ..] ps]
    noContext = Out [] [] ["no mk'" ++ name ++ ": it passes a structure by value and has no context"]
    trampoline k =
      Out
        [ ""
        , "foreign import ccall \"wrapper\" hs'mk'" ++ name ++ " :: (" ++ pointerSig ++ ") -> IO (FunPtr (" ++ pointerSig ++ "))"
        , "foreign import ccall \"&hs_trampoline_" ++ cname ++ "\" p'hs_trampoline_" ++ name ++ " :: " ++ hsTy
        , ""
        , "-- | A Haskell function as a " ++ name ++ ". The type passes a structure by value,"
        , "-- which a Haskell function cannot take: the C function returned is a trampoline of"
        , "-- the binding and the context to give it is the Haskell function, which receives"
        , "-- nullPtr for its context. Free it with freeHaskellFunPtr (castPtrToFunPtr context)"
        , "-- once C holds it no more."
        , "mk'" ++ name ++ " :: (" ++ sig ++ ") -> IO (" ++ hsTy ++ ", Ptr ())"
        , "mk'" ++ name ++ " f = do"
        , "  w <- hs'mk'" ++ name ++ " (\\" ++ unwords lambdaArgs ++ " -> " ++ peeks ++ inner ++ ")"
        , "  pure (p'hs_trampoline_" ++ name ++ ", castFunPtrToPtr w)"
        ]
        ( [ "typedef " ++ renderC (TPtr (TFunc r' ps')) ("hs_" ++ cname) ++ ";"
          , renderC (TFunc r [Param (Just (argName i)) (paramType p) | (i, p) <- zip [1 ..] ps]) ("hs_trampoline_" ++ cname) ++ " {"
          ]
            ++ ["  " ++ renderC (strip r) "out" ++ ";" | byValue m r]
            ++ [ "  " ++ (if strip r == TBase "void" || byValue m r then "" else "return ")
                   ++ "((hs_" ++ cname ++ ")" ++ argName k ++ ")(" ++ intercalate ", " (map callArg (zip [1 ..] ps) ++ ["&out" | byValue m r]) ++ ");"
               ]
            ++ ["  return out;" | byValue m r]
            ++ ["}", ""]
        )
        []
      where
        callArg (i, p)
          | i == k = "NULL"
          | byValue m (paramType p) = '&' : argName i
          | otherwise = argName i
        lambdaArgs = map argName [1 .. length ps] ++ ["out" | byValue m r]
        peeks = concat ["peek " ++ argName i ++ " >>= \\v" ++ show i ++ " -> " | (i, p) <- zip [1 :: Int ..] ps, byValue m (paramType p)]
        inner =
          unwords ("f" : [if byValue m (paramType p) then 'v' : show i else argName i | (i, p) <- zip [1 :: Int ..] ps])
            ++ (if byValue m r then " >>= poke out" else "")
    dynamic =
      Out
        [ ""
        , "foreign import ccall \"hs_call_" ++ cname ++ "\" hs'call'" ++ name ++ " :: " ++ atom hsTy ++ " -> " ++ pointerSig
        , ""
        , "-- | Calls a " ++ name ++ "."
        , "call'" ++ name ++ " :: " ++ atom hsTy ++ " -> " ++ sig
        , shimWrapper m ("call'" ++ name) ("hs'call'" ++ name) ["f"] r ps
        ]
        (cShim m (Shim ("hs_call_" ++ cname) [Param (Just "f") (TPtr (TFunc r ps))] "f" r ps))
        []

emitStruct :: Model -> String -> [Field] -> Out
emitStruct m name fields =
  Out
    ( [ ""
      , "-- | struct " ++ name
      , "data " ++ h ++ " = " ++ h
      ]
        ++ zipWith (\i f -> (if i == (1 :: Int) then "  { " else "  , ") ++ selector f ++ " :: " ++ fieldHs f) [1 ..] fields
        ++ [ "  }"
           , "  deriving (Eq, Show)"
           , ""
           , "instance Storable " ++ h ++ " where"
           , "  sizeOf _ = #{size " ++ name ++ "}"
           , "  alignment _ = #{alignment " ++ name ++ "}"
           , "  peek p = do"
           ]
        ++ zipWith peekField [1 :: Int ..] fields
        ++ [ "    pure (" ++ unwords (h : map (\i -> 'v' : show i) [1 .. length fields]) ++ ")"
           , "  poke p v = do"
           , "    fillBytes p 0 #{size " ++ name ++ "}"
           ]
        ++ map pokeField fields
        ++ [ ""
           , "instance Zero " ++ h ++ " where"
           , "  zero = " ++ unwords (h : map zeroField fields)
           ]
    )
    []
    []
    <> mconcat [emitFunPtr m (name ++ "'" ++ fieldName f) (hsType (fieldType f)) r ps | f <- fields, TPtr (TFunc r ps) <- [strip (fieldType f)]]
  where
    h = hsName name
    selector f = name ++ "'" ++ fieldName f
    flat t = case t of
      TArray e n -> let (base, ns) = flat e in (base, n : ns)
      _ -> (t, [])
    count ns = "#{const " ++ intercalate " * " ["(" ++ n ++ ")" | n <- ns] ++ "}"
    fieldHs f = case flat (fieldType f) of
      (_, []) -> hsType (fieldType f)
      (base, _) -> "[" ++ hsType base ++ "]"
    peekField i f = case flat (fieldType f) of
      (_, []) -> "    v" ++ show i ++ " <- #{peek " ++ name ++ ", " ++ fieldName f ++ "} p"
      (_, ns) -> "    v" ++ show i ++ " <- peekArray (" ++ count ns ++ ") (#{ptr " ++ name ++ ", " ++ fieldName f ++ "} p)"
    pokeField f = case flat (fieldType f) of
      (_, []) -> "    #{poke " ++ name ++ ", " ++ fieldName f ++ "} p (" ++ selector f ++ " v)"
      (_, ns) -> "    pokeArray (#{ptr " ++ name ++ ", " ++ fieldName f ++ "} p) (take (" ++ count ns ++ ") (" ++ selector f ++ " v))"
    zeroField f = case flat (fieldType f) of
      (_, []) -> "zero"
      (_, ns) -> "(replicate (" ++ count ns ++ ") zero)"

emitType :: Model -> String -> Out
emitType m name = case snd (types m Map.! name) of
  KStruct fields -> emitStruct m name fields
  KOpaque -> Out ["", "-- | struct " ++ name ++ ": an object of the library, known by pointer.", "data " ++ hsName name] [] []
  KEnum names ->
    Out
      ( [ ""
        , "-- | enum " ++ name
        , "newtype " ++ hsName name ++ " = " ++ hsName name ++ " #{type " ++ name ++ "}"
        , "  deriving stock (Show)"
        , "  deriving newtype (Eq, Ord, Storable)"
        , ""
        , "instance Zero " ++ hsName name ++ " where"
        , "  zero = " ++ hsName name ++ " 0"
        ]
          ++ concat [["", "pattern " ++ c ++ " :: " ++ hsName name, "pattern " ++ c ++ " = " ++ hsName name ++ " (#{const " ++ c ++ "})"] | c <- names]
      )
      []
      []
  KFunPtr r ps ->
    Out
      ["", "-- | " ++ renderC (TPtr (TFunc r ps)) name, "type " ++ hsName name ++ " = FunPtr (" ++ hsFun r ps ++ ")"]
      []
      []
      <> emitFunPtr m name (hsName name) r ps

emitConstant :: String -> Out
emitConstant c =
  Out ["", "pattern " ++ c ++ " :: (Eq a, Num a) => a", "pattern " ++ c ++ " = (#{const " ++ c ++ "})"] [] []

-- | The function pointers written in the parameter lists of the functions.
emitParameterFunPtrs :: Model -> (String, CType, [Param]) -> Out
emitParameterFunPtrs m (name, _, ps) =
  mconcat
    [ emitFunPtr m (name ++ "'" ++ fromMaybe (show i) (paramName p)) (hsType (paramType p)) r qs
    | (i, p) <- zip [1 :: Int ..] ps
    , TPtr (TFunc r qs) <- [strip (paramType p)]
    ]

banner :: String -> String -> [String]
banner version what =
  [ "-- Generated by bindings/haskell/tools/Generate.hs (myb-generate) from the public headers of"
  , "-- myb " ++ version ++ ": " ++ what ++ ". Do not edit; bindings/haskell/README.md"
  , "-- says how to write it again."
  ]

headerFile :: String -> String
headerFile key = case key of
  'e' : 'x' : 't' : '/' : n -> "myb_" ++ n ++ ".h"
  _ -> "myb/" ++ key ++ ".h"

hscModule :: String -> Model -> Module -> (String, Out)
hscModule version m md = (text, body)
  where
    body =
      mconcat (map (emitType m) (mTypes md))
        <> mconcat (map emitConstant (mConstants md))
        <> mconcat (map (emitFunction m) (mFunctions md))
        <> mconcat (map (emitParameterFunPtrs m) (mFunctions md))
    text =
      unlines $
        banner version ("the declarations of " ++ headerFile (mKey md))
          ++ [ "{-# LANGUAGE DerivingStrategies #-}"
             , "{-# LANGUAGE GeneralizedNewtypeDeriving #-}"
             , "{-# LANGUAGE PatternSynonyms #-}"
             , "{-# OPTIONS_GHC -Wno-unused-imports #-}"
             , ""
             , "-- | The declarations of " ++ headerFile (mKey md) ++ ", by their C names (see \"Myb.C\")."
             , "module " ++ moduleName (mKey md) ++ " where"
             , ""
             , "#include \"myb/myb.h\""
             ]
          ++ ["#include \"" ++ headerFile (mKey md) ++ "\"" | "ext/" `isPrefixOf` mKey md]
          ++ [ ""
             , "import Data.Int"
             , "import Data.Word"
             , "import Foreign.C.String (CString)"
             , "import Foreign.C.Types"
             , "import Foreign.Marshal.Alloc (alloca)"
             , "import Foreign.Marshal.Array (peekArray, pokeArray)"
             , "import Foreign.Marshal.Utils (fillBytes, with)"
             , "import Foreign.Ptr"
             , "import Foreign.Storable"
             , "import Myb.C.Zero"
             ]
          ++ map (\k -> "import " ++ moduleName k) (dependencies m md)
          ++ hsLines body

zeroModule :: String -> String
zeroModule version =
  unlines $
    banner version "the zero values of the C types"
      ++ [ ""
         , "-- | Values of all bits zero, what C writes as @= {0}@. Every structure of the C interface"
         , "-- has one, so that a structure is written as the record update of the fields that matter:"
         , "-- @zero { myb_point_mass_term'gm = 3.986004418e14 }@."
         , "module Myb.C.Zero (Zero (..)) where"
         , ""
         , "import Data.Int"
         , "import Data.Word"
         , "import Foreign.C.Types"
         , "import Foreign.Ptr"
         , ""
         , "class Zero a where"
         , "  zero :: a"
         , ""
         ]
      ++ ["instance Zero " ++ hs ++ " where zero = 0" | (_, hs) <- baseTypes, hs /= "()"]
      ++ [ "instance Zero (Ptr a) where zero = nullPtr"
         , "instance Zero (FunPtr a) where zero = nullFunPtr"
         ]

umbrella :: String -> [String] -> String
umbrella version names =
  unlines $
    banner version "the modules of the core library"
      ++ [ ""
         , "-- | The whole C interface of myb, every function, structure, enumeration and constant of"
         , "-- include/myb, by the names of the headers: a function by its C name, its address as"
         , "-- p'<name>; a structure @myb_x@ as the record @Myb_x@ with the fields @myb_x'<field>@, a"
         , "-- Storable instance and 'zero'; a constant as a pattern synonym of its C name; a type of"
         , "-- function pointer with mk'<name> (a Haskell function as one) and call'<name>. Functions"
         , "-- that take or return a structure by value take and return the Haskell value. Everything"
         , "-- is in IO, as C is: the upper layer, \"Myb\", gives ownership and outcomes. The modules of"
         , "-- ext/ are in the libraries of the package of their names (myb:sgp4, ...)."
         , "module Myb.C"
         ]
      ++ zipWith (\i n -> (if i == (1 :: Int) then "  ( module " else "  , module ") ++ n) [1 ..] names
      ++ ["  ) where", ""]
      ++ map ("import " ++) names

cFile :: String -> String -> [String] -> [String] -> String
cFile version what includes body =
  unlines $
    [ "/* Generated by bindings/haskell/tools/Generate.hs (myb-generate) from the public headers of"
    , " * myb " ++ version ++ ": " ++ what ++ ". Do not edit."
    , " *"
    , " * A foreign import of Haskell passes no structure by value: each function here takes the"
    , " * structures of a call by pointer and calls the library; a trampoline gives a structure to"
    , " * a Haskell function by pointer. Nothing is computed, a structure is only copied."
    , " * MYB_STATIC: the declarations without dllimport; a call reaches a shared library through"
    , " * its import library all the same. */"
    , "#ifndef MYB_STATIC"
    , "#define MYB_STATIC"
    , "#endif"
    , "#include <stddef.h>"
    ]
      ++ map (\i -> "#include \"" ++ i ++ "\"") includes
      ++ [""]
      ++ body

-- | The lines of the package description between the marks of a block.
replaceBlock :: String -> [String] -> String -> Either String String
replaceBlock mark block text =
  case break (== begin) ls of
    (before, _ : rest) -> case break (== finish) rest of
      (_, _ : after) -> Right (unlines (before ++ [begin] ++ block ++ [finish] ++ after))
      _ -> Left ("myb.cabal: no line " ++ finish)
    _ -> Left ("myb.cabal: no line " ++ begin)
  where
    ls = lines (filter (/= '\r') text)
    begin = trimmedMark "BEGIN"
    finish = trimmedMark "END"
    trimmedMark w = indent ++ "-- " ++ w ++ " myb-generate: " ++ mark
    indent = if mark == "modules" then "  " else ""

-- ---- Main ----------------------------------------------------------------------------------

readVersion :: FilePath -> IO String
readVersion root = do
  core <- readText (root </> "include" </> "myb" </> "core.h")
  let part name = case [w | l <- lines core, ["#define", n, w] <- [words l], n == "MYB_VERSION_" ++ name] of
        (w : _) -> pure w
        [] -> die ("core.h: no MYB_VERSION_" ++ name)
  parts <- mapM part ["MAJOR", "MINOR", "PATCH", "TWEAK"]
  pure (intercalate "." parts)

readText :: FilePath -> IO String
readText path = withBinaryFile path ReadMode $ \h -> do
  s <- hGetContents h
  length s `seq` pure (filter (/= '\r') s)

writeText :: FilePath -> String -> IO ()
writeText path text = do
  createDirectoryIfMissing True (takeDirectory path)
  withBinaryFile path WriteMode (`hPutStr` text)

-- | The files under a directory, by paths relative to it with '/' on every system.
filesUnder :: FilePath -> IO [FilePath]
filesUnder dir = do
  exists <- doesDirectoryExist dir
  if not exists
    then pure []
    else do
      names <- listDirectory dir
      fmap concat $ forM (sort names) $ \n -> do
        isDir <- doesDirectoryExist (dir </> n)
        if isDir then map ((n ++ "/") ++) <$> filesUnder (dir </> n) else pure [n]

-- | Removes the directories left empty under a directory; whether it is empty itself.
pruneEmpty :: FilePath -> IO Bool
pruneEmpty dir = do
  names <- listDirectory dir
  kept <- fmap concat $ forM names $ \n -> do
    isDir <- doesDirectoryExist (dir </> n)
    if not isDir
      then pure [n]
      else do
        empty <- pruneEmpty (dir </> n)
        if empty then [] <$ removeDirectory (dir </> n) else pure [n]
  pure (null kept)

readHeaders :: FilePath -> FilePath -> String -> [FilePath] -> String -> IO [(String, [Decl], [String])]
readHeaders clang root source extra wanted = do
  out <- preprocess clang root extra source
  forM [(k, t) | (k, t) <- byHeader out, wanted == "" || k == wanted] $ \(key, text) -> do
    decls <- either (\e -> die (headerFile key ++ ": " ++ e)) pure (declarations text)
    let path = case key of
          'e' : 'x' : 't' : '/' : n -> root </> "ext" </> n </> ("myb_" ++ n ++ ".h")
          _ -> root </> "include" </> "myb" </> (key ++ ".h")
    macros <- either (\e -> die (headerFile key ++ ": " ++ e)) pure . defines =<< readText path
    pure (key, decls, macros)

main :: IO ()
main = do
  args <- getArgs
  fallback <- fromMaybe "clang" <$> lookupEnv "MYB_CLANG"
  let check = "--check" `elem` args
      parse c ps rest = case rest of
        ("--clang" : value : more) -> parse value ps more
        ("--check" : more) -> parse c ps more
        (a : more) -> parse c (ps ++ [a]) more
        [] -> (c, ps)
      (clang, positional) = parse fallback [] args
  root <- case positional of
    [r] -> pure r
    _ -> die "usage: myb-generate [--check] [--clang <clang>] <repository root>"
  version <- readVersion root
  core <- readHeaders clang root "#include \"myb/myb.h\"\n" [] ""
  extNames <- sort <$> listDirectory (root </> "ext")
  exts <- fmap concat $ forM extNames $ \n -> do
    let dir = root </> "ext" </> n
    present <- doesFileExist (dir </> ("myb_" ++ n ++ ".h"))
    if present
      then readHeaders clang root ("#include \"myb/myb.h\"\n#include \"myb_" ++ n ++ ".h\"\n") [dir] ("ext/" ++ n)
      else pure []
  m <- either (die . ("the headers: " ++)) pure (model (core ++ exts))
  either (die . ("the modules: " ++)) pure (checkOrder m)
  let (coreModules, extModules) = span (not . ("ext/" `isPrefixOf`) . mKey) (modules m)
      rendered = [(md, hscModule version m md) | md <- modules m]
      -- the modules of the core library and those of ext/ in separate source directories, so
      -- that a library of ext/ takes the core modules from the library myb, never from source
      sourceDir md = if "ext/" `isPrefixOf` mKey md then "generated/ext/" else "generated/core/"
      hscPath md = sourceDir md ++ map (\c -> if c == '.' then '/' else c) (moduleName (mKey md)) ++ ".hsc"
      coreC = concat [cLines out | (md, (_, out)) <- rendered, not ("ext/" `isPrefixOf` mKey md)]
      extC md = concat [cLines out | (md', (_, out)) <- rendered, mKey md' == mKey md]
      extName md = drop 4 (mKey md)
      cPath n = "generated/cbits/" ++ n ++ ".c"
      generated =
        [ ("generated/core/Myb/C/Zero.hs", zeroModule version)
        , ("generated/core/Myb/C.hs", umbrella version ("Myb.C.Zero" : map (moduleName . mKey) coreModules))
        ]
          ++ [(hscPath md, text) | (md, (text, _)) <- rendered]
          ++ [(cPath "myb_hs", cFile version "what the core library needs" ["myb/myb.h"] coreC) | not (null coreC)]
          ++ [ (cPath ("myb_hs_" ++ extName md), cFile version ("what " ++ headerFile (mKey md) ++ " needs") ["myb/myb.h", headerFile (mKey md)] (extC md))
             | md <- extModules
             , not (null (extC md))
             ]
      moduleBlock =
        ["  exposed-modules:"]
          ++ map ("    " ++) ("Myb.C" : "Myb.C.Zero" : map (moduleName . mKey) coreModules)
          ++ ["  c-sources: " ++ cPath "myb_hs" | not (null coreC)]
      extBlock =
        concat
          [ [ "library " ++ extName md
            , "  import: lower"
            , "  visibility: public"
            , "  hs-source-dirs: generated/ext"
            , "  exposed-modules: " ++ moduleName (mKey md)
            ]
              ++ ["  c-sources: " ++ cPath ("myb_hs_" ++ extName md) | not (null (extC md))]
              ++ [ "  include-dirs: ../../ext/" ++ extName md
                 , "  extra-libraries: myb_" ++ extName md
                 , "  build-depends: myb"
                 , ""
                 ]
          | md <- extModules
          ]
  let package = root </> "bindings" </> "haskell"
  cabal <- readText (package </> "myb.cabal")
  cabal' <- either die pure (replaceBlock "modules" moduleBlock cabal >>= replaceBlock "ext" (init' extBlock))
  let files = generated ++ [("myb.cabal", cabal')]
  present <- map ("generated/" ++) <$> filesUnder (package </> "generated")
  let wanted = Set.fromList (map fst files)
      extra = filter (`Set.notMember` wanted) present
  forM_ (concatMap (\(_, (_, out)) -> notes out) rendered) (hPutStrLn stderr . ("note: " ++))
  if check
    then do
      stale <- fmap concat $ forM files $ \(path, text) -> do
        exists <- doesFileExist (package </> path)
        current <- if exists then readText (package </> path) else pure ""
        pure [path | not exists || current /= text]
      let wrong = stale ++ extra
      unless (null wrong) $ do
        hPutStrLn stderr ("out of date: " ++ unwords wrong ++ "; run myb-generate (bindings/haskell/README.md)")
        exitWith (ExitFailure 1)
      putStrLn ("the lower layer agrees with the headers: " ++ show (length files) ++ " files")
    else do
      forM_ files $ \(path, text) -> writeText (package </> path) text
      forM_ extra $ \path -> removeFile (package </> path)
      _ <- pruneEmpty (package </> "generated")
      let count what = length [() | (_, (n, k)) <- Map.toList (types m), what n k]
      putStrLn
        ( "written " ++ show (length files) ++ " files: "
            ++ show (length (modules m)) ++ " modules, "
            ++ show (count (\_ k -> case k of KStruct _ -> True; _ -> False)) ++ " structures, "
            ++ show (count (\_ k -> case k of KEnum _ -> True; _ -> False)) ++ " enumerations, "
            ++ show (count (\_ k -> case k of KFunPtr _ _ -> True; _ -> False)) ++ " function pointer types, "
            ++ show (sum (map (length . mFunctions) (modules m))) ++ " functions"
        )
  where
    init' xs = if not (null xs) && last xs == "" then init xs else xs
