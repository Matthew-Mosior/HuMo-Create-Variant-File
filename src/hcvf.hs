{-=HuMoCreateVariantFile (HCVF): A Haskell-based solution to=-}
{-=creating a master variant file from bam-readcount results=-}
{-=for a given sample.=-}
{-=Author: Matthew Mosior=-}
{-=Version: 1.0=-}
{-=Synopsis:  This Haskell Script will take in=-}
{-=a .txt file describing the transposed variant=-}
{-=files and a empty sample header file describing=-}
{-=all possible sub-samples in a given sample.=-}
{-=Will return a completed sample header file=-}
{-=filled in with the appropriate variants.=-}


{-Imports-}

import Data.ByteString as DB
import Data.ByteString.Char8 as DBC
import Data.Functor as DF
import Data.List as DL
import Data.List.Split as DLS
import Data.Ord as DO
import Data.Traversable as DT
import System.Console.GetOpt as SCG
import System.Process as SP
import System.Environment as SE
import System.Exit as SX
import System.IO as SIO
import System.IO.Temp as SIOT
import Text.PrettyPrint.Boxes as TPB

{---------}


{-Custom CML Option Datatype.-}

data Flag 
    = Verbose               -- -v
    | Version               -- -V -? 
    | OutputFile String     -- -o
    | Help                  -- --help
    deriving (Eq,Ord,Show) 

{-----------------------------}


{-Custom bool functions for Flag Datatype.-}

--isOutputFile -> This function will
--test for OutputFile flag.
isOutputFile :: Flag -> Bool
isOutputFile (OutputFile _) = True
isOutputFile _              = False

{------------------------------------------}


{-Custom extraction functions for Flag Datatype.-}

--extractOutputFile -> This function will
--extract the string associated with 
--OutputFile.
extractOutputFile :: Flag -> String
extractOutputFile (OutputFile x) = x

{------------------------------------------------}


{-Option Description function relating to datatype above.-}

--options -> This function will
--describe flags.
options :: [OptDescr Flag]
options =
    [ Option ['v']     ["verbose"]             (NoArg Verbose)                "Output on stderr.",
      Option ['V','?'] ["version"]             (NoArg Version)                "Show version number.",
      Option ['o']     ["outputfile"]          (ReqArg OutputFile "OUTFILE")  "The output file.",
      Option []        ["help"]                (NoArg Help)                   "Print this help message."
    ] 

--compilerOpts -> This function will
--parse incoming command line arguments.
compilerOpts :: [String] -> IO ([Flag],[String])
compilerOpts argv =
    case getOpt Permute options argv of
        (args,files,[]) ->
            if DL.elem Help args
                then do SIO.hPutStrLn stderr (greeting ++ SCG.usageInfo header options)
                        SX.exitWith SX.ExitSuccess
                else if DL.elem Version args
                    then do SIO.hPutStrLn stderr (greeting ++ version ++ SCG.usageInfo header options)
                            SX.exitWith SX.ExitSuccess
                    else if (DL.length files > 2 || DL.length files < 2) 
                        then do SIO.hPutStrLn stderr (flerror ++ github ++ SCG.usageInfo header options)
                                SX.exitWith (SX.ExitFailure 1)
                        else return (DL.nub args, files) 
        (_,_,errors) -> do
            SIO.hPutStrLn stderr (DL.concat errors ++ SCG.usageInfo header options)
            SX.exitWith (SX.ExitFailure 1)
        where 
            greeting    = "HuMo Create Variant File, Copyright (c) 2019 Matthew Mosior.\n"
            header      = "Usage: hcvf [-vV?o] [Header File] [Data File]"
            version     = "HuMo Create Variant File (HCVF), Version 1.0.\n"
            github      = "Please see https://github.com/Matthew-Mosior/HuMo-Create-Variant-File/wiki for more information.\n" 
            flerror     = "Incorrect number of input files:  Please provide exactly two input files.\nFirst input file  -> Header file\nSecond input file -> Data file\n"

{---------------------------------------------------------}


{-General Utility Functions.-}

--readFileStrict -> This function will
--read the file in as 
readFileStrict :: FilePath -> IO String
readFileStrict = fmap DBC.unpack . DB.readFile

--lineFeed -> This function will
--read the file in and split on
--whitespace, returning a list
--of lists.
lineFeed :: String -> [[String]]
lineFeed [] = []
lineFeed xs = DL.map DL.words (DL.lines xs)

--mapNotLast -> This function will
--work like the traditional map 
--function in Data.List, but not
--map to the last element of a list.
mapNotLast :: (a -> a) -> [a] -> [a]
mapNotLast fn []     = []
mapNotLast fn [x]    = [x]
mapNotLast fn (x:xs) = fn x : mapNotLast fn xs

--lengthCalculator -> This function will
--determine the length to take and drop 
--within the processArgsandFiles function.
lengthCalculator :: [[String]] -> [Int]
lengthCalculator []     = []
lengthCalculator (x:xs) = [DL.length x - 5] ++ (lengthCalculator xs)

--taker -> This function will take
--n elements from the head of a list.
taker :: [[String]] -> [Int] -> [[String]]
taker []     []     = []
taker (x:xs) (y:ys) = [DL.take y x] ++ (taker xs ys) 

--dropper -> This function will drop
--n elements from the head of a list.
dropper :: [[String]] -> [Int] -> [[String]]
dropper []     []     = []
dropper (x:xs) (y:ys) = [DL.drop y x] ++ (dropper xs ys)

--customChunker -> This function will
--chunk each files contents based on the length
--of each files corresponding list.
customChunker :: [[String]] -> [[[String]]]
customChunker []     = []
customChunker (x:xs) = [DLS.chunksOf (DL.length x `div` 4) x] ++ (customChunker xs)   

--variantAnnotator -> This function will
--determine how to appropriately annotate
--the current variants samples.
variantAnnotator :: [[[String]]] -> [String] -> [[String]] -> [[(String,String)]]
variantAnnotator []     []     []       = []
variantAnnotator []     (_:_)  []       = []
variantAnnotator []     []     (_:_)    = []
variantAnnotator []     (_:_)  (_:_)    = []
variantAnnotator (_:_)  _      []       = []
variantAnnotator (x:xs) ys     (r:rs)   = [(DL.zip (DL.replicate (DL.length r) "N/A") r) ++ [ a | b <- ys , a <- unorderedvariants , fst a == b ]] ++ (variantAnnotator xs ys rs)
    where 
        --Prepare list compreshension.
        unorderedvariants = ((DL.map (\z -> (DL.head z,DL.intercalate "," (DL.init (DL.reverse z)))) x) 
                           ++ (DL.zip (ys \\ (DL.map (DL.head) x)) (DL.replicate (DL.length (ys \\ (DL.map (DL.head) x))) "N/A")))
        ------------------------------

{----------------------------}


{-Printing functions.-}

--tempFileCreation -> This function will
--print the file to stdout using
--readProcess of the unix tool cat.
catFile :: [[String]] -> IO ()
catFile [] = return ()
catFile xs = do
    --Open a temporary file.
    (tempfile,temph) <- SIOT.openTempFile "." "temp.txt"
    --Intercalate a tab, and then a newline into xs.
    let intercalatedxs = DL.intercalate "\n" (DL.map (DL.intercalate "\t") xs)
    --Add intercalatedxs to temp.txt.
    SIO.hPutStrLn temph intercalatedxs
    --Close the temporary file's handle.
    hClose temph
    --Print out the contents of tempfile to the screen using cat unix tool.
    (_,_,_,ph) <- SP.createProcess (SP.proc "cat" [tempfile])
    ec <- SP.waitForProcess ph
    case ec of
        SX.ExitSuccess   -> do _ <- SP.readProcess "rm" [tempfile] []
                               return ()
        SX.ExitFailure _ -> do _ <- error "Could not cat file."
                               _ <- SP.readProcess "rm" [tempfile] []
                               return ()

--printFile -> This function will
--print the file to either stdout
--or to a output file based on
--command-lines options provided.
printFile :: [Flag] -> [[String]] -> IO ()
printFile [] [] = return ()
printFile [] _  = return ()
printFile _  [] = return ()
printFile opts xs = do
    --Grab just "OUTFILE".
    let outfile = DL.head (DL.filter (isOutputFile) opts)
    --Extract the string from FilterFields.
    let outfilestring = extractOutputFile outfile
    --mapNotLast tabs and newlines in xs.
    --let tabsandnewlinesadded = DL.map (\y -> [y]) (mapNotLast (++ "\n") (DL.map (DL.concat) (DL.map (mapNotLast (++ "\t")) xs))) 
    let tabsandnewlinesadded = DL.map (mapNotLast (++ "\t")) xs
    --Write the output to the user-specified filename.
    SIO.writeFile (outfilestring) $
                  (TPB.render $
                  (TPB.hsep 0 TPB.left . DL.map (TPB.vcat TPB.left) . DL.map (DL.map (TPB.text)))
                  (DL.transpose tabsandnewlinesadded)) 

{---------------------}

{-Function to load all files listed in data file.-}

--loadAllDataFiles -> This function will 
--load all files listed in the data file.
loadAllDataFiles :: [String] -> IO [String]
loadAllDataFiles [] = return []
--loadAllDataFiles xs = DT.sequence (DL.map (SIO.readFile) xs) 
loadAllDataFiles xs = DT.sequence (DL.map (readFileStrict) xs)

{-------------------------------------------------}

{-HCVF Specific Function.-}

--processArgsAndFiles -> This function will
--walk through each of the command-line
--arguments and files provided by the user.
processArgsAndFiles :: ([Flag],[String]) -> IO ()
processArgsAndFiles ([],[]) = return () 
processArgsAndFiles (options,files) = do
    --Process the data file first.
    --Read in the data file.
    readdatafile <- SIO.readFile (DL.last files)
    --Apply lineFeed function to readdatafile.
    let linesdatafile = DL.lines readdatafile
    --Read all lines of linesdatafile.
    alldata <- loadAllDataFiles linesdatafile
    --Split alldata on whitespaces.
    let splitalldata = DL.map (DL.words) alldata
    --Create sublists (4 items each) on splitalldata.
    let sublistalldata = taker splitalldata (lengthCalculator splitalldata)
    --Aggregate the lines of sublistalldata based on index.
    let aggregatedalldata = DL.map (DL.transpose) (customChunker sublistalldata)
    --Grab the chr,start,stop,ref,alt from linefeedcurrentf.
    let chrstartstoprefaltalldata = dropper splitalldata (lengthCalculator splitalldata)
    ------------------------------
    --Process header file now.
    readheaderfile <- SIO.readFile (DL.head files)
    --Turn readheaderfile into [[String]].
    let linesheaderfile = DL.lines readheaderfile
    -------------------------- 
    --Line up the current samples and annotate aggregatedf.
    let annotatedalldata = variantAnnotator aggregatedalldata linesheaderfile chrstartstoprefaltalldata
    --Grab just snd for all variants across all data.
    let sndannotatedalldata = DL.map (DL.map snd) annotatedalldata
    --Add header to sndannotatedalldata.
    let finalizedalldata = [["Chromosome","Start","Stop","Ref","Alt"] ++ linesheaderfile] ++ sndannotatedalldata    
    --Print the file to stdout (cat) or to a file.
    if DL.length (DL.filter (isOutputFile) options) > 0 
        then printFile options finalizedalldata
        else catFile finalizedalldata

{-------------------------}


{-Main function.-}

main :: IO ()
main = do
    --Get command line arguments.
    (args,files) <- SE.getArgs >>= compilerOpts
    --Run args and files through processArgsandFiles.
    processArgsAndFiles (args,files)
    
{----------------}
