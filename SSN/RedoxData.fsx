#r @".\bin\Debug\netstandard.dll"

#r @".\bin\Debug\BioFSharp.dll"
#r @".\bin\Debug\BioFSharp.IO.dll"
#r @".\bin\Debug\FSharpAux.dll"
#r @".\bin\Debug\FSharpAux.IO.dll"
#r @".\bin\Debug\Newtonsoft.Json.dll"
#r @".\bin\Debug\FSharp.Stats.dll"
#r @".\bin\Debug\FSharp.Plotly.dll"
#r @".\bin\Debug\FSharpGephiStreamer.dll"

#r @"..\lib\MathNet.Numerics.dll"

#time

open System
open System.IO
open FSharp.Stats
open FSharp.Plotly
open BioFSharp.IO

open FSharpAux
open FSharpAux.IO


open BioFSharp
open BioFSharp
open BioFSharp.Elements

let findMatches chars str =
    /// Returns whether or not the string matches the beginning of the character array
    let rec isStartMatch chars (str: string) =
        match chars with
        | char :: rest when str.Length > 0 ->
            char = str.[0] && (isStartMatch rest str.[1..(str.Length - 1)])
        | _ -> str.Length = 0
    /// The actual function here
    let rec findMatches matchedIndices i chars str =
        match chars with
        | _ :: rest ->
            if isStartMatch chars str
            then findMatches (i :: matchedIndices) (i + 1) rest str
            else findMatches matchedIndices (i + 1) rest str
        | [] -> matchedIndices

    findMatches [] 0 chars str

let adjustCysPosition protSeq cysSeq cysPs a =
    if protSeq |> String.contains cysSeq then
        let c = cysPs |> Array.findIndex (fun i -> i=a)
        let pInSub = cysSeq |> String.toCharArray |> Array.indexed |> Array.filter (fun (i,c) -> c='C') |> Array.item c |> fst
        let startSubseq = findMatches (protSeq |> String.toCharArray |> Array.toList) cysSeq
        startSubseq |> List.map (fun x -> x + pInSub) |> List.minBy (fun t -> abs (t-a))
    else a

//Read file of RedoxDB.A

type CysItem =
    {
    CysID: int;
    Organism: string;
    ProtID: string;
    ProteinSeq: string;
    CysPosition: int;
    CysSeq: string;
    Label: bool;
    }

let fillCysItem cysID org protID protSeq cysP cysSeq label = {CysID=cysID; Organism=org; ProtID=protID; ProteinSeq=protSeq; CysPosition=cysP; CysSeq=cysSeq; Label=label}

let enum id item = {item with CysID=id}

type FastaItem = {
    Header    : string;
    Sequence  : string;       
    }

/// Creates with header line and sequence.
let createFastaItem header sequence =
    { Header = header; Sequence = sequence }

// Conditon of grouping lines
let same_group l =             
    not (String.length l = 0 || l.[0] <> '>')
    
/// Reads FastaItem from file. Converter determines type of sequence by converting seq<char> -> type
let fromFileEnumerator (converter:seq<char>-> string) (fileEnumerator) =
    // Matches grouped lines and concatenates them
    let record d (converter:seq<char>-> string) = 
        match d with
        | [] -> raise (System.Exception "Incorrect FASTA format")
        | (h:string) :: t when h.StartsWith ">" ->  let header = h .Remove(0,1)
                                                    let sequence = (Seq.concat t) |> converter
                                                    createFastaItem header sequence
                                                        
        | h :: _ -> raise (System.Exception "Incorrect FASTA format")        
    
    // main
    fileEnumerator
    |> Seq.filter (fun (l:string) -> not (l.StartsWith ";" || l.StartsWith "#"))
    |> Seq.groupWhen same_group 
    |> Seq.map (fun l -> record (List.ofSeq l) converter)

/// Reads FastaItem from file. Converter determines type of sequence by converting seq<char> -> type
let fromFile converter (filePath) =
    FileIO.readFile filePath
    |> fromFileEnumerator converter


////////////// read RedoxData

let linesA = FileIO.readFile @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\redoxdbA.txt"

let itemsA =
    linesA 
    |> Seq.fold (fun (item,list) i -> 
        if i="" then ([],item::list)
        else (i::item,list)) ([],[])
    |> snd

let readA (line: string list) =
    let org = line |> List.find (fun s -> s |> String.startsWith "ORGANISM") |> String.split ':' |> Array.item 1
    let protID = line.[1] |>  String.replace ">" "" 
    let protSeq = if (line.[0] |> String.startsWith protID) then (line.[0].[protID.Length ..]) else line.[0]
    let cysP = line |> List.filter (fun s -> s |> String.startsWith "CYSTEINE") |> List.map (String.split ' ' >> Array.item 1 >> String.replace "," "" >> String.tryParseIntDefault 0)
    let label = true
    [for a in cysP -> 
        let cP = 
            if protSeq.[a-1]='C' then a-1
            elif protSeq.[a]='C' then a
            elif protSeq.[a+1]='C' then a+1
            elif protSeq.[a+2]='C' then a+2
            elif protSeq.[a-2]='C' then a-2
            elif protSeq.[a+3]='C' then a+3
            elif protSeq.[a-3]='C' then a-3
            else a
        printfn "seq_Length=%i, item at %i is %c" protSeq.Length cP protSeq.[cP]
        fillCysItem 0 org protID protSeq (cP) "" label]

let redoxPart =
    itemsA
    |> List.filter (fun x -> x<>[])
    |> List.map readA
    |> List.concat
    |> List.toArray
    |> Array.filter (fun i -> i.ProteinSeq.[i.CysPosition]='C')


let findNegative (protSeq: string) cysPP = ///// check the function, if the CysPosition is coreecly indexed!
    let cysPs = protSeq |> String.toCharArray |> Array.indexed |> Array.filter (fun (i,res) -> res='C') |> Array.map (fun (x,_) -> x)
    cysPP |> Set.ofArray |> Set.difference (cysPs |> Set.ofArray) |> Set.toArray

let createNegative (dbPart: CysItem []) (seqItem: CysItem)  =
    let cysPP = dbPart |> Array.filter (fun i -> i.ProtID = seqItem.ProtID ) |> Array.map (fun i -> i.CysPosition)
    let cysPN = findNegative seqItem.ProteinSeq cysPP
    cysPN
    |> Array.map (fun i -> fillCysItem 0 seqItem.Organism seqItem.ProtID seqItem.ProteinSeq i "" false)

let negativePool = 
    redoxPart
    |> Array.distinctBy (fun i -> i.ProtID)
    |> Array.map (createNegative redoxPart)
    |> Array.concat

let negativeRedoxPart =
    let r = new Random()
    Array.sampleWithOutReplacement r (negativePool) (redoxPart.Length)

let redoxData =
    [redoxPart;negativeRedoxPart]
    |> Array.concat
    |> Array.mapi enum


redoxData.[3073 ..] |> Array.filter (fun i -> i.ProteinSeq.[i.CysPosition]='C') |> Array.length
redoxData.Length

/// read Lemair data

let fileDirLem = @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\Lemaire_2017.fasta"

//reads from file to an array of FastaItems.
let sequencesLem = 
    fileDirLem
    |> fromFile (Seq.toArray >> String.fromCharArray)
    |> Seq.toArray
    |> Array.map (fun i -> 
            let h = i.Header |> String.split '|' |> Array.item 1 //|> String.tryParseIntDefault 0
            {i with Header=h})

let linesC1 = FileIO.readFile @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\Lemaire_2017.txt"

let lines =
    linesC1
    |> Seq.toList  
    |> List.map (fun x -> (x |> String.split '\t'))
let llines = 
    lines.[5 ..] 
    |> List.map (fun i -> [|i.[14];i.[6];i.[5]|]) 
    |> List.filter (fun i -> i.[1]<>"")
    |> List.mapFold (fun lastProtID i -> 
        if i.[0]="" then
            [|lastProtID;i.[1];i.[2]|], lastProtID
        else
            i, i.[0]
        ) ""
    |> fst

let itemsC =
    llines
    |> List.map (fun x -> 
        let org = "Chlamy"
        let protID = x.[0]
        let protSeq = (sequencesLem |> Array.find (fun i -> i.Header=protID)).Sequence
        let cysSeq = x.[2]
        let cysPs = 
            x.[1] 
            |> String.replace " ; " ";" 
            |> String.replace "C" "" 
            |> String.replace "\"" ""
            |> String.split ';'  
            |> Array.map (String.tryParseIntDefault 0)
        let label = true
        [for a in cysPs -> 
            let cP = 
                if a>=protSeq.Length-1  then adjustCysPosition protSeq cysSeq cysPs a
                elif protSeq.[a-1]='C' then a-1
                elif protSeq.[a]='C' then a
                else adjustCysPosition protSeq cysSeq cysPs a
            
            fillCysItem 0 org protID protSeq cP cysSeq label] )
    |> List.concat
    |> List.toArray
    |> Array.filter (fun i -> i.ProteinSeq.[i.CysPosition]='C')
    
let negativePoolC = 
    itemsC
    |> Array.distinctBy (fun i -> i.ProtID)
    |> Array.map (createNegative itemsC)
    |> Array.concat

let negativePartC =
    let r = new Random()
    Array.sampleWithOutReplacement r (negativePoolC) (itemsC.Length)

let chlamyData =
    [itemsC;negativePartC]
    |> Array.concat
itemsC.Length

/// Read Ara data

let fileDirA = @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\A.fasta"

//reads from file to an array of FastaItems.
let sequencesAra = 
    fileDirA
    |> fromFile (Seq.toArray >> String.fromCharArray)
    |> Seq.toArray
    |> Array.map (fun i -> 
            let h = i.Header |> String.split '|' |> Array.item 0 |> String.subString 0 9//|> String.tryParseIntDefault 0
            {i with Header=h})


let itemsAra = 
    [FileIO.readFile @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\Alix_Fares2011.txt" 
        |> Seq.toList  
        |> List.tail
        |> List.map (fun x -> 
            let i = (x |> String.split '\t')
            [|i.[0];i.[6];i.[2]|]) 
        ;
    FileIO.readFile @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\Alix_Liu2014.txt"  
        |> Seq.toList  
        |> List.tail 
        |> List.tail
        |> List.map (fun x -> 
            let i = (x |> String.split '\t')
            [|i.[0];i.[6];(i.[2] |> String.replace "*" "")|]) 
    ;]
    |> List.concat
    |> List.filter (fun i -> i.[1]<>"")
    |> List.mapFold (fun lastProtID i -> 
        if i.[0]="" then
            [|lastProtID;i.[1];i.[2]|], lastProtID
        else
            i, i.[0]
        ) ""
    |> fst
    |> List.map (fun x -> 
        let org = "Arabi"
        let protID = x.[0]
        let protSeq = (sequencesAra |> Array.find (fun i -> i.Header=protID)).Sequence
        let cysSeq = x.[2]
        let cysPs = 
            x.[1] 
            |> String.split ','  
            |> Array.map (String.tryParseIntDefault 0)
        let label = true
        [for a in cysPs -> 
            let cP = 
                if a>=protSeq.Length-1  then adjustCysPosition protSeq cysSeq cysPs a
                elif protSeq.[a-1]='C' then a-1
                elif protSeq.[a]='C' then a
                else adjustCysPosition protSeq cysSeq cysPs a
            
            fillCysItem 0 org protID protSeq cP cysSeq label] )
    |> List.concat
    |> List.toArray
    |> Array.filter (fun i -> i.ProteinSeq.[i.CysPosition]='C')
    

let negativePoolAra = 
    itemsAra
    |> Array.distinctBy (fun i -> i.ProtID)
    |> Array.map (createNegative itemsAra)
    |> Array.concat

let negativePartAra =
    let r = new Random()
    Array.sampleWithOutReplacement r (negativePoolAra) (itemsAra.Length)

let araData =
    [itemsAra;negativePartAra]
    |> Array.concat

itemsAra.Length

/// the whole data

let data = 
    [redoxPart; itemsC; itemsAra; negativePartAra; negativePartC; negativeRedoxPart; araData] 
    |> Array.concat 
    |> Array.distinctBy (fun x -> (x.ProteinSeq,x.CysPosition))
    |> Array.mapi enum

data.Length // 8702

data |> Array.filter (fun i -> i.Label=true) |> Array.length
data |> Array.filter (fun i -> i.Label=false) |> Array.length

data |> Array.filter (fun i -> (i.Organism |> String.contains "Chlamy") && i.Label=true) |> Array.length
data |> Array.filter (fun i -> (i.Organism |> String.contains "Arab") && i.Label=true) |> Array.length
data |> Array.filter (fun i -> not (i.Organism |> String.contains "Arab") && not (i.Organism |> String.contains "Chlamy")) |> Array.length

/////////////// FEATURES

let one_hot_encoding (sequence: string) =
    let aa = "ARNDCEQGHILKMFPSTWYV"
    let m = Array2D.zeroCreate sequence.Length 20
    for ii in [0 .. sequence.Length-1] do
        for ai in [0 .. aa.Length-1] do
            m.[ii,ai] <- if sequence.[ii]=aa.[ai] then 1 else 0 
    m

let getFlankingRegion fl (item: CysItem) =
    let startP = 
        if item.CysPosition-fl<0 then 0
        else item.CysPosition-fl
    let len =
        if (item.CysPosition+fl)>(item.ProteinSeq.Length-1) then (item.ProteinSeq.Length-item.CysPosition-1)
        else (fl + fl + 1)
    item.ProteinSeq |> String.subString startP len

let getFlankingRegionChars fl (item: CysItem) =
    let startP = 
        if item.CysPosition-fl<0 then 0
        else item.CysPosition-fl
    let len =
        if (item.CysPosition+fl)>(item.ProteinSeq.Length-1) then (item.ProteinSeq.Length-item.CysPosition-1)
        else (fl + fl + 1)
    item.ProteinSeq 
    |> String.subString startP len
    |> String.toCharArray |> Array.fold (fun acc i -> sprintf "%s %c" acc i) ""

let getFlankingRegionTabs fl (item: CysItem) =
    let emptyP = 
        if (item.CysPosition-fl<0) then -(item.CysPosition-fl) + 1
        else 0
    let startP = 
        if item.CysPosition-fl<0 then 0
        else item.CysPosition-fl
    let emptyAfter = 
        if (item.CysPosition+fl)>=(item.ProteinSeq.Length) then ((item.CysPosition + fl) - (item.ProteinSeq.Length - 1) + 1)
        else 0
    let len =
        if (item.CysPosition+fl)>(item.ProteinSeq.Length-1) then (item.ProteinSeq.Length-item.CysPosition+fl)-emptyP
        else (fl + fl + 1)-emptyP
    let ending = String.Join("\t", Array.create emptyAfter "")
    let beginning = String.Join("\t", Array.create emptyP "")
    (item.ProteinSeq 
    |> String.subString startP len
    |> String.toCharArray 
    |> Array.fold (fun acc i -> sprintf "%s\t%c" acc i) beginning)
    + ending

data.[5000] |> getFlankingRegion 2
data.[5000] |> getFlankingRegionChars 2
data.[5000] |> getFlankingRegionTabs 4
data.[35] |> getFlankingRegionTabs 4 |> String.split '\t'


let one_hot_encoding_tabbed (sequence_t: string) =
    let aa = "ARNDCEQGHILKMFPSTWYV"
    let sequence = (sequence_t |> String.split '\t').[1 ..]
    let m = Array2D.zeroCreate sequence.Length 20
    for ii in [0 .. sequence.Length-1] do
        for ai in [0 .. aa.Length-1] do
            m.[ii,ai] <- if sequence.[ii]=(string aa.[ai]) then 1 else 0 
    m
    |> Array2D.array2D_to_seq
    |> Seq.toArray

let itemToLine (item: CysItem) =
    let label = 
        if item.Label then 
            "|labels 1 0 "
        else
            "|labels 0 1 "
    let feat =
        item
        |> getFlankingRegionTabs 12
        |> one_hot_encoding_tabbed
        |> Array.fold (fun acc i -> sprintf "%s %i" acc i)
            "|features"
    label + feat

let dataLines = data |> Array.map itemToLine |> Array.filter (fun x -> x.Length=1021)
dataLines.Length // 8702 with any string length / 8543 with 1021 string length


let dataShaffled = dataLines |> Array.shuffleFisherYates

File.AppendAllLines(@"c:\Users\mikha\Downloads\redoxData_NN_test.txt", dataShaffled.[0 .. 99])  
File.AppendAllLines(@"c:\Users\mikha\Downloads\redoxData_NN_train.txt", dataShaffled.[100 .. ])  


dataShaffled.Length

data.[35] |> getFlankingRegionTabs 4 |> one_hot_encoding_tabbed 

let dataFlanked = data |> Array.map (fun i -> [(string i.Label) ; (getFlankingRegion 2 i)])
let dataFlankedSpaced = data |> Array.map (fun i -> [(string i.Label) ; (getFlankingRegionChars 2 i)])
let dataFlankedTabbed = data |> Array.map (fun i -> (string i.Label) + (getFlankingRegionTabs 12 i))

let redoxDataFlankedTabbed = redoxData |> Array.map (fun i -> (string i.Label) + (getFlankingRegionTabs 12 i))


let writeFileData (dataFlanked: string list []) path =
    File.AppendAllLines(path, dataFlanked |> Array.map (fun i ->  String.Join("\t", i)))

writeFileData dataFlanked @"c:\Users\mikha\Downloads\data.txt"
writeFileData dataFlankedSpaced @"c:\Users\mikha\Downloads\dataSpaced.txt"
File.AppendAllLines(@"c:\Users\mikha\Downloads\data25Tabbed.txt", dataFlankedTabbed)  
File.AppendAllLines(@"c:\Users\mikha\Downloads\redoxData25Tabbed.txt", redoxDataFlankedTabbed)  

/////////////////// ML

#r @"..\lib\FSharpML.dll"
#load @"c:\Users\mikha\source\repos\FSharpML\FSharpML.fsx"

open System;
open Microsoft.ML
open Microsoft.ML.Data;
open FSharpML
open FSharpML.EstimatorModel
open FSharpML.TransformerModel


/// Type representing the Message to run analysis on.
[<CLIMutable>] 
type SpamInput = 
    { 
        [<LoadColumn(0)>] LabelText : string
        [<LoadColumn(1)>] Message : string 
    }

//Create the MLContext to share across components for deterministic results
let mlContext = MLContext(seed = Nullable 1) // Seed set to any number so you
                                             // have a deterministic environment

// STEP 1: Common data loading configuration   
let fullData = 
    @"c:\Users\mikha\Downloads\data.txt"
    |> DataModel.fromTextFileWith<SpamInput> mlContext '\t' false 

let trainingData, testingData = 
    fullData
    |> DataModel.BinaryClassification.trainTestSplit 0.1

//STEP 2: Process data, create and train the model 
let model = 
    EstimatorModel.create mlContext
    // Process data transformations in pipeline
    |> EstimatorModel.appendBy (fun mlc -> mlc.Transforms.Conversion.ValueMap(["ham"; "spam"],[false; true],[| struct (DefaultColumnNames.Label, "LabelText") |]))
    |> EstimatorModel.appendBy (fun mlc -> mlc.Transforms.Text.FeaturizeText(DefaultColumnNames.Features, "Message"))
    |> EstimatorModel.appendCacheCheckpoint
    // Create the model
    |> EstimatorModel.appendBy (fun mlc -> mlc.BinaryClassification.Trainers.StochasticDualCoordinateAscent(DefaultColumnNames.Label, DefaultColumnNames.Features))
    // Train the model
    |> EstimatorModel.fit trainingData.Dataview

// STEP3: Run the prediciton on the test data
let predictions =
    model
    |> TransformerModel.transform testingData.Dataview

// STEP4: Evaluate accuracy of the model
let metrics = 
    model
    |> Evaluation.BinaryClassification.evaluate testingData.Dataview

metrics.Accuracy
metrics.Auc
metrics.NegativeRecall
metrics.PositiveRecall

//// Apply used in the article evaluation metrics: ACC, SN, SP, MCC, AUC




//// STEP5: Create prediction engine function related to the loaded trained model
//let predict = 
//    TransformerModel.createPredictionEngine<_,SpamInput,SpamInput> model

//// Score
//let prediction = predict sampleStatement

//// Test a few examples
//[
//    "That's a great idea. It should work."
//    "free medicine winner! congratulations"
//    "Yes we should meet over the weekend!"
//    "you win pills and free entry vouchers"
//] 
//|> List.iter (classify predictor)
