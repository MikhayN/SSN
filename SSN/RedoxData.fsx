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

open FSharpAux
open FSharpAux.IO
open BioFSharp.BioID

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

//Functions

type CysItem =
    {
    CysID: int;
    Organism: string;
    ProtID: string;
    ProteinSeq: string;
    CysPosition: int;
    CysSeq: string;
    Label: bool;
    IsPlant: bool;
    }

let fillCysItem cysID org protID protSeq cysP cysSeq label plant = {CysID=cysID; Organism=org; ProtID=protID; ProteinSeq=protSeq; CysPosition=cysP; CysSeq=cysSeq; Label=label; IsPlant=plant}

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

let one_hot_encoding (sequence: string) =
    let aa = "ARNDCEQGHILKMFPSTWYV"
    let m = Array2D.zeroCreate sequence.Length 20
    for ii in [0 .. sequence.Length-1] do
        for ai in [0 .. aa.Length-1] do
            m.[ii,ai] <- if sequence.[ii]=aa.[ai] then 1 else 0 
    m
    |> Array2D.array2D_to_seq
    |> Seq.toArray

let getFlankingRegion fl (item: CysItem) =
    let emptyP = 
        if (item.CysPosition-fl<0) then -(item.CysPosition-fl)
        else 0
    let startP = 
        if item.CysPosition-fl<0 then 0
        else item.CysPosition-fl
    let emptyAfter = 
        if (item.CysPosition+fl)>=(item.ProteinSeq.Length) then ((item.CysPosition + fl) - (item.ProteinSeq.Length - 1) )
        else 0
    let len =
        if (item.CysPosition+fl)>(item.ProteinSeq.Length-1) then (item.ProteinSeq.Length-item.CysPosition+fl)-emptyP
        else (fl + fl + 1)-emptyP
    let ending = String.fromCharArray(Array.create emptyAfter '_')
    let beginning = String.fromCharArray(Array.create emptyP '_')
    beginning + (item.ProteinSeq |> String.subString startP len) + ending

let getFlankingRegionTabs fl (item: CysItem) =
    let emptyP = 
        if (item.CysPosition-fl<0) then -(item.CysPosition-fl)
        else 0
    let startP = 
        if item.CysPosition-fl<0 then 0
        else item.CysPosition-fl
    let emptyAfter = 
        if (item.CysPosition+fl)>=(item.ProteinSeq.Length) then ((item.CysPosition + fl) - (item.ProteinSeq.Length - 1) )
        else 0
    let len =
        if (item.CysPosition+fl)>(item.ProteinSeq.Length-1) then (item.ProteinSeq.Length-item.CysPosition+fl)-emptyP
        else (fl + fl + 1)-emptyP
    let ending = String.fromCharArray(Array.create emptyAfter '_')
    let beginning = String.fromCharArray(Array.create emptyP '_')
    let charA =
        beginning + (item.ProteinSeq |> String.subString startP len) + ending
        |> String.toCharArray 
    String.Join("\t", charA)
    

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

///////////// prepare PlantNamesSet

let plantSet =
    ["9POAL"; "9RHOD"; "9ROSI"; "ARATH"; "CARPA"; 
    "CHLIN"; "CHLRE"; "HORVD"; "JACME"; "MAIZE";
    "MEDSA"; "ORYSA"; "PEA"; "PRUAR"; "SOYBN";
    "SPIOL"; "WHEAT"; "DATST"  ]

/// Generate negative samples from the set of positive
let genNeg posDB = 

    let findNegative (protSeq: string) cysPP = ///// check the function, if the CysPosition is coreecly indexed!
        let cysPs = protSeq |> String.toCharArray |> Array.indexed |> Array.filter (fun (i,res) -> res='C') |> Array.map (fun (x,_) -> x)
        cysPP |> Set.ofArray |> Set.difference (cysPs |> Set.ofArray) |> Set.toArray

    let createNegative (dbPart: CysItem []) (seqItem: CysItem)  =
        let cysPP = dbPart |> Array.filter (fun i -> i.ProtID = seqItem.ProtID ) |> Array.map (fun i -> i.CysPosition)
        let cysPN = findNegative seqItem.ProteinSeq cysPP
        cysPN
        |> Array.map (fun i -> fillCysItem 0 seqItem.Organism seqItem.ProtID seqItem.ProteinSeq i "" false seqItem.IsPlant)

    let negativePool = 
        posDB
        |> Array.distinctBy (fun i -> i.ProtID)
        |> Array.map (createNegative posDB)
        |> Array.concat

    let r = new Random(0)
    Array.sampleWithOutReplacement r (negativePool) (posDB.Length)

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
    let plant = plantSet |> List.contains (protID |> String.split '_').[1]
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
        //printfn "seq_Length=%i, item at %i is %c" protSeq.Length cP protSeq.[cP]
        fillCysItem 0 org protID protSeq (cP) "" label plant]

let redoxPart =
    itemsA
    |> List.filter (fun x -> x<>[])
    |> List.map readA
    |> List.concat
    |> List.toArray
    |> Array.filter (fun i -> i.ProteinSeq.[i.CysPosition]='C')
    
let negativeRedoxPart = genNeg redoxPart
    
let redoxData =
    [redoxPart;negativeRedoxPart]
    |> Array.concat
    |> Array.mapi enum

redoxPart.Length
redoxPart |> Array.filter (fun i -> i.IsPlant) |> Array.length
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
        let org = " Chlamydomonas reinhardtii"
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
            
            fillCysItem 0 org protID protSeq cP cysSeq label true] )
    |> List.concat
    |> List.toArray
    |> Array.filter (fun i -> i.ProteinSeq.[i.CysPosition]='C')
    
let negativePartC = genNeg itemsC

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
        let org = " Arabidopsis thaliana (Mouse-ear cress)"
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
            
            fillCysItem 0 org protID protSeq cP cysSeq label true] )
    |> List.concat
    |> List.toArray
    |> Array.filter (fun i -> i.ProteinSeq.[i.CysPosition]='C')
    
let negativePartAra = genNeg itemsAra

let araData =
    [itemsAra;negativePartAra]
    |> Array.concat

itemsAra.Length

/// the whole data

let data = 
    [redoxPart; itemsC; itemsAra; negativePartAra; negativePartC; negativeRedoxPart] 
    |> Array.concat 
    |> Array.distinctBy (fun x -> (x.ProteinSeq,x.CysPosition) ) // avoid repetition of samples
    |> Array.mapi enum

data.Length // 8702

data |> Array.filter (fun i -> i.Label=true) |> Array.length
data |> Array.filter (fun i -> i.Label=false) |> Array.length

data |> Array.filter (fun i -> (i.Organism |> String.contains "Chlamy") && i.Label=true) |> Array.length
data |> Array.filter (fun i -> (i.Organism |> String.contains "Arab") && i.Label=true) |> Array.length
data |> Array.filter (fun i -> not (i.Organism |> String.contains "Arab") && not (i.Organism |> String.contains "Chlamy")) |> Array.length

data 
|> Array.distinctBy (fun i -> i.Organism) 
|> Array.map (fun i -> i.ProtID |> String.split '_' |> (fun arr ->  if arr.Length>1 then arr.[1] else "") , i.Organism)
|> Array.sortBy (fun (x,y) -> x)
|> Array.iter (fun x -> printfn "%A" x)

/// Evaluation data: RSC and BALOSCTdb

/// RSC dataset
let dataRSC =
    FileIO.readFile @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\dataFiles\RSC758.txt" 
    |> Seq.toList  
    |> List.tail
    |> List.map (fun x -> 
        let i = (x |> String.split '\t')
            
        i.[0 .. 3]) 
    |> List.map (fun x -> 

        let org = (x.[0] |> String.split '_').[1]
        let protID = x.[0]
        let protSeq = x.[3] //""
        let cysSeq = x.[3]
        let cysP = 10 //x.[1] |> String.tryParseIntDefault 0
        let label = x.[2]="1"
        let plant = plantSet |> List.contains org
        
        fillCysItem 0 org protID protSeq cysP cysSeq label plant)

    |> List.toArray
    |> Array.filter (fun i -> i.CysSeq.[10]='C')

dataRSC |> Array.filter (fun i -> i.IsPlant) |> Array.length
dataRSC.Length
let itemsRSC_filtered =
    dataRSC
    |> Array.filter (fun x -> not (redoxData |> Array.exists (fun i -> i.ProtID=x.ProtID)))
itemsRSC_filtered.Length

///
let mappingFile =
    FileIO.readFile @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\pdbtosp.txt" 
    |> Seq.toList  
    |> fun i -> i.[24 ..]
    |> List.filter (fun i -> i.Length>1)
    |> List.map (fun x -> (x |> String.subString 0 4), (x |> String.subString  28 12 |> String.filter (fun c -> c<>' ')) )

/// Balanced OSCTdb
let dataBALOSCT =
    FileIO.readFile @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\dataFiles\BALOSCTdb.txt" 
    |> Seq.toList  
    |> List.tail
    |> List.map (fun x -> 
        let i = (x |> String.split '\t')
            
        i.[0 .. 4]) 
    |> List.map (fun x -> 

        //let org = (x.[0] |> String.split '_').[1]
        let protID' = (x.[0] |> String.split '_').[0]
        let protID = 
            mappingFile 
            |> List.tryFind (fun (k,v) -> k=protID') 
            |> (fun xx -> if xx.IsSome then xx.Value |> snd else x.[4])
        let org = (protID |> String.split '_').[1]
        let protSeq = x.[3] // ""
        let cysSeq = x.[3]
        let cysP = 10 // x.[1] |> String.tryParseIntDefault 0
        let label = x.[2]="1"
        let plant = plantSet |> List.contains org
        
        fillCysItem 0 org protID protSeq cysP cysSeq label plant)

    |> List.toArray
    |> Array.filter (fun i -> i.CysSeq.[10]='C')

dataBALOSCT |> Array.filter (fun i -> i.IsPlant) |> Array.length
dataBALOSCT.Length
let dataBALOSCT_filtered =
    dataBALOSCT
    |> Array.filter (fun x -> not (redoxData |> Array.exists (fun i -> i.ProtID=x.ProtID)))

dataBALOSCT_filtered |> Array.map (fun i -> i.Organism) |> Array.distinct

/// Arabidopsis Dataset from Puyaubert article
let itemsAraTest =
    FileIO.readFile @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\Ara_TEST.txt" 
    |> Seq.toList  
    |> List.tail
    |> List.map (fun x -> 
        let i = (x |> String.split '\t')
        [|i.[1] |> String.map (fun c -> c |> Char.ToUpper)
        ;i.[2]
        ;i.[5] |> String.filter (fun c -> c |> Char.IsLetter) |])

    |> List.map (fun x -> 
        let org = "Arabi"
        let protID = x.[0]
        let protSeq = (sequencesAra |> Array.find (fun i -> i.Header=protID)).Sequence
        let cysSeq = x.[2]
        let label = true
        let a = int x.[1] 
        let cP = 
                if a>=protSeq.Length-1  then adjustCysPosition protSeq cysSeq [|a|] a
                elif protSeq.[a-1]='C' then a-1
                elif protSeq.[a]='C' then a
                else adjustCysPosition protSeq cysSeq [|a|] a
            
        fillCysItem 0 org protID protSeq cP cysSeq label true)

    |> List.toArray
    |> Array.filter (fun i -> i.ProteinSeq.[i.CysPosition]='C')

itemsAraTest.Length // 61 out of 62

let itemsAraTest_filtered =
    itemsAraTest
    |> Array.filter (fun x -> not (data |> Array.exists (fun i -> i.ProteinSeq=x.ProteinSeq && i.CysPosition=x.CysPosition )))

itemsAraTest_filtered.Length // 35 or 37

let linesPlants = itemsAraTest_filtered |> Array.map (fun x -> [sprintf ">%s" x.ProtID; x.ProteinSeq; ""]) |> List.concat

FileIO.writeToFile false @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\PlantTest.fa" linesPlants

/// All PLANTS DB

let plantDB =
    let fromRedoxDB = redoxPart |> Array.filter (fun i -> i.IsPlant)
    let positive = 
        [|fromRedoxDB; itemsC; itemsAra;|] 
        |> Array.concat
        |> Array.distinctBy (fun x -> (x.ProteinSeq,x.CysPosition))
    let negative = genNeg positive
    [positive; negative]
    |> Array.concat
    |> Array.mapi enum

plantDB.Length 

///////////////////////////////// Final TRAINING & TESTING sets

/// Dataset A : plants only
let setA_pos =
    let fromRedoxDB = redoxPart |> Array.filter (fun i -> i.IsPlant)
    let positive = 
        [|fromRedoxDB; itemsC; itemsAra; itemsAraTest_filtered|] 
        |> Array.concat
        |> Array.distinctBy (fun x -> (x.ProteinSeq,x.CysPosition))
    positive
    |> Array.mapi enum

/// Dataset B : everything not plants
let setB_pos =
    redoxPart 
    |> Array.filter (fun i -> not i.IsPlant)
    |> Array.mapi enum

let (setT_test, setA_rest) =
    setA_pos 
    |> Array.shuffleFisherYates
    |> Array.splitAt 100
    
let setC_train =
    [|setB_pos; genNeg setB_pos; setA_rest; genNeg setA_rest|]
    |> Array.concat
    |> Array.shuffleFisherYates
         

/////////////// FEATURES

//data.[35] |> getFlankingRegionTabs 4 |> one_hot_encoding_tabbed 

//let dataFlanked = data |> Array.map (fun i -> [(string i.Label) ; (getFlankingRegion 4 i)])

let dataFlankedTabbed = 
    setC_train
    |> Array.map (fun i -> 
        (string i.Label) + "\t" + (getFlankingRegionTabs 4 i) 
        |> String.replace "_" "")
        
File.AppendAllLines(@"c:\Users\mikha\Downloads\TRAIN_9.txt", dataFlankedTabbed)

//// for testing Redox predictor

let linesPlants' = itemsAraTest_filtered |> Array.map (fun x -> [sprintf ">%s" x.ProtID; getFlankingRegion 4 x; ""]) |> List.concat

FileIO.writeToFile false @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\PlantTestCysSeq.fa" linesPlants'

// to NN

let itemToLine_forNN fl (item: CysItem) =
    let label = 
        if item.Label then 
            "|labels 1 0 "
        else
            "|labels 0 1 "
    let feat =
        item
        |> getFlankingRegion fl
        |> one_hot_encoding
        |> Array.fold (fun acc i -> sprintf "%s %i" acc i)
            "|features"
    label + feat

let dataLines = 
    dataBALOSCT 
    |> Array.map (itemToLine_forNN 4) 

let dataShaffled = dataLines |> Array.shuffleFisherYates

File.AppendAllLines(@"c:\Users\mikha\Downloads\dataBALOSCT_NN_test.txt", dataShaffled.[0 ..])  
File.AppendAllLines(@"c:\Users\mikha\Downloads\dataRSC_NN_train.txt", dataShaffled.[100 .. ])  

/////////////////// ML -> look into FSharpML.sln
